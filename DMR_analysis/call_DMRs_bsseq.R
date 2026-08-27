#!/usr/bin/env Rscript
# =============================================================================
# 02_call_dmrs.R
#
# DMR calling for cfDNA Nanopore methylation (modkit pileup input),
# replicating the analytical framework of the cfMeDiP BRCA1/2 study
# (Differential methylation of cfDNA via cfMeDiP) with data-type-appropriate
# modifications:
#   - Direct 5mC calls from modkit (not MeDIP density) -> real methylation
#     fractions, so MEDIPS / relH / GoGe are dropped.
#   - 5hmC ("h") rows excluded; DMRs called on 5mC ("m") only.
#   - Coverage-aware beta-binomial test with shrinkage (DSS), which substitutes
#     for the density-pooling stability MEDIPS relied on.
#   - HYBRID output: DSS callDMR data-driven regions as primary biology,
#     PLUS a 300-bp tiled table for comparability with the origin study.
#   - Study thresholds preserved: |log2FC| >= 2 (on group-mean fractions)
#     AND p < 0.01.
#   - User coverage protocol enforced as an explicit post-hoc filter:
#     a CpG is "covered" if Nvalid >= 1; a sample is kept in a region only if
#     it has >= MIN_COVERED_CPGS covered CpGs there; a region is kept only if
#     >= MIN_VALID_SAMPLES_PER_GROUP samples pass in EACH group. This step is
#     what defends against the DSS sparse-smoothing artifact where CpGs covered
#     in only one group receive spuriously low p-values.
#
# Usage:
#   Rscript 02_call_dmrs.R \
#       --manifest samples.tsv \      # columns: sample_id  group  path
#       --outdir   results/dmr \
#       [--min-covered-cpgs 3] [--min-valid-per-group 4] \
#       [--tile 300] [--p 0.01] [--log2fc 2] [--smoothing-span 500]
#
# samples.tsv (tab-separated, header required):
#   sample_id   group   path
#   BRCA2_01    carrier /path/to/BRCA2_01.bed
#   WT_01       control /path/to/WT_01.bed
#   ...
#   group must contain exactly two levels; "carrier" is treated as group1 so
#   that positive log2FC / hypermethylation means carrier > control.
# =============================================================================

suppressMessages({
  library(optparse)
  library(DSS)
  library(bsseq)
  library(data.table)
  library(GenomicRanges)
})

# ------------------------- arguments -----------------------------------------
opt_list <- list(
  make_option("--manifest", type = "character", default = NULL,
              help = "Optional TSV (sample_id, group, path). If omitted, the directory is scanned."),
  make_option("--dir", type = "character",
              default = "/mnt/d/BRCA2-misc-files/BRCA2-ctrl-mod",
              help = "Directory of .bed files to scan when no manifest is given [%default]"),
  make_option("--group1-prefix", type = "character", default = "BRCA2",
              dest = "group1_prefix",
              help = "Filename prefix marking group1 (carrier) [%default]"),
  make_option("--group2-prefix", type = "character", default = "ctrl",
              dest = "group2_prefix",
              help = "Filename prefix marking group2 (control) [%default]"),
  make_option("--outdir", type = "character", default = "results/dmr"),
  make_option("--min-covered-cpgs", type = "integer", default = 3L,
              dest = "min_covered_cpgs",
              help = "Min covered CpGs (Nvalid>=1) per sample within a region [%default]"),
  make_option("--min-valid-per-group", type = "integer", default = 4L,
              dest = "min_valid_per_group",
              help = "Min passing samples required in EACH group per region [%default]"),
  make_option("--tile", type = "integer", default = 300L,
              help = "Tile width in bp for the comparability layer [%default]"),
  make_option("--p", type = "double", default = 0.01,
              help = "p-value threshold (origin study) [%default]"),
  make_option("--log2fc", type = "double", default = 2.0,
              help = "abs(log2FC) threshold on group-mean fractions [%default]"),
  make_option("--smoothing-span", type = "integer", default = 500L,
              dest = "smoothing_span",
              help = "DSS smoothing span in bp [%default]"),
  make_option("--min-cg", type = "integer", default = 3L, dest = "min_cg",
              help = "callDMR minCG (CpGs per DMR) [%default]"),
  make_option("--group1", type = "character", default = "carrier",
              help = "Label for group1 (positive FC = this group higher) [%default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))

# ------------------------- build sample table --------------------------------
# Two ways to define samples:
#   (1) --manifest file.tsv  (columns sample_id, group, path), OR
#   (2) scan --dir for *.bed and assign group by filename prefix.
# Group labels are normalised to opt$group1 ("carrier") and "control" so the
# rest of the script and the FC direction stay consistent.
group2_label <- "control"

if (!is.null(opt$manifest)) {
  man <- fread(opt$manifest)
  stopifnot(all(c("sample_id", "group", "path") %in% names(man)))
} else {
  if (!dir.exists(opt$dir))
    stop(sprintf("Directory not found: %s\n(On WSL, D:\\... is /mnt/d/...)", opt$dir))
  beds <- list.files(opt$dir, pattern = "\\.bed$", full.names = TRUE)
  if (length(beds) == 0)
    stop(sprintf("No .bed files found in %s", opt$dir))
  bn <- basename(beds)
  grp <- ifelse(startsWith(bn, opt$group1_prefix), opt$group1,
         ifelse(startsWith(bn, opt$group2_prefix), group2_label, NA_character_))
  if (any(is.na(grp))) {
    bad <- bn[is.na(grp)]
    stop(sprintf("%d file(s) matched neither prefix '%s' nor '%s': %s",
                 length(bad), opt$group1_prefix, opt$group2_prefix,
                 paste(head(bad, 5), collapse = ", ")))
  }
  man <- data.table(
    sample_id = sub("\\.bed$", "", bn),
    group     = grp,
    path      = beds
  )
  log_msg("Scanned %s: %d files (%d %s, %d %s)", opt$dir, nrow(man),
          sum(grp == opt$group1), opt$group1,
          sum(grp == group2_label), group2_label)
}

groups <- unique(man$group)
if (length(groups) != 2)
  stop(sprintf("Need exactly 2 groups; found: %s", paste(groups, collapse = ", ")))
if (!opt$group1 %in% groups)
  stop(sprintf("group1 label '%s' not present after grouping", opt$group1))
group2 <- setdiff(groups, opt$group1)
log_msg("group1 (positive FC = higher here): %s", opt$group1)
log_msg("group2: %s", group2)
log_msg("Loading %d samples", nrow(man))

# ------------------------- read one modkit bedMethyl -------------------------
# modkit pileup bedMethyl columns (0-based description, 1-based index here):
#   1 chrom  2 start  3 end  4 mod_code  5 score  6 strand
#   7 tstart 8 tend   9 color 10 Nvalid  11 percent_mod  12 Nmod
#   13 Ncanonical 14 Nother 15 Ndelete 16 Nfail 17 Ndiff 18 Nnocall
# We use: chrom(1), start(2), mod_code(4), Nvalid(10), Nmod(12).
# "Covered" = Nvalid >= 1. We keep ONLY mod_code == "m" (5mC); "h" (5hmC) dropped.
read_modkit <- function(path, sample_id) {
  # Filter to mod_code == "m" (5mC) at the OS level with awk, so 5hmC ("h")
  # rows are never allocated in R. Column 4 is mod_code in modkit bedMethyl.
  # Fall back to plain fread if awk is unavailable (e.g. non-WSL).
  awk_ok <- nzchar(Sys.which("awk"))
  if (awk_ok) {
    dt <- fread(cmd = sprintf("awk -F'\\t' '$4==\"m\"' %s", shQuote(path)),
                header = FALSE, sep = "\t",
                select = c(1, 2, 10, 12),
                col.names = c("chr", "pos", "N", "X"))
    n_kept_5mc <- nrow(dt)
  } else {
    dt <- fread(path, header = FALSE, sep = "\t",
                select = c(1, 2, 4, 10, 12),
                col.names = c("chr", "pos", "mod", "N", "X"))
    dt <- dt[mod == "m"]
    dt[, mod := NULL]
    n_kept_5mc <- nrow(dt)
  }
  dt <- dt[N >= 1]                          # "covered" = Nvalid >= 1
  dt[, pos := pos + 1L]                     # bedMethyl 0-based -> DSS 1-based
  log_msg("  %s: %d 5mC rows, %d covered (Nvalid>=1)", sample_id,
          n_kept_5mc, nrow(dt))
  # DSS requires: chr, pos, N (coverage), X (methylated count)
  dt[, .(chr, pos, N, X)]
}

# ------------------------- split every sample to per-chromosome files --------
# ROOT-CAUSE FIX: never hold genome-wide data for all samples in RAM. Earlier
# versions kept all 48 sample tables resident (bs_list), which at this coverage
# (some samples 14-18M CpGs) is a multi-GB floor that pushed chr10 into swap and
# got OOM-killed. Instead we read each sample ONCE, split it by chromosome, and
# write per-(sample,chromosome) RDS files to a scratch dir. The loop then loads
# only the current chromosome's slice across samples. Peak RAM = one chromosome.
canonical <- paste0("chr", c(1:22, "X", "Y"))
scratch <- file.path(opt$outdir, "scratch_bychrom")
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)

# Split each sample to per-chromosome RDS files. A per-sample ".done" marker
# lets a re-run skip samples already split, so resume doesn't pay the full read
# cost again. Delete the scratch dir to force a clean re-split.
for (i in seq_len(nrow(man))) {
  done_marker <- file.path(scratch, sprintf("%s.done", man$sample_id[i]))
  if (file.exists(done_marker)) {
    log_msg("  %s already split, skipping", man$sample_id[i]); next
  }
  d <- read_modkit(man$path[i], man$sample_id[i])
  d <- d[chr %in% canonical]                     # drop non-canonical contigs early
  for (ch in unique(d$chr)) {
    saveRDS(d[chr == ch, .(chr, pos, N, X)],
            file.path(scratch, sprintf("%s__%s.rds", man$sample_id[i], ch)))
  }
  file.create(done_marker)                       # mark this sample complete
  rm(d); invisible(gc())                          # free this sample before the next
}
# Recover the chromosome set from whatever RDS files exist on disk (works
# whether we just split or resumed from a prior split).
existing <- list.files(scratch, pattern = "__chr.*\\.rds$")
seen_chroms <- unique(sub("^.*__(chr[^.]+)\\.rds$", "\\1", existing))
keep_chroms <- canonical[canonical %in% seen_chroms]
log_msg("Scratch ready in %s", scratch)
log_msg("Processing %d chromosomes: %s", length(keep_chroms),
        paste(keep_chroms, collapse = ", "))

g1_samples <- man$sample_id[man$group == opt$group1]
g2_samples <- man$sample_id[man$group == group2]

# ============================================================================
# PER-CHROMOSOME PROCESSING LOOP
# For each chromosome: subset every sample to that chromosome, build a small
# BSseq object, run DMLtest + callDMR, apply the coverage post-filter, and
# build the tile layer. Accumulate results across chromosomes, then write once.
# Peak memory is bounded by the largest single chromosome (~chr1), not the
# whole-genome union that OOM-killed the all-at-once build.
# ============================================================================

safe_log2fc <- function(f1, f2, eps = 1e-3) {
  f1 <- min(max(f1, eps), 1 - eps)
  f2 <- min(max(f2, eps), 1 - eps)
  log2(f1 / f2)
}

# Accumulators across chromosomes
dmr_candidates_all <- list()   # per-chrom data.tables of ALL candidate DMRs
tile_all           <- list()   # per-chrom data.tables of ALL tiles

for (chrom in keep_chroms) {
  log_msg("=== chromosome %s ===", chrom)

  # ---- resume: if this chromosome's checkpoint already exists, load it and
  #      skip recomputation. Lets a re-run continue from where a crash stopped.
  #      Delete the scratch dir (or pass --fresh handling) to force a clean run.
  ckpt_dmr  <- file.path(scratch, sprintf("dmrs_%s.tsv", chrom))
  ckpt_tile <- file.path(scratch, sprintf("tiles_%s.tsv.gz", chrom))
  if (file.exists(ckpt_dmr)) {
    log_msg("  checkpoint found for %s, loading and skipping recompute", chrom)
    dd <- fread(ckpt_dmr)
    if (nrow(dd) > 0) dmr_candidates_all[[chrom]] <- dd
    if (file.exists(ckpt_tile)) {
      tt <- fread(ckpt_tile)
      if (nrow(tt) > 0) tile_all[[chrom]] <- tt
    }
    next
  }

  # ---- load this chromosome's slice for every sample from scratch ----
  chr_files <- file.path(scratch, sprintf("%s__%s.rds", man$sample_id, chrom))
  have <- file.exists(chr_files)
  chr_list <- lapply(chr_files[have], readRDS)
  names(chr_list) <- man$sample_id[have]
  n_cpg_here <- sum(vapply(chr_list, nrow, integer(1)))
  if (n_cpg_here == 0 || length(chr_list) < 2) {
    log_msg("  no usable CpGs on %s, skipping", chrom); next }
  # samples present for THIS chromosome (some low-cov samples may lack a chrom)
  chr_samples <- names(chr_list)
  g1_here <- intersect(g1_samples, chr_samples)
  g2_here <- intersect(g2_samples, chr_samples)

  # guard: need at least 2 samples in EACH group present on this chromosome,
  # otherwise a two-group test is not defined here.
  if (length(g1_here) < 2 || length(g2_here) < 2) {
    log_msg("  %s: insufficient samples per group (g1=%d, g2=%d), skipping",
            chrom, length(g1_here), length(g2_here))
    rm(chr_list); invisible(gc()); next
  }

  # ---- build BSseq for this chromosome only ----
  BSobj <- makeBSseqData(chr_list, chr_samples)
  rm(chr_list); invisible(gc())

  # ---- DMLtest (smoothed) ----
  # Smoothing ON: neighboring CpGs act as pseudo-replicates for power at low
  # coverage. Its cost (low p-values at CpGs covered in only one group) is
  # neutralized by the coverage post-filter below, not by disabling smoothing.
  dml <- tryCatch(
    DMLtest(BSobj, group1 = g1_here, group2 = g2_here,
            smoothing = TRUE, smoothing.span = opt$smoothing_span),
    error = function(e) { log_msg("  DMLtest failed on %s: %s", chrom,
                                  conditionMessage(e)); NULL })
  if (is.null(dml)) { rm(BSobj); invisible(gc()); next }

  dml_pval <- dml$pval                       # indexed by CpG row, same order as loci

  # ---- data-driven DMRs ----
  dmrs <- callDMR(dml, p.threshold = opt$p, minCG = opt$min_cg, delta = 0)
  rm(dml); invisible(gc())
  if (is.null(dmrs)) dmrs <- data.frame()
  log_msg("  callDMR: %d candidate regions", nrow(dmrs))

  # ---- coverage matrices for this chromosome ----
  cov_mat  <- getCoverage(BSobj, type = "Cov")
  meth_mat <- getCoverage(BSobj, type = "M")
  loci     <- granges(BSobj)
  sample_order <- sampleNames(BSobj)
  rm(BSobj); invisible(gc())
  g1_idx <- which(sample_order %in% g1_here)
  g2_idx <- which(sample_order %in% g2_here)

  # region_stats closes over this chromosome's matrices
  region_stats <- function(rows) {
    covg1 <- cov_mat[rows, g1_idx, drop = FALSE]
    covg2 <- cov_mat[rows, g2_idx, drop = FALSE]
    covered_g1 <- colSums(covg1 >= 1)
    covered_g2 <- colSums(covg2 >= 1)
    pass_g1 <- sum(covered_g1 >= opt$min_covered_cpgs)
    pass_g2 <- sum(covered_g2 >= opt$min_covered_cpgs)
    keep_g1 <- g1_idx[covered_g1 >= opt$min_covered_cpgs]
    keep_g2 <- g2_idx[covered_g2 >= opt$min_covered_cpgs]
    frac_g1 <- if (length(keep_g1))
      sum(meth_mat[rows, keep_g1]) / sum(cov_mat[rows, keep_g1]) else NA_real_
    frac_g2 <- if (length(keep_g2))
      sum(meth_mat[rows, keep_g2]) / sum(cov_mat[rows, keep_g2]) else NA_real_
    list(pass_g1 = pass_g1, pass_g2 = pass_g2,
         frac_g1 = frac_g1, frac_g2 = frac_g2)
  }

  # ---- apply post-hoc coverage + log2FC filter to candidate DMRs ----
  # FIX: direction and log2FC now come from DSS's own meanMethy1/meanMethy2,
  # which are coverage-weighted and smoothing-aware and therefore robust. The
  # earlier version derived them from pooled Nmod/Nvalid fractions, which at low
  # absolute methylation collapse toward 0 and, once eps-floored, produce a
  # large POSITIVE log2FC regardless of true direction -- flipping hypo regions
  # to "hyper". Pooled fractions are retained as informational columns only.
  if (nrow(dmrs) > 0) {
    dmr_gr <- GRanges(dmrs$chr, IRanges(dmrs$start, dmrs$end))
    hits <- findOverlaps(loci, dmr_gr)
    rows_by_dmr <- split(queryHits(hits), subjectHits(hits))
    keep <- logical(nrow(dmrs))
    frac1 <- frac2 <- l2fc <- rep(NA_real_, nrow(dmrs))
    p1 <- p2 <- rep(0L, nrow(dmrs))
    for (j in seq_len(nrow(dmrs))) {
      rows <- rows_by_dmr[[as.character(j)]]
      if (is.null(rows)) next
      st <- region_stats(rows)
      p1[j] <- st$pass_g1; p2[j] <- st$pass_g2
      frac1[j] <- st$frac_g1; frac2[j] <- st$frac_g2
      # log2FC from DSS's robust group means (meanMethy1 = group1 = carrier).
      # These are bounded (0,1); eps-floor only guards exact 0/1 endpoints.
      l2fc[j] <- safe_log2fc(dmrs$meanMethy1[j], dmrs$meanMethy2[j])
      coverage_ok <- st$pass_g1 >= opt$min_valid_per_group &&
                     st$pass_g2 >= opt$min_valid_per_group
      keep[j] <- coverage_ok && abs(l2fc[j]) >= opt$log2fc
    }
    dmrs$log2FC <- l2fc
    # direction from DSS diff.Methy (meanMethy1 - meanMethy2): positive = group1
    # (carrier) higher = hyper. Using diff.Methy keeps direction consistent with
    # DSS's own statistic by construction.
    dmrs$direction <- ifelse(is.na(dmrs$diff.Methy), NA,
                      ifelse(dmrs$diff.Methy > 0, "hyper", "hypo"))
    dmrs$pooled_frac_group1 <- frac1   # informational only (noisy at low cov)
    dmrs$pooled_frac_group2 <- frac2   # informational only
    dmrs$n_pass_group1 <- p1
    dmrs$n_pass_group2 <- p2
    dmrs$passes_filter <- keep
    dmr_candidates_all[[chrom]] <- as.data.table(dmrs)
    # checkpoint: write this chromosome's candidates immediately, so a later
    # crash never loses completed chromosomes.
    fwrite(as.data.table(dmrs),
           file.path(scratch, sprintf("dmrs_%s.tsv", chrom)), sep = "\t")
  }

  # ---- tile layer for this chromosome ----
  max_pos <- max(start(loci))
  starts <- seq(1, max_pos + opt$tile, by = opt$tile)
  tiles <- GRanges(chrom, IRanges(starts, width = opt$tile))
  tile_hits <- findOverlaps(loci, tiles)
  rows_by_tile <- split(queryHits(tile_hits), subjectHits(tile_hits))
  tstart <- start(tiles); tend <- end(tiles)
  chr_tiles <- rbindlist(lapply(names(rows_by_tile), function(tk) {
    ti <- as.integer(tk)
    rows <- rows_by_tile[[tk]]
    if (length(rows) < opt$min_cg) return(NULL)
    st <- region_stats(rows)
    if (is.na(st$frac_g1) || is.na(st$frac_g2)) return(NULL)
    if (st$pass_g1 < opt$min_valid_per_group ||
        st$pass_g2 < opt$min_valid_per_group) return(NULL)
    data.table(chr = chrom, start = tstart[ti], end = tend[ti],
               n_cpg = length(rows),
               mean_frac_group1 = st$frac_g1, mean_frac_group2 = st$frac_g2,
               # Tiles are the secondary comparability layer; DSS provides group
               # means only for called DMRs, so tile FC necessarily uses pooled
               # fractions. Sign is meaningful (higher group -> correct sign);
               # magnitude is noisier than the DMR log2FC. Treat tiles as a
               # cross-check, not the primary hyper/hypo call.
               log2FC = safe_log2fc(st$frac_g1, st$frac_g2),
               direction = ifelse(st$frac_g1 >= st$frac_g2, "hyper", "hypo"),
               min_pval = suppressWarnings(min(dml_pval[rows], na.rm = TRUE)),
               n_pass_group1 = st$pass_g1, n_pass_group2 = st$pass_g2)
  }), fill = TRUE)
  if (nrow(chr_tiles) > 0) {
    tile_all[[chrom]] <- chr_tiles
    fwrite(chr_tiles, file.path(scratch, sprintf("tiles_%s.tsv.gz", chrom)),
           sep = "\t")   # checkpoint
  }

  rm(cov_mat, meth_mat, loci, dml_pval, region_stats); invisible(gc())
  log_msg("  %s done", chrom)
}

# ============================================================================
# MERGE ACROSS CHROMOSOMES AND WRITE OUTPUTS
# ============================================================================
write_bed <- function(df, path) {
  if (nrow(df) == 0) { file.create(path); return(invisible()) }
  bed <- data.frame(df$chr, df$start - 1L, df$end)   # BED start is 0-based
  fwrite(bed, path, sep = "\t", col.names = FALSE)
}

# ---- DMRs ----
dmrs_all <- if (length(dmr_candidates_all))
  rbindlist(dmr_candidates_all, fill = TRUE) else data.table()
if (nrow(dmrs_all) > 0) {
  fwrite(dmrs_all, file.path(opt$outdir, "dmrs_all_candidates.tsv"), sep = "\t")
  final <- dmrs_all[passes_filter == TRUE]
  log_msg("DMRs passing coverage + log2FC>=%.1f filter: %d / %d",
          opt$log2fc, nrow(final), nrow(dmrs_all))
  fwrite(final, file.path(opt$outdir, "dmrs_filtered.tsv"), sep = "\t")
  write_bed(final[direction == "hyper"], file.path(opt$outdir, "dmrs_hyper.bed"))
  write_bed(final[direction == "hypo"],  file.path(opt$outdir, "dmrs_hypo.bed"))
  log_msg("Wrote hyper (%d) / hypo (%d) BED files",
          nrow(final[direction == "hyper"]), nrow(final[direction == "hypo"]))
} else {
  log_msg("No candidate DMRs on any chromosome. Writing empty outputs.")
  file.create(file.path(opt$outdir, "dmrs_filtered.tsv"))
  file.create(file.path(opt$outdir, "dmrs_hyper.bed"))
  file.create(file.path(opt$outdir, "dmrs_hypo.bed"))
}

# ---- tiles ----
tiles_all <- if (length(tile_all)) rbindlist(tile_all, fill = TRUE) else data.table()
if (nrow(tiles_all) > 0) {
  tiles_all[, passes_filter := abs(log2FC) >= opt$log2fc & min_pval < opt$p]
  # direction already set at tile creation (from pooled fraction ordering)
  fwrite(tiles_all, file.path(opt$outdir, "tiles_300bp_all.tsv.gz"), sep = "\t")
  fwrite(tiles_all[passes_filter == TRUE],
         file.path(opt$outdir, "tiles_300bp_filtered.tsv"), sep = "\t")
  log_msg("Tiles passing filter: %d / %d",
          sum(tiles_all$passes_filter), nrow(tiles_all))
}

log_msg("Done. Outputs in %s", opt$outdir)
