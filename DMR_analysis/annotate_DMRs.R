#!/usr/bin/env Rscript
suppressMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(annotatr)
})

opt_list <- list(
  make_option("--dmrs", type = "character",
              default = "results/dmr/dmrs_filtered.tsv"),
  make_option("--outdir", type = "character", default = "results/annotation"),
  make_option("--genome", type = "character", default = "hg38"),
  make_option("--flag-regions", type = "character", default = NULL,
              dest = "flag_regions",
              help = "Optional TSV (chr,start,end) of DMRs to mark as low-confidence"),
  make_option("--go", action = "store_true", default = FALSE,
              help = "Run GO enrichment (needs clusterProfiler + org.Hs.eg.db)"),
  make_option("--go-pvalue", type = "double", default = 0.05, dest = "go_pvalue")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))

# ------------------------- dependency preflight ------------------------------
# annotatr's genic annotations for hg38 pull from these; build_annotations()
# errors cryptically if they're absent, so check up front with a clear message.
needed <- c("TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss) > 0)
  stop(sprintf(paste0("Missing Bioconductor package(s): %s\n",
                      "Install with:\n  BiocManager::install(c(%s))"),
               paste(miss, collapse = ", "),
               paste(sprintf('\"%s\"', miss), collapse = ", ")))

# ------------------------- load DMRs -> GRanges ------------------------------
dmrs <- fread(opt$dmrs)
if (nrow(dmrs) == 0) stop("No DMRs in input.")
log_msg("Loaded %d DMRs", nrow(dmrs))

# mark low-confidence DMRs if a flag file is given (matched on chr:start:end)
dmrs[, flag := "ok"]
if (!is.null(opt$flag_regions) && file.exists(opt$flag_regions)) {
  fl <- fread(opt$flag_regions, header = TRUE)
  key_dmr <- paste(dmrs$chr, dmrs$start, dmrs$end)
  key_fl  <- paste(fl[[1]], fl[[2]], fl[[3]])
  dmrs[key_dmr %in% key_fl, flag := "low_confidence"]
  log_msg("Flagged %d DMRs as low_confidence", sum(dmrs$flag == "low_confidence"))
}

# annotatr wants a GRanges with an 'id' and any data columns preserved as mcols.
regions <- GRanges(
  seqnames = dmrs$chr,
  ranges   = IRanges(start = dmrs$start, end = dmrs$end),
  DMR_id   = sprintf("DMR_%03d", seq_len(nrow(dmrs))),
  direction  = dmrs$direction,
  diff_methy = dmrs$`diff.Methy`,
  log2FC     = dmrs$log2FC,
  flag       = dmrs$flag
)

# ------------------------- build annotations ---------------------------------
# CpG-context + genic features + enhancers for hg38, matching the study's use of
# annotatr's standard annotation set.
annot_codes <- c(
  paste0(opt$genome, "_cpgs"),          # islands, shores, shelves, inter (shortcut)
  paste0(opt$genome, "_basicgenes"),    # 1-5kb, promoters, 5UTR, exons, introns, 3UTR
  paste0(opt$genome, "_genes_intergenic"),
  paste0(opt$genome, "_enhancers_fantom")
)
log_msg("Building annotations: %s", paste(annot_codes, collapse = ", "))
annotations <- build_annotations(genome = opt$genome, annotations = annot_codes)

# ------------------------- annotate ------------------------------------------
log_msg("Annotating regions")
ann <- annotate_regions(regions = regions, annotations = annotations,
                        ignore.strand = TRUE, quiet = TRUE)
ann_dt <- as.data.table(as.data.frame(ann))

# annotate_regions flattens annotation mcols under 'annot.*'
# Key columns: annot.type (e.g. hg38_cpg_islands, hg38_genes_promoters),
# annot.symbol, annot.gene_id (Entrez).
setnames(ann_dt,
         old = grep("^annot\\.", names(ann_dt), value = TRUE),
         new = sub("^annot\\.", "annot_", grep("^annot\\.", names(ann_dt), value = TRUE)))

fwrite(ann_dt, file.path(opt$outdir, "dmr_annotations_full.tsv"), sep = "\t")
log_msg("Wrote full annotation table (%d DMR-annotation pairs)", nrow(ann_dt))

# ------------------------- per-DMR summary -----------------------------------
# Collapse to one row per DMR: which feature types it hits, and which gene(s).
summ <- ann_dt[, .(
  feature_types = paste(sort(unique(annot_type)), collapse = ";"),
  genes = paste(sort(unique(na.omit(annot_symbol))), collapse = ";"),
  entrez = paste(sort(unique(na.omit(as.character(annot_gene_id)))), collapse = ";")
), by = .(DMR_id, seqnames, start, end, direction, flag)]
fwrite(summ, file.path(opt$outdir, "dmr_annotation_summary.tsv"), sep = "\t")
log_msg("Wrote per-DMR summary")

# ------------------------- feature-type counts (hyper vs hypo) ---------------
type_counts <- ann_dt[, .(n_DMRs = uniqueN(DMR_id)),
                      by = .(annot_type, direction)]
type_counts <- dcast(type_counts, annot_type ~ direction,
                     value.var = "n_DMRs", fill = 0)
fwrite(type_counts, file.path(opt$outdir, "annotation_type_counts.tsv"), sep = "\t")
log_msg("Wrote feature-type counts")

# ------------------------- gene list for GO ----------------------------------
genes <- unique(na.omit(as.character(ann_dt$annot_gene_id)))
gene_tbl <- unique(ann_dt[!is.na(annot_gene_id),
                          .(entrez = as.character(annot_gene_id),
                            symbol = annot_symbol)])
fwrite(gene_tbl, file.path(opt$outdir, "dmr_gene_list.tsv"), sep = "\t")
log_msg("Genes hit by DMRs: %d unique Entrez IDs", length(genes))

# ------------------------- GO enrichment (optional) --------------------------
# Off by default because with only ~27 DMRs the gene list is small and GO power
# is limited; enrichment here is exploratory, not confirmatory. Enable with --go.
if (opt$go) {
  ok <- requireNamespace("clusterProfiler", quietly = TRUE) &&
        requireNamespace("org.Hs.eg.db", quietly = TRUE)
  if (!ok) {
    log_msg("clusterProfiler / org.Hs.eg.db not installed; skipping GO.")
  } else if (length(genes) < 5) {
    log_msg("Only %d genes -- too few for meaningful GO enrichment; skipping.",
            length(genes))
  } else {
    log_msg("Running GO enrichment (clusterProfiler) on %d genes", length(genes))
    for (ont in c("BP", "MF", "CC")) {
      eg <- clusterProfiler::enrichGO(
        gene = genes, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID", ont = ont,
        pvalueCutoff = opt$go_pvalue, qvalueCutoff = 0.2,
        readable = TRUE)
      if (!is.null(eg) && nrow(as.data.frame(eg)) > 0) {
        fwrite(as.data.table(as.data.frame(eg)),
               file.path(opt$outdir, sprintf("go_enrichment_%s.tsv", ont)),
               sep = "\t")
        log_msg("  GO %s: %d terms at p<%.2g", ont, nrow(as.data.frame(eg)),
                opt$go_pvalue)
      } else {
        log_msg("  GO %s: no enriched terms", ont)
      }
    }
  }
}

log_msg("Done. Outputs in %s", opt$outdir)
