#!/usr/bin/env Rscript
# =============================================================================
# 03b_go_enrichment.R
# Run with: 
#       Rscript 03b_go_enrichment.R \
#       --genes results/annotation/dmr_gene_list.tsv \
#       --outdir results/annotation \
#       [--pvalue 0.05] [--qvalue 0.2]
#
# =============================================================================

suppressMessages({
  library(optparse)
  library(data.table)
})

opt_list <- list(
  make_option("--genes", type = "character",
              default = "results/annotation/dmr_gene_list.tsv",
              help = "TSV with an 'entrez' column (from 03_annotate_dmrs.R) [%default]"),
  make_option("--outdir", type = "character", default = "results/annotation"),
  make_option("--pvalue", type = "double", default = 0.05,
              help = "enrichGO pvalueCutoff [%default]"),
  make_option("--qvalue", type = "double", default = 0.2,
              help = "enrichGO qvalueCutoff [%default]"),
  make_option("--min-genes", type = "integer", default = 5L, dest = "min_genes",
              help = "Refuse to run below this many genes (underpowered) [%default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     sprintf(...)))

# ------------------------- dependency preflight ------------------------------
for (p in c("clusterProfiler", "org.Hs.eg.db")) {
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("Package '%s' not found. This script must run in the conda env:\n  conda run -n go Rscript 03b_go_enrichment.R ...", p))
}

# ------------------------- load gene list ------------------------------------
if (!file.exists(opt$genes)) stop(sprintf("Gene list not found: %s", opt$genes))
gl <- fread(opt$genes)
if (!"entrez" %in% names(gl))
  stop(sprintf("Expected an 'entrez' column in %s; found: %s",
               opt$genes, paste(names(gl), collapse = ", ")))
genes <- unique(as.character(gl$entrez))
genes <- genes[!is.na(genes) & genes != ""]
log_msg("Loaded %d unique Entrez IDs", length(genes))

if (length(genes) < opt$min_genes)
  stop(sprintf("Only %d genes (< --min-genes %d). Too few for meaningful GO; not running.",
               length(genes), opt$min_genes))

# ------------------------- run enrichGO per ontology -------------------------
suppressMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

any_terms <- FALSE
for (ont in c("BP", "MF", "CC")) {
  eg <- tryCatch(
    enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
             ont = ont, pvalueCutoff = opt$pvalue, qvalueCutoff = opt$qvalue,
             readable = TRUE),
    error = function(e) { log_msg("  enrichGO %s errored: %s", ont,
                                  conditionMessage(e)); NULL })
  df <- if (!is.null(eg)) as.data.frame(eg) else data.frame()
  if (nrow(df) == 0) {
    log_msg("  GO %s: no enriched terms (p<%.2g, q<%.2g)", ont, opt$pvalue, opt$qvalue)
    next
  }
  any_terms <- TRUE
  dt <- as.data.table(df)
  # 'Count' = number of input genes supporting the term. Add a fraction and a
  # flag for single-gene terms, which are the ones to distrust.
  if ("Count" %in% names(dt)) {
    dt[, single_gene_term := Count <= 1L]
    n_single <- sum(dt$single_gene_term)
    if (n_single > 0)
      log_msg("  GO %s: %d/%d terms rest on a SINGLE gene -- treat as noise",
              ont, n_single, nrow(dt))
  }
  outp <- file.path(opt$outdir, sprintf("go_enrichment_%s.tsv", ont))
  fwrite(dt, outp, sep = "\t")
  log_msg("  GO %s: %d terms at p<%.2g -> %s", ont, nrow(dt), opt$pvalue, outp)
}

if (!any_terms) {
  log_msg("No enriched GO terms in any ontology. This is the expected result")
  log_msg("for a small, scattered gene set and is reportable as such.")
  # write a marker so the null result is explicit on disk, not just absent files
  writeLines(
    sprintf("No GO terms enriched at p<%.3g, q<%.3g for %d genes (underpowered).",
            opt$pvalue, opt$qvalue, length(genes)),
    file.path(opt$outdir, "go_enrichment_NONE.txt"))
}

log_msg("Done.")
