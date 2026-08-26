library(data.table)

chr_list <- commandArgs(trailingOnly=TRUE)[1]
stats_path <- commandArgs(trailingOnly=TRUE)[2]
motif_path <- commandArgs(trailingOnly=TRUE)[3]
out_path <- commandArgs(trailingOnly=TRUE)[4]
sample_name <- commandArgs(trailingOnly=TRUE)[5]
debug <- c()

stats_path <- file.path(stats_path,paste(sample_name,".stats", sep=""))
motif_path <- file.path(motif_path,paste(sample_name,".motif", sep=""))

#### loading chromosome list
chr <- read.table(chr_list, stringsAsFactors = F)

#### loading stat files using data.table for memory efficiency
stats <- fread(stats_path, sep="\t", header=FALSE)
setnames(stats, paste0("V", seq_len(ncol(stats))))
debug <- c(debug, paste("Stats total reads:", nrow(stats)))

#### loading motif files
motif <- fread(motif_path, sep=" ", header=FALSE)
motif[, V3 := toupper(V3)]
setnames(motif, c("V1", "V2", "V3"), c("read_name", "chr_motif", "motif_seq"))
setnames(stats, "V4", "read_name")

# Merge using data.table (much faster and memory efficient)
setkey(stats, read_name)
setkey(motif, read_name)
stats <- stats[motif, nomatch=0]
debug <- c(debug, paste("Reads with both stats and motif:", nrow(stats)))

#### keeping reads mapped on selected chromosomes, with MAPQ > 20, readlength < 700 (no soft-clip restriction)
stats <- stats[V8 > 20 & V1 %in% chr[,1] & V11 < 700 & V11 > 50]
debug <- c(debug, paste("Stats filtered reads:", nrow(stats)))

####obtaining 4bp motif counts
counts <- table(stats$motif_seq)

if (!(file.exists(out_path))){
  dir.create(out_path)
}

# Save as RDS object
saveRDS(counts, file.path(out_path,paste(sample_name,".motif.R", sep="")))

# Save as CSV for easy viewing
motif_df <- data.frame(
  motif = names(counts),
  count = as.numeric(counts),
  frequency = as.numeric(counts) / sum(counts) * 100,
  stringsAsFactors = FALSE
)
motif_df <- motif_df[order(motif_df$count, decreasing = TRUE), ]
write.csv(motif_df, file.path(out_path,paste(sample_name,".motif.csv", sep="")), row.names = FALSE)

# Save debug log
write.table(debug, file.path(out_path,paste(sample_name,".motif.log", sep="")), row.names = F, col.names = F, quote = F)
