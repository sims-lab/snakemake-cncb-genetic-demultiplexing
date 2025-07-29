library(DropletUtils)
library(rtracklayer)
library(tidyverse)

cellranger_h5 <- snakemake@input[["cellranger"]]
gtf_file <- snakemake@input[["gtf"]]

message("Loading H5 file ...")
sce <- read10xCounts(cellranger_h5, col.names = TRUE)
message("Done")

message("Loading GTF file ...")
gtf_data <- BiocIO::import(gtf_file)
message("Done")

message("Identifying mitochondrial genes ...")
mt_gene_ids <- unique(mcols(subset(gtf_data, seqnames == "mitochondrion_genome"))[["gene_id"]])
rm(gtf_data)
message(sprintf("Done (%i genes)", length(mt_gene_ids)))

message("Computing stats ...")
tmp_barcodes <- colnames(sce)
tmp_sum <- colSums(counts(sce))
tmp_detected <- colSums(counts(sce) > 0)
tmp_mt_fraction <- colSums(counts(sce)[mt_gene_ids, ]) / colSums(counts(sce))
message("Done")

rm(sce)

message("Assembling tibble ...")
sce_stats <- tibble(
  barcode = tmp_barcodes,
  sum = tmp_sum,
  detected = tmp_detected,
  mt_fraction = tmp_mt_fraction
)
message("Done")

message("Writing to file")
write_tsv(sce_stats, snakemake@output[[1]])
message("Done")
