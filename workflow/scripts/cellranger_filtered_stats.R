library(DropletUtils)
library(tidyverse)

cellranger_h5 <- snakemake@input[[1]]

print(snakemake@input)
print(snakemake@output)

message("Loading H5 file ...")
sce <- read10xCounts(cellranger_h5, col.names = TRUE)
message("Done")

message("Computing stats ...")
tmp_barcodes <- colnames(sce)
tmp_sum <- colSums(counts(sce))
tmp_detected <- colSums(counts(sce) > 0)
message("Done")

rm(sce)

message("Assembling tibble ...")
sce_stats <- tibble(
  barcode = tmp_barcodes,
  sum = tmp_sum,
  detected = tmp_detected
)
message("Done")

message("Writing to file")
write_tsv(sce_stats, snakemake@output[[1]])
message("Done")
