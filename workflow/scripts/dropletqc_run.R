.libPaths(snakemake@input[["rlib"]])

library(DropletQC)

nf1 <- nuclear_fraction_tags(
  outs = file.path(snakemake@input[["cellranger"]], "outs"),
  tiles = snakemake@threads,
  cores = snakemake@threads, verbose = FALSE
)

write.csv(nf1, snakemake@output[["csv"]], quote = FALSE)
