.libPaths(snakemake@input[["rlib"]])

library(DropletQC)

nf1 <- nuclear_fraction_tags(
    outs = dirname(snakemake@input[["bam"]]), # nolint: indentation_linter.
    tiles = snakemake@threads,
    cores = snakemake@threads, verbose = FALSE)

write.csv(nf1, snakemake@output[["csv"]], quote = FALSE)
