# BiocManager::install(version = "3.20")

# options(repos = BiocManager::repositories())

# install.packages("devtools")

dir.create(snakemake@output[["path"]])
.libPaths(snakemake@output[["path"]])

devtools::install_github("powellgenomicslab/DropletQC", Ncpus = snakemake@threads)

# library(DropletQC)
# nf1 <- nuclear_fraction_tags(
#     outs = system.file("extdata", "outs", package = "DropletQC"),
#      tiles = 1, cores = 1, verbose = FALSE)

# cat(as.character(packageVersion("DropletQC")), file = file.path(snakemake@output, ""))
