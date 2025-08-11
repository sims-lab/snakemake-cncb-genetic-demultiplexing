options(repos = c(CRAN = "https://cran.rstudio.com/"))
dir.create(snakemake@output[["path"]])
.libPaths(snakemake@output[["path"]])
install.packages("devtools")
devtools::install_github(
  "powellgenomicslab/DropletQC",
  Ncpus = snakemake@threads
)
