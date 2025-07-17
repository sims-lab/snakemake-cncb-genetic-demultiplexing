rule report_vireo_all:
    input:
        summary=expand("results/vireo/{id}", id=scrnaseq.index.unique()),
    output:
        "reports/report_vireo_all.html",
    conda:
        "../envs/tidyverse.yml"
    threads: 1
    resources:
        mem="2G",
        runtime="5m",
    script:
        "../../notebooks/report_vireo_all.Rmd"
