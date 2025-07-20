rule cellranger_web_summary:
    input:
        "results/cellranger_count/{id}/outs/web_summary.html",
    output:
        "reports/cellranger_web_summary/{id}.html",
    shell:
        "cp {input} {output}"

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

rule report_droplet_filtering_metrics:
    input:
        vireo="results/vireo/{id}",
        dropletqc="results/dropletqc_run/{id}.csv",
        cellranger="results/cellranger_stats/{id}.tsv",
    output:
        "reports/report_droplet_filtering_metrics/{id}.html",
    container:
        "docker://continuumio/miniconda3:24.9.2-0"
    conda:
        "../envs/dropletutils.yml"
    threads: 1
    resources:
        mem="8G",
        runtime="10m",
    script:
        "../../notebooks/report_droplet_filtering_metrics.Rmd"

rule report_index:
    input:
        expand("reports/cellranger_web_summary/{id}.html", id=scrnaseq.index.unique()),
        "reports/report_vireo_all.html",
        expand("reports/report_droplet_filtering_metrics/{id}.html", id=scrnaseq.index.unique()),
    output:
        "reports/index.html",
    container:
        "docker://continuumio/miniconda3:24.9.2-0"
    conda:
        "../envs/dropletutils.yml"
    threads: 1
    resources:
        mem="8G",
        runtime="10m",
    script:
        "../../notebooks/report_index.Rmd"
