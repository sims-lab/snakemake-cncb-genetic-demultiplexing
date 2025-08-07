rule cellranger_web_summary:
    input:
        "results/cellranger_count/{id}",
    output:
        "reports/cellranger_web_summary/{id}.html",
    shell:
        "cp {input}/outs/web_summary.html {output}"

rule report_vireo_all:
    input:
        vireo=expand("results/vireo/{id}", id=scrnaseq.index.unique()),
        cellranger=expand("results/cellranger_count/{id}", id=scrnaseq.index.unique()),
    output:
        "reports/report_vireo_all.html",
    conda:
        "../envs/report_vireo_all.yml"
    threads: 1
    resources:
        mem="4G",
        runtime="5m",
    script:
        "../../notebooks/report_vireo_all.Rmd"

rule report_dropletqc_all:
    input:
        dropletqc=expand("results/dropletqc_run/{id}.csv", id=scrnaseq.index.unique()),
        vireo=expand("results/vireo/{id}", id=scrnaseq.index.unique()),
        cellranger=expand("results/cellranger_stats/{id}", id=scrnaseq.index.unique()),
        samples=config["samplesheet"],
    output:
        "reports/report_dropletqc_all.html",
    conda:
        "../envs/reports.yml"
    threads: 1
    resources:
        mem="4G",
        runtime="10m",
    script:
        "../../notebooks/report_dropletqc_all.Rmd"

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
        "../envs/report_tidyverse.yml"
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
