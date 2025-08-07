rule cellranger_filtered_stats:
    input:
        cellranger="results/cellranger_count/{id}",
        gtf=config["genome"]["gtf"],
    output:
        "results/cellranger_stats/{id}.tsv",
    conda:
        "../envs/cellranger_filtered_stats.yml"
    threads: 2
    resources:
        mem="8G",
        runtime="10m",
    script:
        "../scripts/cellranger_filtered_stats.R"
