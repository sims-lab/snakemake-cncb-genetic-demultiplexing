rule cellranger_filtered_stats:
    input:
        cellranger="results/cellranger_count/{id}/outs/filtered_feature_bc_matrix.h5",
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
