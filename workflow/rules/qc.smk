rule cellranger_filtered_stats:
    input:
        "results/cellranger_count/{id}/outs/filtered_feature_bc_matrix.h5",
    output:
        "results/cellranger_stats/{id}.tsv",
    conda:
        "../envs/dropletutils.yml"
    threads: 2
    resources:
        mem="8G",
        runtime="10m",
    script:
        "../scripts/cellranger_filtered_stats.R"
