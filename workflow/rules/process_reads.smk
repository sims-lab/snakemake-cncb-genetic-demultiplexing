# ----------------------------------------------------- #
# EXAMPLE WORKFLOW                                      #
# ----------------------------------------------------- #


def get_sample_fastqs(wildcards):
    sample_fastqs_info=dnaseq[dnaseq['sample'] == wildcards.sample]
    fq1=sample_fastqs_info['fq1'].tolist()
    fq2=sample_fastqs_info['fq2'].tolist()
    return {
        'fq1': fq1,
        'fq2': fq2,
    }


# concatenate fastq files for each sample
# -----------------------------------------------------
rule concatenate_fastqs:
    input:
        unpack(get_sample_fastqs),
    output:
        f1="results/concatenate_fastqs/{sample}/read1.fastq.gz",
        f2="results/concatenate_fastqs/{sample}/read2.fastq.gz",
    log:
        "logs/concatenate_fastqs/{sample}.log",
    threads: 1
    shell:
        "cat {input.fq1} > {output.f1} && "
        "cat {input.fq2} > {output.f2}"


# map reads to a pre-built index using BWA-MEM
# -----------------------------------------------------
rule bwa_mem:
    input:
        reads=multiext(
            "results/concatenate_fastqs/{sample}/read",
            "1.fastq.gz",
            "2.fastq.gz"
        ),
        idx=multiext(config["genome"]["index"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "results/bwa_mem/{sample}.unsorted.bam",
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    threads: 8
    resources:
        mem=lookup(within=config, dpath="bwa_mem/mem"),
        runtime=lookup(within=config, dpath="bwa_mem/runtime"),
    wrapper:
        "v7.2.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        bams="results/bwa_mem/{sample}.unsorted.bam",
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="results/mark_duplicates/{sample}.unsorted.bam",
        metrics="results/mark_duplicates/{sample}.metrics.txt",
    log:
        "logs/mark_duplicates/{sample}.log",
    params:
        extra="--VALIDATION_STRINGENCY SILENT "
              "--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 "
              "--ASSUME_SORT_ORDER queryname "
              "--CLEAR_DT false "
              "--ADD_PG_TAG_TO_READS false",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    threads: 8
    resources:
        mem=lookup(within=config, dpath="mark_duplicates/mem"),
        runtime=lookup(within=config, dpath="mark_duplicates/runtime"),
    wrapper:
        "v7.2.0/bio/picard/markduplicates"

rule samtools_sort:
    input:
        "results/mark_duplicates/{sample}.unsorted.bam",
    output:
        "results/samtools_sort/{sample}.duplicate_marked.bam",
    log:
        "logs/samtools_sort/{sample}.log",
    params:
        extra="-m 2G",
    threads: 8
    resources:
        mem=lookup(within=config, dpath="samtools_sort/mem"),
        runtime=lookup(within=config, dpath="samtools_sort/runtime"),
    wrapper:
        "v7.2.0/bio/samtools/sort"


rule samtools_index:
    input:
        "results/samtools_sort/{sample}.duplicate_marked.bam",
    output:
        "results/samtools_sort/{sample}.duplicate_marked.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v7.2.0/bio/samtools/index"


rule samtools_flagstat:
    input:
        "results/samtools_sort/{sample}.duplicate_marked.bam",
    output:
        "results/samtools_sort/{sample}.flagstat",
    log:
        "logs/samtools_flagstat/{sample}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v7.2.0/bio/samtools/flagstat"


rule samtools_idxstats:
    input:
        bam="results/samtools_sort/{sample}.duplicate_marked.bam",
        idx="results/samtools_sort/{sample}.duplicate_marked.bam.bai",
    output:
        "results/samtools_sort/{sample}.idxstats",
    log:
        "logs/samtools_idxstats/{sample}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v7.2.0/bio/samtools/idxstats"


rule multiqc:
    input:
        expand("results/samtools_sort/{sample}.flagstat", sample=dnaseq.index.unique()),
        expand("results/samtools_sort/{sample}.idxstats", sample=dnaseq.index.unique()),
        expand("results/mark_duplicates/{sample}.metrics.txt", sample=dnaseq.index.unique()),
    output:
        "results/multiqc/qc/multiqc.html",
        directory("results/multiqc/qc_data/multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "v7.2.0/bio/multiqc"
