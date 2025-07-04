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
# https://snakemake-wrappers.readthedocs.io/en/v7.2.0/wrappers/bio/bwa/mem.html
rule bwa_mem:
    input:
        reads=multiext(
            "results/concatenate_fastqs/{sample}/read",
            "1.fastq.gz",
            "2.fastq.gz"
        ),
        idx=multiext(config["genome"]["fasta"], ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "results/bwa_mem/{sample}.bam",
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


# https://snakemake-wrappers.readthedocs.io/en/v7.2.0/wrappers/bio/picard/markduplicates.html
rule mark_duplicates:
    input:
        bams="results/bwa_mem/{sample}.bam",
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="results/mark_duplicates/{sample}.bam",
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


# https://snakemake-wrappers.readthedocs.io/en/v7.2.0/wrappers/bio/samtools/sort.html
rule samtools_sort:
    input:
        "results/mark_duplicates/{sample}.bam",
    output:
        "results/samtools_sort/{sample}.bam",
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


# https://snakemake-wrappers.readthedocs.io/en/v7.2.0/wrappers/bio/samtools/index.html
rule samtools_index:
    input:
        "results/samtools_sort/{sample}.bam",
    output:
        "results/samtools_sort/{sample}.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v7.2.0/bio/samtools/index"


# https://snakemake-wrappers.readthedocs.io/en/v7.2.0/wrappers/bio/samtools/flagstat.html
rule samtools_flagstat:
    input:
        "results/samtools_sort/{sample}.bam",
    output:
        "results/samtools_sort/{sample}.flagstat",
    log:
        "logs/samtools_flagstat/{sample}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v7.2.0/bio/samtools/flagstat"


# https://snakemake-wrappers.readthedocs.io/en/v7.2.0/wrappers/bio/samtools/idxstats.html
rule samtools_idxstats:
    input:
        bam="results/samtools_sort/{sample}.bam",
        idx="results/samtools_sort/{sample}.bam.bai",
    output:
        "results/samtools_sort/{sample}.idxstats",
    log:
        "logs/samtools_idxstats/{sample}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v7.2.0/bio/samtools/idxstats"


# https://snakemake-wrappers.readthedocs.io/en/v7.2.0/wrappers/bio/gatk/haplotypecaller.html
rule haplotype_caller:
    input:
        # single or list of bam files
        bam="results/samtools_sort/{sample}.bam",
        ref=config["genome"]["fasta"],
    output:
        vcf="calls/{sample}.vcf",
    log:
        "logs/haplotype_caller/{sample}.log",
    params:
        extra="",  # optional
        java_opts="",  # optional
    threads: 8
    resources:
        mem=lookup(within=config, dpath="haplotype_caller/mem"),
        runtime=lookup(within=config, dpath="haplotype_caller/runtime"),
    wrapper:
        "v7.2.0/bio/gatk/haplotypecaller"


# https://snakemake-wrappers.readthedocs.io/en/v7.2.0/wrappers/bio/multiqc.html
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
