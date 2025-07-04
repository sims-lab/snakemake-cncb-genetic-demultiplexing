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
