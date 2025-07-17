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
        vcf="results/haplotype_caller/{sample}.vcf.gz",
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


rule mrna_bed:
    input:
        gtf=config["genome"]["gtf"],
    output:
        bed="results/mrna_bed/mRNAs.bed",
    message:
        """--- Filter GTF to exons and UTRs."""
    threads: 1
    log:
        "logs/mrna_bed.log",
    shell:
        """awk '$3 == "exon" || $3 == "UTR" {{print $1 "\\t" $4 "\\t" $5}}' {input.gtf} > {output.bed} 2> {log}"""


# index VCF files
# necessary for filtering VCF files by region later
# -----------------------------------------------------
rule bcftools_index_variants_call:
    input:
        "results/haplotype_caller/{sample}.vcf.gz",
    output:
        "results/haplotype_caller/{sample}.vcf.gz.csi",
    log:
        "logs/bcftools_index_variants_call/{sample}.log",
    params:
        extra="",  # optional parameters for bcftools index
    wrapper:
        "v7.2.0/bio/bcftools/index"


# filter reference VCF file for exonic SNVs
# -----------------------------------------------------
rule filter_coding_snvs:
    input:
        "results/haplotype_caller/{sample}.vcf.gz",
        index="results/haplotype_caller/{sample}.vcf.gz.csi",
        regions="results/mrna_bed/mRNAs.bed",
    output:
        "results/filter_coding_snvs/{sample}.vcf.gz",
    message:
        """--- Filtering common variants."""
    log:
        "logs/filter_coding_snvs/{sample}.log",
    params:
        extra="--include 'QUAL>=30' --exclude-types indels",
    threads: 8,
    resources:
        mem=lookup(within=config, dpath="bcftools_view/mem"),
        runtime=lookup(within=config, dpath="bcftools_view/runtime"),
    wrapper:
        "v7.2.0/bio/bcftools/view"


def get_cellranger_count_fastqs(wildcards):
    # throw error if sample not found
    if wildcards.id not in scrnaseq.index:
        raise ValueError(f"Sample {wildcards.id} not found in the sample sheet.")
    sample_fastqs_info=scrnaseq[scrnaseq['id'] == wildcards.id]
    fastqs=sample_fastqs_info['path'].iloc[0]
    return fastqs


def get_cellranger_count_sample_names(wildcards):
    # throw error if sample not found
    if wildcards.id not in scrnaseq.index:
        raise ValueError(f"Sample {wildcards.id} not found in the sample sheet.")
    sample_fastqs_info=scrnaseq[scrnaseq['id'] == wildcards.id]
    sample=sample_fastqs_info['sample'].iloc[0]
    return sample


rule cellranger_count:
    input:
        lambda wildcards: get_cellranger_count_fastqs(wildcards),
    output:
        bam="results/cellranger_count/{id}/outs/possorted_genome_bam.bam",
    params:
        genome=config["genome"]["cellranger"],
        sample=lambda wildcards: get_cellranger_count_sample_names(wildcards),
    message:
        """--- Run cellranger count."""
    threads: lookup(within=config, dpath="cellranger_count/cpus")
    resources:
        mem_gb=lookup(within=config, dpath="cellranger_count/mem_gb"),
        mem=str(lookup(within=config, dpath="cellranger_count/mem_gb")) + "G",
        runtime=lookup(within=config, dpath="cellranger_count/runtime"),
        tmpdir="results/cellranger_count/.tmp",
    shell:
        "cd $TMPDIR && "
        " cellranger count --id={wildcards.id} "
        "   --transcriptome={params.genome} "
        "   --fastqs={input} "
        "   --sample={params.sample} "
        "   --create-bam=true "
        "   --localcores={threads} "
        "   --localmem={resources.mem_gb} && "
        " cd - && "
        " mv $TMPDIR/{wildcards.id} results/cellranger_count/ "


# index VCF files
# necessary for filtering VCF files by region later
# -----------------------------------------------------
rule bcftools_index_filtered_variants:
    input:
        "results/filter_coding_snvs/{sample}.vcf.gz",
    output:
        "results/filter_coding_snvs/{sample}.vcf.gz.csi",
    log:
        "logs/bcftools_index_filtered_variants/{sample}.log",
    params:
        extra="",  # optional parameters for bcftools index
    wrapper:
        "v7.2.0/bio/bcftools/index"


rule bcftools_merge_filtered_variants:
    input:
        calls=expand("results/filter_coding_snvs/{sample}.vcf.gz", sample=dnaseq.index.unique()),
        idx=expand("results/filter_coding_snvs/{sample}.vcf.gz.csi", sample=dnaseq.index.unique()),
    output:
        "results/bcftools_merge_filtered_variants/all.vcf.gz",
    log:
        "logs/bcftools_merge_filtered_variants.log",
    params:
        uncompressed_bcf=False,
        extra="",  # optional parameters for bcftools concat (except -o)
    wrapper:
        "v7.2.0/bio/bcftools/merge"


rule cellsnp_lite:
    input:
        bam="results/cellranger_count/{id}/outs/possorted_genome_bam.bam",
        barcode="results/cellranger_count/{id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
        vcf="results/bcftools_merge_filtered_variants/all.vcf.gz",
    output:
        dir=directory("results/cellsnp_lite/{id}"),
    params:
        max_depth=lookup(within=config, dpath="cellsnp_lite/max_depth")
    log:
        "logs/cellsnp_lite/{id}.log",
    conda:
        "../envs/cellsnp_lite.yml",
    threads: lookup(within=config, dpath="cellsnp_lite/cpus")
    resources:
        mem=lookup(within=config, dpath="cellsnp_lite/mem"),
        runtime=lookup(within=config, dpath="cellsnp_lite/runtime"),
    shell:
        "cellsnp-lite "
        " -s {input.bam} "
        " -b {input.barcode} "
        " -O {output.dir} "
        " -R {input.vcf} "
        " --maxDEPTH {params.max_depth}"
        " -p {threads} "
        " --minMAF 0.1 "
        " --minCOUNT 20 "
        " --gzip "
        " 2> {log} "


# TODO: add the outputs of the rule below to the Snakefile target
rule vireo:
    input:
        cellsnp="results/cellsnp_lite/{id}",
        vcf="results/bcftools_merge_filtered_variants/all.vcf.gz",
    output:
        directory("results/vireo/{id}"),
    log:
        "results/vireo/{id}.log",
    conda:
        "../envs/vireo.yml"
    threads: lookup(within=config, dpath="vireo/cpus")
    resources:
        mem=lookup(within=config, dpath="vireo/mem"),
        runtime=lookup(within=config, dpath="vireo/runtime"),
    shell:
        "vireo "
        " -c {input.cellsnp}"
        " -d {input.vcf}"
        " -o {output}"
        " 2> {log} "


# rule popscle_dsc:
#     input:
#         bam="results/cellranger_count/{id}/outs/possorted_genome_bam.bam",
#         vcf="results/bcftools_merge_filtered_variants/all.vcf.gz",
#     output:
#         pileup="results/popscle_dsc/{id}.pileup",
#     conda:
#         "../envs/popscle.yml"
#     message:
#         """--- Running popscle dsc-pileup."""
#     log:
#         "logs/popscle_dsc/{id}.log",
#     threads: 1
#     resources:
#         mem=lookup(within=config, dpath="popscle_dsc/mem"),
#         runtime=lookup(within=config, dpath="popscle_dsc/runtime"),
#     shell:
#         "popscle dsc-pileup"
#         " --sam {input.bam}"
#         " --vcf {input.vcf}"
#         " --out {output.pileup} > {log} 2>&1 &&"
#         " touch {output.pileup}"


# https://snakemake-wrappers.readthedocs.io/en/v7.2.0/wrappers/bio/multiqc.html
rule multiqc:
    input:
        expand("results/samtools_sort/{sample}.flagstat", sample=dnaseq.index.unique()),
        expand("results/samtools_sort/{sample}.idxstats", sample=dnaseq.index.unique()),
        expand("results/mark_duplicates/{sample}.metrics.txt", sample=dnaseq.index.unique()),
        expand("results/filter_coding_snvs/{sample}.vcf.gz", sample=dnaseq.index.unique()),
    output:
        "results/multiqc/qc/multiqc.html",
        directory("results/multiqc/qc_data/multiqc_data"),
    params:
        extra="--verbose",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "v7.2.0/bio/multiqc"
