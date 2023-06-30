###############################################################################
# @author:      Carl-Eric Wegner
# @affiliation: Küsel Lab - Aquatic Geomicrobiology
#              Friedrich Schiller University of Jena
#
#              carl-eric.wegner@uni-jena.de
#              https://github.com/wegnerce
#              https://www.exploringmicrobes.science
###############################################################################


rule bismark_genome_prep:
    # in silico bisulfite treatment of the ref genome
    input:
        genome_dir="resources/",
    output:
        "resources/genomic_nucleotide_frequencies.txt",
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        bismark_genome_preparation --genomic_composition {input.genome_dir}
        """


rule bismark:
    # based on: https://github.com/seb-mueller/snakemake-bisulfite
    input:
        genome_dir="resources/",
        read1="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[0] + ".fastq.gz",
        read2="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[1] + ".fastq.gz",
    output:
        bam=temp("results/02_MAPPED/{sample}_mapped_{REF_SHORT}_trim_bismark_pe.bam"),
    resources:
        mem_mb=16000,
    log:
        log1="logs/bismark/{sample}_mapped_{REF_SHORT}.log",
        log2="logs/bismark/{sample}_mapped_{REF_SHORT}_trim_bismark_PE_report.txt",
        log3="logs/bismark/{sample}_mapped_{REF_SHORT}_trim_bismark_pe.nucleotide_stats.txt",
    benchmark:
        "logs/bismark/{sample}_mapped_{REF_SHORT}_bismark.benchmark"
    threads: 4
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        bismark -p {threads} --nucleotide_coverage {REF_GENOME_DIR} -1 {input.read1} -2 {input.read2} --basename {wildcards.sample}_mapped_{REF_SHORT}_trim_bismark --output_dir results/02_MAPPED 2> {log.log1}
        mv results/02_MAPPED/{wildcards.sample}_mapped_{REF_SHORT}_trim_bismark_PE_report.txt {log.log2}
        mv results/02_MAPPED/{wildcards.sample}_mapped_{REF_SHORT}_trim_bismark_pe.nucleotide_stats.txt {log.log3}
        """


rule deduplicate:
    # remove alignments with identical mapping positions
    input:
        bam="results/02_MAPPED/{sample}_mapped_{REF_SHORT}_trim_bismark_pe.bam",
    output:
        bam=temp(
            "results/02_MAPPED/{sample}_mapped_{REF_SHORT}_trim_bismark_pe.deduplicated.bam"
        ),
    log:
        "logs/bismark_deduplication/{sample}_mapped_{REF_SHORT}.deduplication.log",
    benchmark:
        "logs/bismark_deduplication/{sample}_mapped_{REF_SHORT}.deduplication.benchmark"
    conda:
        "../envs/bismark.yaml"
    shell:
        """
        deduplicate_bismark --paired --bam {input.bam} --output_dir results/02_MAPPED 2> {log}
        """


rule extract_methylation:
    # generate cytosine methylation calls
    input:
        bam="results/02_MAPPED/{sample}_mapped_{REF_SHORT}_trim_bismark_pe.dedup.sorted.bam",
    output:
        CHH="results/03_METHYL_EX/CHH_context_{sample}_mapped_{REF_SHORT}_trim_bismark_pe.dedup.sorted.txt.gz",
        CHG="results/03_METHYL_EX/CHG_context_{sample}_mapped_{REF_SHORT}_trim_bismark_pe.dedup.sorted.txt.gz",
        CpG="results/03_METHYL_EX/CpG_context_{sample}_mapped_{REF_SHORT}_trim_bismark_pe.dedup.sorted.txt.gz",
    resources:
        mem_mb=16000,
    threads: 4
    conda:
        "../envs/bismark.yaml"
    params:
        ref=REF_GENOME_DIR,
        extra=EXTRA_PARAMS_METH_EX,
    log:
        "logs/bismark_methylation_extraction/{sample}_mapped_{REF_SHORT}.methylation_extractor.log",
    benchmark:
        "logs/bismark_methylation_extraction/{sample}_mapped_{REF_SHORT}.methylation_extractor.benchmark"
    shell:
        """
        bismark_methylation_extractor {params.extra} --gzip --paired-end --multicore {threads} --genome_folder {params.ref} -s {input.bam} --output results/03_METHYL_EX 2> {log}
        """


rule bismark2bedGraph:
    # generate bedGraph plus coverage output
    input:
        methylex="results/03_METHYL_EX/{context}_context_{sample}_mapped_{REF_SHORT}_trim_bismark_pe.dedup.sorted.txt.gz",
    output:
        bedGraph="results/04_COVERAGE/{sample}_mapped_{REF_SHORT}_{context}.bedGraph.gz",
    log:
        "logs/bismark_coverage/{sample}_mapped_{REF_SHORT}_{context}.bismark2bedGraph.log",
    benchmark:
        "logs/bismark_coverage/{sample}_mapped_{REF_SHORT}_{context}.bismark2bedGraph.benchmark"
    resources:
        mem_mb=16000,
    threads: 4
    conda:
        "../envs/bismark.yaml"
    params:
        bedGraph="{sample}_mapped_{REF_SHORT}_{context}.bedGraph",
    shell:
        """
        bismark2bedGraph --CX --dir results/04_COVERAGE -o {params.bedGraph} {input.methylex} 2> {log}
        """
