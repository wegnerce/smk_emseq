###############################################################################
# @author:      Carl-Eric Wegner
# @affiliation: Küsel Lab - Aquatic Geomicrobiology
#              Friedrich Schiller University of Jena
#
#              carl-eric.wegner@uni-jena.de
#              https://github.com/wegnerce
#              https://www.exploringmicrobes.science
###############################################################################


rule fastqc_raw:
    # generate QC reports for the raw data
    input:
        read=raw_data_dir + "/{sample}_{pair}.fastq.gz",
    output:
        qual="logs/fastqc/raw/{sample}_{pair}_fastqc.html",
        zip="logs/fastqc/raw/{sample}_{pair}_fastqc.zip",
    resources:
        mem_mb=2000,
    conda:
        "../envs/fastqc.yaml"
    threads: 2
    shell:
        """
        fastqc {input.read} -t {threads} -f fastq --outdir logs/fastqc/raw
        """


rule bbduk_trim:
    # adapter removal and sequence trimming with bbduk
    input:
        read1=raw_data_dir + "/{sample}_" + PAIRS[0] + ".fastq.gz",
        read2=raw_data_dir + "/{sample}_" + PAIRS[1] + ".fastq.gz",
    output:
        read1="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[0] + ".fastq.gz",
        read2="results/01_TRIMMED/{sample}_trimmed_" + PAIRS[1] + ".fastq.gz",
        trim_stats="logs/bbduk/{sample}_stats_QC.txt",
    resources:
        mem_mb=8000,
    conda:
        "../envs/bbmap.yaml"
    params:
        adapter=ADAPTER,
    threads: 4
    shell:
        """
        bbduk.sh -Xmx8g in1={input.read1} in2={input.read2} out1={output.read1} out2={output.read2} stats={output.trim_stats} minlen=75 qtrim=rl trimq=20 ktrim=r k=25 mink=11 ref={params.adapter} hdist=1
        """


rule fastqc_trim:
    # generate QC reports for the trimmed data
    input:
        read="results/01_TRIMMED/{sample}_trimmed_{pair}.fastq.gz",
    output:
        qual="logs/fastqc/trimmed/{sample}_trimmed_{pair}_fastqc.html",
        zip="logs/fastqc/trimmed/{sample}_trimmed_{pair}_fastqc.zip",
    resources:
        mem_mb=2000,
    conda:
        "../envs/fastqc.yaml"
    threads: 2
    shell:
        """
        fastqc {input.read} -t {threads} -f fastq --outdir logs/fastqc/trimmed
        """
