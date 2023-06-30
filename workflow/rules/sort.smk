###############################################################################
# @author:      Carl-Eric Wegner
# @affiliation: Küsel Lab - Aquatic Geomicrobiology
#              Friedrich Schiller University of Jena
#
#              carl-eric.wegner@uni-jena.de
#              https://github.com/wegnerce
#              https://www.exploringmicrobes.science
###############################################################################

rule sort_individual:
    input:
        bam="results/02_MAPPED/{sample}_mapped_{REF_SHORT}_trim_bismark_pe.deduplicated.bam",
    output:
        sort="results/02_MAPPED/{sample}_mapped_{REF_SHORT}_trim_bismark_pe.dedup.sorted.bam",
        index="results/02_MAPPED/{sample}_mapped_{REF_SHORT}_trim_bismark_pe.dedup.sorted.bam.bai",
    conda: "../envs/samtools.yaml"
    shell:
        """
        samtools sort {input.bam} > {output.sort}
        samtools index {output.sort}
        """
