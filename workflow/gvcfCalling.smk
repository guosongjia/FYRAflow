##############################################################
#  workflow: gvcfCalling.smk
#  author: Guo-Song Jia
#  last edited: 2021.07.10
#  description: Snakemake workflow for in FYRAflow-varaint submodule. Generate gVCF files using GATK HaplotypeCaller
##############################################################

import csv
configfile: "config/variant_config.yaml"

input_path = config["READSPATH"]
working_dir = config["WORKPATH"]
species = config["SPECIES"]
## Generate reference genome index files
species_index = "genome/" + species + "/" + species + ".fasta"

with open(str(config["SAMPLE"]),'r') as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = [row[0] for row in reader]

rule all:
    input:
        expand(working_dir + "/04_variantCalling/{sample}.raw.snps.indels.g.vcf",sample=sampleList),
        expand(working_dir + "/04_variantCalling/{sample}.raw.snps.indels.g.vcf.idx",sample=sampleList)

# Do variant calling under HaplotypeCaller "-ERC GVCF" mode to genreate gvcf files per sample. 
rule generate_gvcf:
    input:
        mapped_sorted_rmdup_addGroup_bam_file = working_dir + "/03_readsMapping/{sample}.kat_filtered.bwa.sorted.rmdup.RGadded.bam"
    output:
        gvcf_file = working_dir + "/04_variantCalling/{sample}.raw.snps.indels.g.vcf",
        gvcf_index = working_dir + "/04_variantCalling/{sample}.raw.snps.indels.g.vcf.idx"
    threads: 48
    run:
        shell("gatk HaplotypeCaller --java-options '-Xmx128G' --native-pair-hmm-threads {threads}  -I {input.mapped_sorted_rmdup_addGroup_bam_file} -R {species_index} -ERC GVCF -O {output.gvcf_file}")