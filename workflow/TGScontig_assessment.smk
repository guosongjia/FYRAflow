##############################################################
#  workflow: TGScontig_assessment.smk
#  author: Guo-Song Jia
#  last edited: 2021.07.10
#  description: Snakemake workflow for in FYRAflow-assembly submodule.
#               Do Correctness assessment and Continuity assessment to TGS assembly.
##############################################################

import csv
configfile: "config/assembly_config.yaml"

species = config['SPECIES']
assembly_path = config['RAWassemblyPATH']
working_dir = config["WORKPATH"]
PREFIX = "report"

with open(str(config["SAMPLES"]),'r') as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = sorted([row[0] for row in reader])

rule all:
    input: 

# Correctness assessment of TGS assembly using NGS reads
# Generate BWA index, do mapping, sorting and mark-duplicate reads to NGS data
rule bwa_mapping_sorting_indexing_Markduplicates:
    input:
        TGS_contigs = assembly_path + "/{sample}.fasta",
        gziped_filtered_fastp_kat_reads1 = katReads_path + "/{sample}.in.R1.fq.gz",
        gziped_filtered_fastp_kat_reads2 = katReads_path + "/{sample}.in.R2.fq.gz"
    output:
        TGS_contig_amb = temp(assembly_path + "/{sample}.fasta.amb"),
        TGS_contig_ann = temp(assembly_path + "/{sample}.fasta.ann"),
        TGS_contig_bwt = temp(assembly_path + "/{sample}.fasta.bwt"),
        TGS_contig_pac = temp(assembly_path + "/{sample}.fasta.pac"),
        TGS_contig_sa = temp(assembly_path + "/{sample}.fasta.sa"),
        mapped_sam_file = temp(working_dir + "/{sample}.kat_filtered.bwa.sam"),
        mapped_sorted_bam_file = temp(working_dir + "/01C_correctness/{sample}.kat_filtered.bwa.sorted.bam"),
        mapped_sorted_bam_index_file = temp(working_dir + "/01C_correctness/{sample}.kat_filtered.bwa.sorted.bam.bai"),

    threads: 32


# Continuity assessment of TGS assembly 
rule TGS_contigs_quast_assessment:
    input:
        TGS_contigs = assembly_path + "/{sample}.fasta"
    params:
        outdir = working_dir + "/Assessment/02C_continuity_assessment/{sample}_QUAST",
        prefix = PREFIX
    output:
        NGS_raw_contigs_quast_report = working_dir + "/Assessment/02_quast_assessment/{sample}_NGS.raw/{prefix}.tsv"
    threads: 36
    run:
        shell("quast.py -r genome/{species}/{species}.fasta -o {params.outdir} -t {threads} {input.NGS_raw_contigs}")
