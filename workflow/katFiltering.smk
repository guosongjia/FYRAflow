##############################################################
#  workflow: katFiltering.smk
#  author: Guo-Song Jia
#  last edited: 2021.07.10
#  description: Snakemake workflow for in FYRAflow-variant submodule. 
#               Filter out reads containing low-frequency kmer in FASTP processed fastq files using KAT.
##############################################################

import csv
import json

configfile: "config/variant_config.yaml"

input_path = config["READSPATH"]
working_dir = config["WORKPATH"]
intermediate_dir = config["WORKPATH"] + "/intermediate_" + config["PROJECT"]
PREFIX = "in.jf27"

with open(str(config["SAMPLE"]),'r') as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = sorted([row[0] for row in reader])

rule all:
    input:
        expand(working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist",sample=sampleList),
        expand(working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist.dist_analysis.json",sample=sampleList),
        expand(working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist.png",sample=sampleList),
        expand(working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.filter.kmer-{prefix}",sample=sampleList,prefix=PREFIX),
        expand(working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R1.fq.gz",sample=sampleList),
        expand(working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R2.fq.gz",sample=sampleList),
        expand(working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist_beforeKAT.pdf",sample=sampleList),
        expand(working_dir + "/02_kmerStatFiltering_kat/{sample}_afterFiltering.kat.hist",sample=sampleList),
        expand(working_dir + "/02_kmerStatFiltering_kat/{sample}_afterFiltering.kat.hist.dist_analysis.json",sample=sampleList),
        expand(working_dir + "/02_kmerStatFiltering_kat/{sample}_afterFiltering.kat.hist.png",sample=sampleList),
        expand(working_dir + "/02_kmerStatFiltering_kat/{sample}_afterFiltering.kat.hist_afterKAT.pdf",sample=sampleList)

# Do decompression to fastq files, because "kat filter seq" need decompressed fastq files as input
rule fq_gunzip:
    input: 
        fastp_read1 = working_dir + "/01_QC_fastp/trimmed/{sample}_fastp_1.fq.gz",
        fastp_read2 = working_dir + "/01_QC_fastp/trimmed/{sample}_fastp_2.fq.gz"
    output: 
        fastp_read1_unzip = temp(intermediate_dir + "/{sample}_fastp_1.fq"),
        fastp_read2_unzip = temp(intermediate_dir + "/{sample}_fastp_2.fq")
    run: 
        shell("gunzip -dc {input.fastp_read1} > {output.fastp_read1_unzip}")
        shell("gunzip -dc {input.fastp_read2} > {output.fastp_read2_unzip}")

# Do kat hist statistic, generate kmer analysis results
rule kat_hist:
    input: 
        fastp_read1 = working_dir + "/01_QC_fastp/trimmed/{sample}_fastp_1.fq.gz",
        fastp_read2 = working_dir + "/01_QC_fastp/trimmed/{sample}_fastp_2.fq.gz"
    output: 
        kat_summary = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist",
        kat_json = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist.dist_analysis.json",
        kat_fig = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist.png",
        modified_kmer_fig = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist_beforeKAT.pdf"
    threads: 32
    log:
        working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist.log"
    run:
        shell("kat hist -o {output.kat_summary} -t {threads} {input.fastp_read1} {input.fastp_read2} > {log} 2>&1")
        with open(output.kat_json,'r') as jsonFile:
            jsonParser = json.load(jsonFile)
            xmax = int(jsonParser['global_maxima']['freq']) * 3
            ymax = int(jsonParser['global_maxima']['count']) + 100000
        shell("Rscript scripts/kat_draw_kat_hist.R {output.kat_summary} beforeKAT {xmax} {ymax}")

# Do kat filter kmer, generate binary kmer index file with ".kmer-in.jf27" suffrix.
rule kat_filtering_index:
    input:
        fastp_read1_unzip = intermediate_dir + "/{sample}_fastp_1.fq",
        fastp_read2_unzip = intermediate_dir + "/{sample}_fastp_2.fq",
        kat_json = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist.dist_analysis.json"
    output:
        kmer_index = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.filter.kmer-{prefix}"
    params:
        outdir_index = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.filter.kmer"
    threads: 32
    run:
        with open(input.kat_json,'r') as jsonFile:
            jsonParser = json.load(jsonFile)
            lcfactor = int(int(jsonParser['global_maxima']['freq']) * 0.2)
        shell("kat filter kmer -o {params.outdir_index} -t {threads} --low_count {lcfactor} {input.fastp_read1_unzip} {input.fastp_read2_unzip}")
        
# Do kat filter seq, generate filtered fastq files based on binary kmer index 
rule kat_filtering_seq:
    input:
        fastp_read1_unzip = intermediate_dir + "/{sample}_fastp_1.fq",
        fastp_read2_unzip = intermediate_dir + "/{sample}_fastp_2.fq",
        kmer_index = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.filter.kmer-in.jf27"
    output:
        filtered_fastp_read1 = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R1.fq",
        filtered_fastp_read2 = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R2.fq"
    params:
        kmer_seq = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}"
    threads: 32
    log:
        working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.filtering_seq.log"
    run: 
        shell("kat filter seq -o {params.kmer_seq} -T 1 -t {threads} --seq {input.fastp_read1_unzip} --seq2 {input.fastp_read2_unzip} {input.kmer_index} > {log} 2>&1")

# Compress the kmer filtered fastq files
rule fq_gzip:
    input:
        filtered_fastp_read1 = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R1.fq",
        filtered_fastp_read2 = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R2.fq"
    output:
        gziped_filtered_fastp_read1 = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R1.fq.gz",
        gziped_filtered_fastp_read2 = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R2.fq.gz"
    run:
        shell("gzip {input.filtered_fastp_read1}")
        shell("gzip {input.filtered_fastp_read2}")

# Do kat hist statistic after kat filtering, generate kmer analysis results
rule kat_hist_again:
    input: 
        gziped_filtered_fastp_read1 = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R1.fq.gz",
        gziped_filtered_fastp_read2 = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R2.fq.gz"
    output:
        kat_summary_after = working_dir + "/02_kmerStatFiltering_kat/{sample}_afterFiltering.kat.hist",
        kat_json_after = working_dir + "/02_kmerStatFiltering_kat/{sample}_afterFiltering.kat.hist.dist_analysis.json",
        kat_fig_after = working_dir + "/02_kmerStatFiltering_kat/{sample}_afterFiltering.kat.hist.png",
        modified_kmer_fig_after = working_dir + "/02_kmerStatFiltering_kat/{sample}_afterFiltering.kat.hist_afterKAT.pdf"
    threads: 32
    run:
        shell("kat hist -o {output.kat_summary_after} -t {threads} {input.gziped_filtered_fastp_read1} {input.gziped_filtered_fastp_read2}")
        with open(output.kat_json_after,'r') as jsonFile:
            jsonParser = json.load(jsonFile)
            xmax = int(jsonParser['global_maxima']['freq']) * 3
            ymax = int(jsonParser['global_maxima']['count']) + 100000
        shell("Rscript scripts/kat_draw_kat_hist.R {output.kat_summary_after} afterKAT {xmax} {ymax}")