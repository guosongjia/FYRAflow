import csv
import json

configfile: "config/main_config.yaml"

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
        expand(working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R2.fq.gz",sample=sampleList)

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


rule kat_hist:
    input: 
        fastp_read1 = working_dir + "/01_QC_fastp/trimmed/{sample}_fastp_1.fq.gz",
        fastp_read2 = working_dir + "/01_QC_fastp/trimmed/{sample}_fastp_2.fq.gz"
    output: 
        kat_summary = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist",
        kat_json = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist.dist_analysis.json",
        kat_fig = working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist.png"
    threads: 32
    log:
        working_dir + "/02_kmerStatFiltering_kat/{sample}.kat.hist.log"
    shell:
        "kat hist -o {output.kat_summary} -t {threads} {input.fastp_read1} {input.fastp_read2} > {log} 2>&1"


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


