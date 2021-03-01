import csv
configfile: "config/main_config.yaml"

input_path = config["READSPATH"]
working_dir = config["WORKPATH"]


with open(str(config["SAMPLE"]),'r') as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = [row[0] for row in reader]

rule all:
    input:
        expand(working_dir + "/01_QC_fastp/trimmed/{sample}_fastp_1.fq.gz",sample=sampleList),
        expand(working_dir + "/01_QC_fastp/trimmed/{sample}_fastp_2.fq.gz",sample=sampleList),
        expand(working_dir + "/01_QC_fastp/report/{sample}.fastp.html",sample=sampleList),
        expand(working_dir + "/01_QC_fastp/report/{sample}.fastp.json",sample=sampleList)

# Do fastp quality control and trimming.
rule qualityControl:
    input: 
        read1 = input_path + "/{sample}_1.fq.gz",
        read2 = input_path + "/{sample}_2.fq.gz"
    output: 
        trimmed_read1 = working_dir + "/01_QC_fastp/trimmed/{sample}_fastp_1.fq.gz", 
        trimmed_read2 = working_dir + "/01_QC_fastp/trimmed/{sample}_fastp_2.fq.gz",
        html = working_dir + "/01_QC_fastp/report/{sample}.fastp.html",
        json = working_dir + "/01_QC_fastp/report/{sample}.fastp.json"
    log:
        working_dir + "/01_QC_fastp/{sample}.fastp_run.log"
    params:
        "--length_required 70"
    threads: 36
    shell:
        "fastp -w {threads} {params} -i {input.read1} -I {input.read2} -o {output.trimmed_read1} -O {output.trimmed_read2} --html {output.html} --json {output.json} > {log} 2>&1" 