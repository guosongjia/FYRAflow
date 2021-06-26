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
        expand(assembly_path + "/02_assessment/02.2_quast_assessment/{sample}_NGS.raw/{prefix}.tsv",sample=sampleList,prefix=PREFIX),
        expand(assembly_path + "/02_assessment/02.2_quast_assessment/{sample}_NGS.filtered/{prefix}.tsv",sample=sampleList,prefix=PREFIX),
        expand(assembly_path + "/02_assessment/{sample}_assembly_summary.txt",sample=sampleList)

# Do quast assessment to Raw NGS contigs
rule quast_assessment_raw_contigs:
    input:
        NGS_raw_contigs = assembly_path + "/{sample}_NGS.contigs.fasta"
    params:
        outdir = assembly_path + "/02_assessment/02.2_quast_assessment/{sample}_NGS.raw",
        prefix = PREFIX
    output:
        NGS_raw_contigs_quast_report = assembly_path + "/02_assessment/02.2_quast_assessment/{sample}_NGS.raw/{prefix}.tsv"
    threads: 36
    run:
        shell("quast.py -r genome/{species}/{species}.fasta -o {params.outdir} -t {threads} {input.NGS_raw_contigs}")

# Do quast assessment to KAT-filtered reads filtered NGS contigs
rule quast_assessment_filtered_contigs:
    input:
        NGS_filtered_contigs = assembly_path + "/01_katReads_filtering/filtered_contigs/{sample}_NGS.filtered.contigs.fasta"
    params:
        outdir = assembly_path + "/02_assessment/02.2_quast_assessment/{sample}_NGS.filtered",
        prefix = PREFIX
    output:
        NGS_filtered_contigs_quast_report = assembly_path + "/02_assessment/02.2_quast_assessment/{sample}_NGS.filtered/{prefix}.tsv"
    threads: 36
    run:
        shell("quast.py -r genome/{species}/{species}.fasta -o {params.outdir} -t {threads} {input.NGS_filtered_contigs}")

# raw and filtered NGS contigs QUAST statistics
rule assembly_statistics:
    input: 
        NGS_raw_contigs_quast_report = assembly_path + "/02_assessment/02.2_quast_assessment/{sample}_NGS.raw/report.tsv",
        NGS_filtered_contigs_quast_report = assembly_path + "/02_assessment/02.2_quast_assessment/{sample}_NGS.filtered/report.tsv"
    output:
        assembly_statictics = assembly_path + "/02_assessment/{sample}_assembly_summary.txt"
    run:
        shell("python scripts/assembly_assessment_summary.py quast -r {input.NGS_raw_contigs_quast_report} -f {input.NGS_filtered_contigs_quast_report} -o {output.assembly_statictics}")
