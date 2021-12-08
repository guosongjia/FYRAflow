##############################################################
#  workflow: liftoffTransfer.smk
#  author: Guo-Song Jia
#  last edited: 2021.12.05
#  description: Snakemake workflow for in FYRAflow-annotation submodule. 
#               Do annotation transfer from the reference genome to genome assemblies using Liftoff.
#               Reliable annotations (labeled as "valid_ORFs=1") will be retained. 
##############################################################

import csv
import json

configfile: "config/annotation_config.yaml"

species = config['SPECIES']
working_dir = config["WORKPATH"]
assembly_path = config['AssemblyPATH']
assembly_type = config['AssemblyTYPE']

with open(str(config["SAMPLES"]),'r') as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = [row[0] for row in reader]

## Generate reference genome index files
species_genome = "genome/" + species + "/" + species + ".fasta"
species_annotation = "genome/" + species + "/" + species + ".gff3"

rule all:
    input:
        expand(working_dir + "/Liftoff/{sample}.liftoff.gff3",sample=sampleList),
        expand(working_dir + "/Liftoff/{sample}.liftoff.reliable.gff3",sample=sampleList),
        expand(working_dir + "/Liftoff/{sample}.unmapped_features.txt",sample=sampleList)

# Copy the genome assembly fasta file to the working directory and unify assembly file name.
rule copy_assembly:
    output: 
        genome_assembly = temp(working_dir + "/{sample}.fasta")
    run: 
        shell("cp {assembly_path}/{wildcards.sample}*.fasta {output.genome_assembly}")

# Perform liftoff annotation
rule do_liftoff_transfer:
    input: 
        genome_assembly = working_dir + "/{sample}.fasta"
    output: 
        liftoff_annotation = working_dir + "/Liftoff/{sample}.liftoff.gff3",
        liftoff_reliable_annotation = working_dir + "/Liftoff/{sample}.liftoff.reliable.gff3",
        unmapped_file = working_dir + "/Liftoff/{sample}.unmapped_features.txt",
        fai_file = temp(working_dir + "/{sample}.fasta.fai"),
        mmi_file = temp(working_dir + "/{sample}.fasta.mmi")
    threads:
        16
    shell:
        """ 
        liftoff -g {species_annotation} -o {output.liftoff_annotation} -p {threads} -dir /tmp -u {output.unmapped_file} {input.genome_assembly} {species_genome}
        if [[ "{species}" != "S_osmophilus_final" ]];then
        list=$(cat {output.liftoff_annotation}|grep valid_ORFs=0|awk '{{print $9}}'|cut -d ";" -f 1|cut -d "=" -f 2|cut -d ":" -f 2|tr "\\n" "|"|sed 's/.$//')
        cat {output.liftoff_annotation}|grep -E -v ${{list}} > {output.liftoff_reliable_annotation}
        elif [[ "{species}" == "S_osmophilus_final" ]];then
        list=$(cat {output.liftoff_annotation}|grep valid_ORFs=0|awk '{{print $9}}'|cut -d ";" -f 1|cut -d "=" -f 2|tr "\\n" "|"|sed 's/.$//')
        cat {output.liftoff_annotation}|grep -E -v ${{list}} > {output.liftoff_reliable_annotation}
        fi
        """