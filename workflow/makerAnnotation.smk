##############################################################
#  workflow: makerAnnotation.smk
#  author: Guo-Song Jia
#  last edited: 2021.12.08
#  description: Snakemake workflow for in FYRAflow-annotation submodule. 
#               Do de novo genome annotation using maker to genome assemblies. 
##############################################################

import csv
import json
import os

configfile: "config/annotation_config.yaml"

species = config['SPECIES']
working_dir = config["WORKPATH"]
assembly_path = config['AssemblyPATH']

# Generate sample list file
with open(str(config["SAMPLES"]),'r') as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = [row[0] for row in reader]

# Generate annotation preparation files.
protein_evidence_files = os.getcwd() + "/data/uniprot_sprot_fungi.fasta"
teprotein = os.getcwd() + "/data/te_proteins.fasta"
if str(species) == "S_pombe":
    rmlib = os.getcwd() + "/data/S_pombe.telib.fa"
    est_evidence = os.getcwd() + "/data/S_pombe.EST.fa"
    snaphmm = os.getcwd() + "/data/S_pombe.hmm"
    augustus_species = "schizosaccharomyces_pombe"
# elif str(species) == "S_osmophilus":

rule all:
    input:
        expand(working_dir + "/Maker/{sample}/{sample}.fasta",sample=sampleList),
        expand(working_dir + "/Maker/{sample}/maker_exe.ctl",sample=sampleList),
        expand(working_dir + "/Maker/{sample}/maker_bopts.ctl",sample=sampleList),
        expand(working_dir + "/Maker/{sample}/maker_evm.ctl",sample=sampleList),
        expand(working_dir + "/Maker/{sample}/maker_opts.ctl",sample=sampleList),
        expand(working_dir + "/Maker/{sample}/{sample}.maker.gff3",sample=sampleList),
        expand(working_dir + "/Maker/{sample}/{sample}.maker.nuclear.gff3",sample=sampleList),
        expand(working_dir + "/Maker/{sample}/{sample}.all.gff3",sample=sampleList),
        expand(working_dir + "/Maker/{sample}/{sample}.protein_evidence.gff3",sample=sampleList),
        expand(working_dir + "/Maker/{sample}/{sample}.tRNA.gff3",sample=sampleList)


# Copy the genome assembly fasta file to the working directory and unify assembly file name.
rule copy_assembly:
    output: 
        genome_assembly = working_dir + "/Maker/{sample}/{sample}.fasta"
    run: 
        shell("cp {assembly_path}/{wildcards.sample}*.fasta {output.genome_assembly}")

# Generate common maker config files
rule generateMakerCommonConfig:
    input:
        genome_assembly = working_dir + "/Maker/{sample}/{sample}.fasta"
    output:
        maker_exe_file = working_dir + "/Maker/{sample}/maker_exe.ctl",
        maker_bopts_file = working_dir + "/Maker/{sample}/maker_bopts.ctl",
        maker_evm_file = working_dir + "/Maker/{sample}/maker_evm.ctl"
    shell:
        """
        cd {working_dir}/Maker/{wildcards.sample}/
        maker -CTL
        rm {working_dir}/Maker/{wildcards.sample}/maker_opts.ctl
        """
# Generate strain specific maker config files
rule generateMakerStrainSpecificConfig:
    input:
        maker_exe_file = working_dir + "/Maker/{sample}/maker_exe.ctl",
        maker_bopts_file = working_dir + "/Maker/{sample}/maker_bopts.ctl",
        maker_evm_file = working_dir + "/Maker/{sample}/maker_evm.ctl",
        genome_assembly = working_dir + "/Maker/{sample}/{sample}.fasta"
    output:
        maker_opts_file = working_dir + "/Maker/{sample}/maker_opts.ctl"
    run:
        with open(output.maker_opts_file,'w') as optsFile:
            optsFile.write("genome="+str(working_dir)+"/Maker/"+str(wildcards.sample)+"/"+str(wildcards.sample)+".fasta\n")
            optsFile.write("organism_type=eukaryotic\n")
            optsFile.write("est="+str(est_evidence)+"\n")
            optsFile.write("protein="+str(protein_evidence_files)+"\n")
            optsFile.write("model_org=\n")
            optsFile.write("rmlib="+str(rmlib)+"\n")
            optsFile.write("repeat_protein="+str(teprotein)+"\n")
            optsFile.write("snaphmm="+str(snaphmm)+"\n")
            optsFile.write("augustus_species="+str(augustus_species)+"\n")
            optsFile.write("est2genome=0\n")
            optsFile.write("protein2genome=1\n")
            optsFile.write("trna=1\n")
            optsFile.write("correct_est_fusion=1\n")
            optsFile.write("keep_preds=1\n")

# Run maker
rule runMaker:
    input:
        genome_assembly = working_dir + "/Maker/{sample}/{sample}.fasta",
        maker_exe_file = working_dir + "/Maker/{sample}/maker_exe.ctl",
        maker_bopts_file = working_dir + "/Maker/{sample}/maker_bopts.ctl",
        maker_evm_file = working_dir + "/Maker/{sample}/maker_evm.ctl",
        maker_opts_file = working_dir + "/Maker/{sample}/maker_opts.ctl"
    output:
        maker_gff_output = working_dir + "/Maker/{sample}/{sample}.maker.gff3",
        all_gff_output = working_dir + "/Maker/{sample}/{sample}.all.gff3",
        nuclear_maker_gff_output = working_dir + "/Maker/{sample}/{sample}.maker.nuclear.gff3",
        protein_evidence_output = working_dir + "/Maker/{sample}/{sample}.protein_evidence.gff3",
        trna_gff_output = working_dir + "/Maker/{sample}/{sample}.tRNA.gff3"
    threads:
        36
    shell:
        """
        cd {working_dir}/Maker/{wildcards.sample}/
        mpiexec -n 36 maker -base {wildcards.sample} -fix_nucleotides
        gff3_merge -g -n -d {wildcards.sample}.maker.output/{wildcards.sample}_master_datastore_index.log -o {output.maker_gff_output}
        gff3_merge -n -d {wildcards.sample}.maker.output/{wildcards.sample}_master_datastore_index.log -o {output.all_gff_output}
        cat {output.maker_gff_output} | egrep -v "=trnascan" > {output.nuclear_maker_gff_output}
        cat {output.all_gff_output} |awk '{{if ($2=="protein2genome") {{print $0}} }}' > {output.protein_evidence_output}
        cat {output.maker_gff_output} |egrep "=trnascan" > {output.trna_gff_output}
        """