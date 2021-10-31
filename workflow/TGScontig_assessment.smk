##############################################################
#  workflow: TGScontig_assessment.smk
#  author: Guo-Song Jia
#  last edited: 2021.07.13
#  description: Snakemake workflow for in FYRAflow-assembly submodule.
#               Do Correctness assessment and Continuity assessment to TGS assembly.
##############################################################

import csv
import os
configfile: "config/assembly_config.yaml"

species = config['SPECIES']
assembly_path = config['RAWassemblyPATH']
katReads_path = config['KATREADS']
intermediate_dir = config["WORKPATH"] + "/intermediate_" + config["PROJECT"]
working_dir = config["WORKPATH"]
PREFIX = "report"

with open(str(config["SAMPLES"]),'r') as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = list()
    appendixList = list()
    for row in reader:
        sampleList.append(row[0])
        appendixList.append(row[1])
        sampleList = list(set(sampleList))

os.system("mkdir {0}".format(intermediate_dir))
for strains in sampleList:
    for appendix in appendixList:
        os.system("ln -s {0}/{1}.in.R1.fq.gz {2}/{1}_{3}.in.R1.fq.gz".format(katReads_path,strains,intermediate_dir,appendix))
        os.system("ln -s {0}/{1}.in.R2.fq.gz {2}/{1}_{3}.in.R2.fq.gz".format(katReads_path,strains,intermediate_dir,appendix))

rule all:
    input:
        expand(assembly_path + "/{sample}_{appendix}.fasta.fai",sample=sampleList,appendix=appendixList),
        expand(assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sorted.rmdup.RGadded.bam",sample=sampleList,appendix=appendixList),
        expand(assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sorted.rmdup.RGadded.bam.bai",sample=sampleList,appendix=appendixList),
        expand(working_dir + "/Assessment/01C_correctness/01.1_GATK/{sample}_{appendix}.gatk.vcf",sample=sampleList,appendix=appendixList),
        expand(working_dir + "/Assessment/01C_correctness/01.1_GATK/{sample}_{appendix}.gatk.final.vcf",sample=sampleList,appendix=appendixList),
        expand(working_dir + "/Assessment/01C_correctness/01.3_samtools/{sample}_{appendix}.samtools.vcf",sample=sampleList,appendix=appendixList),
        expand(working_dir + "/Assessment/01C_correctness/01.3_samtools/{sample}_{appendix}.samtools.final.vcf",sample=sampleList,appendix=appendixList)

# Correctness assessment of TGS assembly using NGS reads
## Generate BWA index, do mapping, sorting and mark-duplicate reads to NGS data
rule bwa_mapping_sorting_indexing_Markduplicates:
    input:
        TGS_contigs = assembly_path + "/{sample}_{appendix}.fasta",
        gziped_filtered_fastp_kat_reads1 = intermediate_dir + "/{sample}_{appendix}.in.R1.fq.gz",
        gziped_filtered_fastp_kat_reads2 = intermediate_dir + "/{sample}_{appendix}.in.R2.fq.gz"
    output:
        TGS_contig_fai = assembly_path + "/{sample}_{appendix}.fasta.fai",
        TGS_contig_amb = temp(assembly_path + "/{sample}_{appendix}.fasta.amb"),
        TGS_contig_ann = temp(assembly_path + "/{sample}_{appendix}.fasta.ann"),
        TGS_contig_bwt = temp(assembly_path + "/{sample}_{appendix}.fasta.bwt"),
        TGS_contig_pac = temp(assembly_path + "/{sample}_{appendix}.fasta.pac"),
        TGS_contig_sa = temp(assembly_path + "/{sample}_{appendix}.fasta.sa"),
        mapped_sam_file = temp(assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sam"),
        mapped_sorted_bam_file = temp(assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sorted.bam"),
        mapped_sorted_bam_index_file = temp(assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sorted.bam.bai"),
        mapped_sorted_rmdup_bam_file = temp(assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sorted.rmdup.bam"),
        rmdup_metrics = temp(assembly_path + "/{sample}_{appendix}_kat_filtered.marked_dup_matrics.txt"),
        mapped_sorted_rmdup_addGroup_bam_file = assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sorted.rmdup.RGadded.bam",
        mapped_sorted_rmdup_addGroup_bam_index_file = assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sorted.rmdup.RGadded.bam.bai"
    threads: 32
    params:
        sampleName = "{sample}"
    run:
        shell("samtools faidx {input.TGS_contigs}")
        shell("bwa index {input.TGS_contigs}")
        shell("bwa mem -t {threads} {input.TGS_contigs} {input.gziped_filtered_fastp_kat_reads1} {input.gziped_filtered_fastp_kat_reads2} > {output.mapped_sam_file}")
        shell("samtools sort -o {output.mapped_sorted_bam_file} {output.mapped_sam_file}")
        shell("samtools index {output.mapped_sorted_bam_file}")
        shell("gatk MarkDuplicates --java-options '-Xmx16G' -I {output.mapped_sorted_bam_file} -O {output.mapped_sorted_rmdup_bam_file} -M {output.rmdup_metrics}")
        shell("gatk AddOrReplaceReadGroups -I {output.mapped_sorted_rmdup_bam_file} -O {output.mapped_sorted_rmdup_addGroup_bam_file} -RGID {params.sampleName} -RGPU unit1 -RGSM {params.sampleName} -RGLB {params.sampleName} -RGPL ILLUMINA")
        shell("samtools index {output.mapped_sorted_rmdup_addGroup_bam_file}")
        shell("rm {input.gziped_filtered_fastp_kat_reads1} {input.gziped_filtered_fastp_kat_reads2}")

# ## Do GATK variant calling
rule Do_GATK_calling:
    input:
        TGS_contigs = assembly_path + "/{sample}_{appendix}.fasta",
        mapped_sorted_rmdup_addGroup_bam_file = assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sorted.rmdup.RGadded.bam"
    output:
        TGS_contig_dict = assembly_path + "/{sample}_{appendix}.dict",
        gatk_vcf_file = working_dir + "/Assessment/01C_correctness/01.1_GATK/{sample}_{appendix}.gatk.vcf",
        gatk_reliable_vcf_file = working_dir + "/Assessment/01C_correctness/01.1_GATK/{sample}_{appendix}.gatk.final.vcf"
    params:
        sampleName = "{sample}_{appendix}",
        saving_PATH = working_dir + "/Assessment/01C_correctness/01.1_GATK"
    run:
        shell("gatk CreateSequenceDictionary -R {input.TGS_contigs} -O {output.TGS_contig_dict}")
        shell("gatk HaplotypeCaller -R {input.TGS_contigs} -I {input.mapped_sorted_rmdup_addGroup_bam_file} -O {output.gatk_vcf_file}")
        shell("bash scripts/reliable_variant_select.sh gatk {output.gatk_vcf_file} {params.saving_PATH} {params.sampleName}")

## Do Deepvariant variant calling
# rule Do_deepvariant_calling:
#     input:
#         TGS_contigs = assembly_path + "/{sample}_{appendix}.fasta",
#         mapped_sorted_rmdup_addGroup_bam_file = assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sorted.rmdup.RGadded.bam"
#     output:
#         # deepvariant_vcf_file = working_dir + "/Assessment/01C_correctness/01.2_Deepvariant/{sample}_{appendix}.deepvariant.vcf",
#         # deepvariant_report_file = working_dir + "/Assessment/01C_correctness/01.2_Deepvariant/{sample}_{appendix}.deepvariant.visual_report.html"
#         # deepvariant_reliable_vcf_file = working_dir + "/Assessment/01C_correctness/01.2_Deepvariant/{sample}_{appendix}.deepvariant.final.vcf"
#     params:
#         sampleName = "{sample}_{appendix}",
#         saving_PATH = working_dir + "/Assessment/01C_correctness/01.2_Deepvariant"
#     threads: 32
#     run:
#         # shell('docker run -v "{assembly_path}/":"/input" -v "{params.saving_PATH}/":"/output" \
#         #             google/deepvariant:"1.1.0" /opt/deepvariant/bin/run_deepvariant --model_type=WGS \
#         #             --ref=/input/{params.sampleName}.fasta --reads=/input/{params.sampleName}.kat_filtered.bwa.sorted.rmdup.RGadded.bam \
#         #             --output_vcf=/output/{params.sampleName}.deepvariant.vcf ')
#         # shell("bash scripts/reliable_variant_select.sh deepvariant {output.gatk_vcf_file} {params.saving_PATH} {params.sampleName}")

# ## Do samtools variant calling
rule Do_samtools_calling:
    input:
        TGS_contigs = assembly_path + "/{sample}_{appendix}.fasta",
        mapped_sorted_rmdup_addGroup_bam_file = assembly_path + "/{sample}_{appendix}.kat_filtered.bwa.sorted.rmdup.RGadded.bam"
    output:
        samtools_vcf_file = working_dir + "/Assessment/01C_correctness/01.3_samtools/{sample}_{appendix}.samtools.vcf",
        samtools_reliable_vcf_file = working_dir + "/Assessment/01C_correctness/01.3_samtools/{sample}_{appendix}.samtools.final.vcf"
    params:
        sampleName = "{sample}_{appendix}",
        saving_PATH = working_dir + "/Assessment/01C_correctness/01.3_samtools"
    run:
        shell('bcftools mpileup -Ov  -f {input.TGS_contigs} {input.mapped_sorted_rmdup_addGroup_bam_file} | bcftools call -mv -Ov -o {output.samtools_vcf_file}')
        shell("bash scripts/reliable_variant_select.sh samtools {output.samtools_vcf_file} {params.saving_PATH} {params.sampleName}")

# # Continuity assessment of TGS assembly 
rule TGS_contigs_quast_assessment:
    input:
        TGS_contigs = assembly_path + "/{sample}_{appendix}.fasta"
    params:
        outdir = working_dir + "/Assessment/02C_continuity_assessment/{sample}_{appendix}_QUAST",
        prefix = PREFIX
    output:
        NGS_raw_contigs_quast_report = working_dir + "/Assessment/02_quast_assessment/{sample}_{appendix}_QUAST/{prefix}.tsv"
    threads: 36
    run:
        shell("quast.py -r genome/{species}/{species}.fasta -o {params.outdir} -t {threads} {input.NGS_raw_contigs}")
