import csv
import os
configfile: "config/variant_config.yaml"

project = config["PROJECT"]
working_dir = config["WORKPATH"]
species = config["SPECIES"]
## Generate reference genome index files
species_index = "genome/" + species + "/"+ species +".fasta"

gvcfPath = working_dir + "/04_variantCalling"
gvcfList = [files for files in os.listdir(gvcfPath) if os.path.splitext(files)[1]==".vcf"]

rule all:
    input:
        working_dir + "/05_jointCalling/gvcf.list",
        working_dir + "/05_jointCalling/" + project + ".raw.snps.indels.vcf",
        working_dir + "/05_jointCalling/SNP_analysis/" + project + ".filtered.snps.vcf",
        working_dir + "/05_jointCalling/INDEL_analysis/" + project + ".filtered.indels.vcf",
        working_dir + "/05_jointCalling/SNP_analysis/" + project + ".filtered.PASS.snps.vcf",
        working_dir + "/05_jointCalling/INDEL_analysis/" + project + ".filtered.PASS.indels.vcf",
        working_dir + "/05_jointCalling/SNP_analysis/" + project + ".filtered.PASS.unphased.snps.vcf"

# Obtain gvcf file list for gvcf combine using custom python scripts
rule get_gvcf_list:
    output:
        gvcf_list_file = working_dir + "/05_jointCalling/gvcf.list"
    run:
        with open(str(output.gvcf_list_file),'w') as gvcf:
            for items in gvcfList:
                gvcf.write(working_dir + "/04_variantCalling/" + str(items) + "\n")

# Combine gvcf files for joint calling using "gatk CombineGVCFs"
rule combine_gvcf:
    input:
        gvcf_list_file = working_dir + "/05_jointCalling/gvcf.list"
    output:
        combined_gvcf = temp(working_dir + "/05_jointCalling/" + project + ".raw.snps.indels.g.vcf"),
        combined_gvcf_index = temp(working_dir + "/05_jointCalling/" + project + ".raw.snps.indels.g.vcf.idx")
    run:
        shell("gatk CombineGVCFs --java-options '-Xmx128G' -R {species_index} -V {input.gvcf_list_file} -O {output.combined_gvcf}")

# Do joint calling using "gatk GenotypeGVCFs"
rule gvcf_joint_calling:
    input:
        combined_gvcf = working_dir + "/05_jointCalling/" + project + ".raw.snps.indels.g.vcf"
    output:
        raw_joint_vcf = working_dir + "/05_jointCalling/" + project + ".raw.snps.indels.vcf"
    run:
        shell("gatk GenotypeGVCFs --java-options '-Xmx128G' -R {species_index} -V {input.combined_gvcf} -O {output.raw_joint_vcf}")

# Select and separate SNPs and INDELs into different files
rule separate_snp_indel:
    input:
        raw_joint_vcf = working_dir + "/05_jointCalling/" + project + ".raw.snps.indels.vcf"
    output:
        raw_snp_vcf = working_dir + "/05_jointCalling/SNP_analysis/" + project + ".raw.snps.vcf",
        raw_indel_vcf = working_dir + "/05_jointCalling/INDEL_analysis/" + project + ".raw.indels.vcf"
    run:
        shell("gatk SelectVariants --java-options '-Xmx128G' -R {species_index} -V {input.raw_joint_vcf} -O {output.raw_snp_vcf} --select-type-to-include SNP")
        shell("gatk SelectVariants --java-options '-Xmx128G' -R {species_index} -V {input.raw_joint_vcf} -O {output.raw_indel_vcf} --select-type-to-include INDEL")

# Do hard filtering to SNP and indels using different criteria
rule hard_filtering_snp_indel:
    input:
        raw_snp_vcf = working_dir + "/05_jointCalling/SNP_analysis/" + project + ".raw.snps.vcf",
        raw_indel_vcf = working_dir + "/05_jointCalling/INDEL_analysis/" + project + ".raw.indels.vcf"
    output:
        filtered_snp_vcf = working_dir + "/05_jointCalling/SNP_analysis/" + project + ".filtered.snps.vcf",
        filtered_indel_vcf = working_dir + "/05_jointCalling/INDEL_analysis/" + project + ".filtered.indels.vcf"
    run:
        shell("gatk VariantFiltration -R {species_index} -V {input.raw_snp_vcf} -O {output.filtered_snp_vcf} -filter-name 'QD_filter' -filter 'QD < 2.0' -filter-name 'FS_filter' -filter 'FS > 60.0' -filter-name 'MQ_filter' -filter 'MQ < 40.0' -filter-name 'SOR_filter' -filter 'SOR > 4.0' -filter-name 'MQRankSum_filter' -filter 'MQRankSum < -12.5' -filter-name 'ReadPosRankSum_filter' -filter 'ReadPosRankSum < -8.0'")
        shell("gatk VariantFiltration -R {species_index} -V {input.raw_indel_vcf} -O {output.filtered_indel_vcf} -filter-name 'QD_filter' -filter 'QD < 2.0' -filter-name 'FS_filter' -filter 'FS > 200.0' -filter-name 'SOR_filter' -filter 'SOR > 10.0'")

# Filter out variants that PASS all filters
rule filter_PASS_sites:
    input:
        filtered_snp_vcf = working_dir + "/05_jointCalling/SNP_analysis/" + project + ".filtered.snps.vcf",
        filtered_indel_vcf = working_dir + "/05_jointCalling/INDEL_analysis/" + project + ".filtered.indels.vcf"
    output:
        filtered_PASS_snp_vcf = working_dir + "/05_jointCalling/SNP_analysis/" + project + ".filtered.PASS.snps.vcf",
        filtered_PASS_indel_vcf = working_dir + "/05_jointCalling/INDEL_analysis/" + project + ".filtered.PASS.indels.vcf"
    run:
        shell("gatk SelectVariants --exclude-filtered -V {input.filtered_snp_vcf} -O {output.filtered_PASS_snp_vcf}")
        shell("gatk SelectVariants --exclude-filtered -V {input.filtered_indel_vcf} -O {output.filtered_PASS_indel_vcf}")

# Do unphasing to SNP vcf file
rule unphase_snp_vcf:
    input:
        filtered_PASS_snp_vcf = working_dir + "/05_jointCalling/SNP_analysis/" + project + ".filtered.PASS.snps.vcf"
    output:
        filtered_PASS_unphased_snp_vcf = working_dir + "/05_jointCalling/SNP_analysis/" + project + ".filtered.PASS.unphased.snps.vcf"
    run:
        shell("whatshap unphase {input.filtered_PASS_snp_vcf} > {output.filtered_PASS_unphased_snp_vcf}")