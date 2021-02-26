import csv
configfile: "config/main_config.yaml"

input_path = config["READSPATH"]
working_dir = config["WORKPATH"]
intermediate_dir = config["WORKPATH"] + "/intermediate_" + config["PROJECT"]
species = config["SPECIES"]
## Generate reference genome index files
species_index = "genome/" + species + "/"+ species +".fasta"

print(species_index)
with open(str(config["SAMPLE"]),'r') as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = [row[0] for row in reader]

rule all:
    input:
        expand(working_dir + "/03_readsMapping/{sample}.bwa.sorted.rmdup.bam",sample=sampleList),
        expand(working_dir + "/03_readsMapping/{sample}.bwa.sorted.rmdup.bam.bai",sample=sampleList),
        expand(working_dir + "/03_readsMapping/03.1_rmdup_metrics/{sample}_marked_dup_matrics.txt",sample=sampleList),
        expand(working_dir + "/03_readsMapping/03.2_insert_size/{sample}_insert_size_metrics.txt",sample=sampleList),
        expand(working_dir + "/03_readsMapping/03.2_insert_size/{sample}.insert_size_histogram.pdf",sample=sampleList)

## Do bwa mem reads mapping and obtain sam files. 
rule bwa_mapping:
    input:
        gziped_filtered_fastp_read1 = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R1.fq.gz",
        gziped_filtered_fastp_read2 = working_dir + "/02_kmerStatFiltering_kat/filtered_reads/{sample}.in.R2.fq.gz"
    output:
        mapped_sam_file = temp(intermediate_dir + "{sample}.bwa.sam")
    threads: 32
    run: 
        shell("bwa mem -t {threads} -R '@RG\tID:{sample}\tPL:ILLUMINA\tLB:{sample}\tSM:{sample}' {species_index} {input.gziped_filtered_fastp_read1} {input.gziped_filtered_fastp_read2} > {output.mapped_sam_file} ")

## Do samtools sorting and generate bam index files 
rule sam_to_sotredbam:
    input:
        mapped_sam_file = intermediate_dir + "{sample}.bwa.sam"
    output:
        mapped_sorted_bam_file = temp(working_dir + "/03_readsMapping/{sample}.bwa.sorted.bam"),
        mapped_sorted_bam_file_index = temp(working_dir + "/03_readsMapping/{sample}.bwa.sorted.bam.bai")
    run:
        shell("samtools sort {input.mapped_sam_file} -o {output.mapped_sorted_bam_file}")
        shell("samtools index {output.mapped_sorted_bam_file}")

## Do duplicates reads mark using gatk MarkDuplicates and generate bam index files
rule mark_duplicates:
    input:
        mapped_sorted_bam_file = working_dir + "/03_readsMapping/{sample}.bwa.sorted.bam"
    output:
        mapped_sorted_rmdup_bam_file = working_dir + "/03_readsMapping/{sample}.bwa.sorted.rmdup.bam",
        mapped_sorted_rmdup_bam_file_index = working_dir + "/03_readsMapping/{sample}.bwa.sorted.rmdup.bam.bai",
        rmdup_metrics = working_dir + "/03_readsMapping/03.1_rmdup_metrics/{sample}_marked_dup_matrics.txt"
    run:
        shell("gatk MarkDuplicates --java-options '-Xmx16G' -I {input.mapped_sorted_bam_file} -O {output.mapped_sorted_rmdup_bam_file} -M {output.rmdup_metrics}")
        shell("samtools index {output.mapped_sorted_rmdup_bam_file_index}")

## Collect insert size distribution and alignment summary using gatk tools
rule collect_metrics:
    input:
        mapped_sorted_rmdup_bam_file = working_dir + "/03_readsMapping/{sample}.bwa.sorted.rmdup.bam"
    output:
        insert_size_metrics = working_dir + "/03_readsMapping/03.2_insert_size/{sample}_insert_size_metrics.txt"
        insert_size_histogrom = working_dir + "/03_readsMapping/03.2_insert_size/{sample}.insert_size_histogram.pdf",
        alignemnt_metrics = working_dir + "/03_readsMapping/03.3_alignment_summary/{sample}.alignment_summary.txt"
    run:
        shell("gatk CollectInsertSizeMetrics --java-options '-Xmx16G' -I {input.mapped_sorted_rmdup_bam_file} -O {output.insert_size_metrics}")
        shell("gatk CollectAlignmentSummaryMetrics --java-options '-Xmx16G' -R {species_index} -I {input.mapped_sorted_rmdup_bam_file} -O ${output.alignemnt_metrics}")
        