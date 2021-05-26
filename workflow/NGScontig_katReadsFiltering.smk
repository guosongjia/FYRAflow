import csv
import os
import numpy as np

configfile: "config/assembly_config.yaml"

species = config['SPECIES']
assembly_path = config['RAWassemblyPATH']
katReads_path = config['KATREADS']
working_dir = config["WORKPATH"]

with open(str(config["SAMPLES"]),'r') as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = sorted([row[0] for row in reader])

rule all:
    input:
        expand(working_dir + "/01_katReads_filtering/{sample}.passed_contigs.txt",sample=sampleList),
        expand(working_dir + "/01_katReads_filtering/{sample}.failed_contigs.txt",sample=sampleList),
        expand(working_dir + "/01_katReads_filtering/filtered_contigs/{sample}_NGS.filtered.contigs.fasta",sample=sampleList)

# Generate contig windows files
rule generate_contig_windows:
    input:
        NGS_raw_contigs = assembly_path + "/{sample}_NGS.contigs.fasta"
    output:
        NGS_raw_contigs_windows = temp(working_dir + "/01_katReads_filtering/{sample}.1k.windows")
    run:
        with open(input.NGS_raw_contigs) as raw_assembly:
            seqDict = {}
            for lines in raw_assembly:
                if lines[0] == ">":
                    seqName = lines[1:].rstrip()
                    seqDict[seqName] = ''
                else:
                    seqDict[seqName] += lines.replace("\n","")
        with open(output.NGS_raw_contigs_windows,'w') as output:
            for contigs in seqDict.keys():
                length = len(seqDict[contigs])
                times = length // 1000
                for i in range(times):
                    if i == 0:
                        output.write(str(contigs) + "\t1\t1000\n")
                    elif i > 0:
                        output.write(str(contigs) + "\t" + str(i*1000 + 1 ) + "\t" + str((i+1)*1000) + "\n")
                output.write(str(contigs) + "\t" + str(times*1000+1) + "\t" + str(length) + "\n")

# Generate bwa index files and do KAT reads mapping 
rule bwa_mapping_sorting_indexing:
    input: 
        NGS_raw_contigs = assembly_path + "/{sample}_NGS.contigs.fasta",
        gziped_filtered_fastp_kat_reads1 = katReads_path + "/{sample}.in.R1.fq.gz",
        gziped_filtered_fastp_kat_reads2 = katReads_path + "/{sample}.in.R2.fq.gz"
    output: 
        NGS_raw_contigs_amb = temp(assembly_path + "/{sample}_NGS.contigs.fasta.amb"),
        NGS_raw_contigs_ann = temp(assembly_path + "/{sample}_NGS.contigs.fasta.ann"),
        NGS_raw_contigs_bwt = temp(assembly_path + "/{sample}_NGS.contigs.fasta.bwt"),
        NGS_raw_contigs_pac = temp(assembly_path + "/{sample}_NGS.contigs.fasta.pac"),
        NGS_raw_contigs_sa = temp(assembly_path + "/{sample}_NGS.contigs.fasta.sa"),
        mapped_sam_file = temp(working_dir + "/{sample}.kat_filtered.bwa.sam"),
        mapped_sorted_bam_file = working_dir + "/01_katReads_filtering/{sample}.kat_filtered.bwa.sorted.bam",
        mapped_sorted_bam_index_file = working_dir + "/01_katReads_filtering/{sample}.kat_filtered.bwa.sorted.bam.bai"
    threads: 32
    run: 
        shell("bwa index {input.NGS_raw_contigs}")
        shell("bwa mem -t {threads} {input.NGS_raw_contigs} {input.gziped_filtered_fastp_kat_reads1} {input.gziped_filtered_fastp_kat_reads2} > {output.mapped_sam_file}")
        shell("samtools sort -o {output.mapped_sorted_bam_file} {output.mapped_sam_file}")
        shell("samtools index {output.mapped_sorted_bam_file}")

# Calculate binning base coverage
rule binning_coverage_calculation:
    input: 
        NGS_raw_contigs_windows = working_dir + "/01_katReads_filtering/{sample}.1k.windows",
        mapped_sorted_bam_file = working_dir + "/01_katReads_filtering/{sample}.kat_filtered.bwa.sorted.bam"
    output:
        binning_coverage = working_dir + "/01_katReads_filtering/{sample}.bining.coverage.txt"
    run:
        with open(input.NGS_raw_contigs_windows) as biningFile:
            biningList = [i.rstrip() for i in biningFile.readlines()]
        with open(output.binning_coverage,'w') as output:
            for bins in biningList:
                length = int(bins.split("\t")[2]) - int(bins.split("\t")[1]) +1
                depthRange = str(bins.split("\t")[0]) + ":" + str(bins.split("\t")[1]) + "-" + str(bins.split("\t")[2])
                depth = os.popen("samtools depth -aa -r {0} {1}|cut -f 3".format(depthRange, input.mapped_sorted_bam_file))
                total_depth = sum([int(i.rstrip()) for i in depth.readlines()])
                if int(total_depth) != 0:
                    coverage = int(total_depth) // int(length)
                    output.write(str(depthRange) + "\t" + str(coverage) + "\n")
                else:
                    output.write(str(depthRange) + "\t0\n")

# Filter contigs based on coverage
rule contig_filtering:
    input: 
        NGS_raw_contigs = assembly_path + "/{sample}_NGS.contigs.fasta",
        binning_coverage = working_dir + "/01_katReads_filtering/{sample}.bining.coverage.txt"
    output:
        passed_contigs_list = working_dir + "/01_katReads_filtering/{sample}.passed_contigs.txt",
        failed_contigs_list = working_dir + "/01_katReads_filtering/{sample}.failed_contigs.txt",
        NGS_filtered_contigs = working_dir + "/01_katReads_filtering/filtered_contigs/{sample}_NGS.filtered.contigs.fasta"
    run: 
        with open(input.NGS_raw_contigs) as biningFile:
            fasta_reader = biningFile.readlines()
            contigList = [i.rstrip()[1:] for i in fasta_reader if i[0] == ">"]
        with open(input.binning_coverage) as coverage:
            passed = open(output.passed_contigs_list,'w')
            failed = open(output.failed_contigs_list,'w')
            reader = coverage.readlines()
            median_coverage = np.median([int(i.split("\t")[1].rstrip()) for i in reader])
            passedContigList = list()
            for contig in contigList:
                contigCoverage = [int(i.split("\t")[1].rstrip()) / median_coverage for i in reader if i.split(":")[0] == contig]
                contigNumber = int(len(contigCoverage))
                passedBins = [i for i in contigCoverage if i > 0.3]
                passedBinsNumber = int(len(passedBins))
                failedBins = [i for i in contigCoverage if i <= 0.3]
                failedBinsNumber = int(len(failedBins))
                if passedBinsNumber >= contigNumber / 2:
                    passed.write(str(contig) + "\n")
                    passedContigList.append(contig)
                elif passedBinsNumber < contigNumber / 2:
                    failed.write(str(contig) + "\n")
            passed.close()
            failed.close()
        with open(output.NGS_filtered_contigs,'w') as outputFasta:
            flag = 0
            for passedContigs in passedContigList:
                for lines in fasta_reader:
                    if lines[0] == ">":
                        if lines[1:].rstrip() == passedContigs:
                            flag = 1
                            outputFasta.write(">" + str(passedContigs) + "\n")
                        else:
                            flag = 0
                    else:
                        if flag == 0:
                            pass
                        else:
                            outputFasta.write(lines)

