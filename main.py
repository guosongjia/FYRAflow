# This is the main page of FYRAflow
import yaml
import sys
import os
import time
import argparse

# Input main functions
if len(sys.argv) > 1 and sys.argv[1] == 'variant':
    mode = 'variant'
elif len(sys.argv) > 1 and sys.argv[1] == 'assembly':
    mode = 'assembly'
else:
    sys.stderr.write('This is the main python script for Fission Yeast Re-sequencing Analysis workflow.\n')
    sys.stderr.write('usage: python main.py <mode> --help \n')
    sys.stderr.write('modes: variant, assembly\n')
    sys.exit(1)

# 'variant' mode argparser setting
if mode == 'variant':
    parser = argparse.ArgumentParser(description='fyraflow variant mode options',usage='python main.py variant -s <running step> -l <PATH to store log files>')
    parser.add_argument('variant')
    parser.add_argument("-s","--step",required=True,type=str,dest='step',help='Input the step the workflow: qualityControl, kmerFiltering, readsMapping_duplicatesMark, gvcf_calling, jointCallingFiltering')
    parser.add_argument("-l","--log",required=True,type=str,dest='log',help='Input the path to save the running time log files.')
    args = parser.parse_args()
    if args.step:
        Runstep = str(args.step)
    if args.log:
        TimeLog = str(args.log)

# 'assembly' mode argparser setting
elif mode == 'assembly':
    parser = argparse.ArgumentParser(description='fyraflow assemble mode options',usage='python main.py assemble -s <running step> -l <PATH to store log files>')
    parser.add_argument('assembly')
    parser.add_argument("-s","--step",required=True,type=str,dest='step',help='Input the step the workflow: NGSfiltering, QUAST_assessment')
    parser.add_argument("-l","--log",required=True,type=str,dest='log',help='Input the path to save the running time log files.')
    args = parser.parse_args()
    if args.step:
        Runstep = str(args.step)
    if args.log:
        TimeLog = str(args.log)

# Define running time calculator.
def spend_time(start_time, end_time):
    seconds = end_time - start_time
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    return "%d:%02d:%02d" % (hours, minutes, seconds)


# "Running fyraflow variant sub-module"
if mode == 'variant':
# Import global config file from given PATH
    with open("config/variant_config.yaml") as yamlFile:
        config = yaml.safe_load(yamlFile)
        project = config["PROJECT"]
        species = config["SPECIES"]
    # Start FYRAflow on the project
    print("Start FYRAflow variant sub-module on project: " + project)
    # Do qualityControl when needed
    if Runstep == "qualityControl":
        file_log_time = open(str(TimeLog) + "01_qualityControl_running_time.txt", "a+")
        file_log_time.write("\nProject name: " + project + "\n")
        file_log_time.write("Running step: quality Control." + "\n")
        file_log_time.write("Start time: " + time.ctime() + "\n")
        print("Let's start QC!")
        start_time = time.time()
        os.system("snakemake -s workflow/qualityControl.smk --cores 32")
        end_time = time.time()
        file_log_time.write("Time of running quality Control: " + spend_time(start_time, end_time) + "\n")
        file_log_time.close()
        print("Quality control is done!")
        os._exit(0)
    # Do kmerFiltering when meeded
    if Runstep == "kmerFiltering":
        file_log_time = open(str(TimeLog) + "02_kmerFiltering_running_time.txt", "a+")
        file_log_time.write("\nProject name: " + project + "\n")
        file_log_time.write("Running step: kmer statistic and filtering." + "\n")
        file_log_time.write("Start time: " + time.ctime() + "\n")
        print("Let's start KAT kmer statistic and filtering!")
        start_time = time.time()
        os.system("snakemake -s workflow/katFiltering.smk --cores 16")
        end_time = time.time()
        file_log_time.write("Time of running kmer statistic and filtering: " + spend_time(start_time, end_time) + "\n")
        file_log_time.close()
        print("KAT kmer statistic and filtering is done!")
    # Do readsMapping when meeded
    if Runstep == "readsMapping_duplicatesMark":
        file_log_time = open(str(TimeLog) + "03_readsMapping_running_time.txt", "a+")
        file_log_time.write("\nProject name: " + project + "\n")
        file_log_time.write("Running step: reads mapping, bam sorting and duplicates reads mark." + "\n")
        file_log_time.write("Start time: " + time.ctime() + "\n")
        print("Let's start reads mapping and sorting!")
        start_time = time.time()
        os.system("snakemake -s workflow/readsMapping_duplicateMark.smk --cores 32")
        end_time = time.time()
        file_log_time.write("Time of running reads mapping and bam sorting: " + spend_time(start_time, end_time) + "\n")
        file_log_time.close()
        print("BWA MEM reads mapping, samtools sorting and gatk markduplicates is done!")
    # Do gvcf file calling when needed
    # This step is time comsuming!
    if Runstep == "gvcf_calling":
        file_log_time = open(str(TimeLog) + "04_gvcfCalling_running_time.txt", "a+")
        file_log_time.write("\nProject name: " + project + "\n")
        file_log_time.write("Running step: gvcf file calling." + "\n")
        file_log_time.write("Start time: " + time.ctime() + "\n")
        print("Let's start gvcf file calling!")
        start_time = time.time()
        os.system("snakemake -s workflow/gvcfCalling.smk --cores 54")
        end_time = time.time()
        file_log_time.write("Time of running gvcf file calling: " + spend_time(start_time, end_time) + "\n")
        file_log_time.close()
        print("gatk haplotypecaller calling is done!")
    # Do joint calling and variant filtering when needed
    if Runstep == "jointCallingFiltering":
        file_log_time = open(str(TimeLog) + "05_jointCallingFiltering_running_time.txt", "a+")
        file_log_time.write("\nProject name: " + project + "\n")
        file_log_time.write("Running step: joint calling and variant filtering." + "\n")
        file_log_time.write("Start time: " + time.ctime() + "\n")
        print("Let's start joint calling and variant filtering!")
        start_time = time.time()
        os.system("snakemake -s workflow/jointCallingFiltering.smk --cores 54")
        end_time = time.time()
        file_log_time.write("Time of running joint calling and filtering: " + spend_time(start_time, end_time) + "\n")
        file_log_time.close()
        print("Joint calling and variant filtering is done!")

# "Running fyraflow assembly sub-module"
if mode == 'assembly':
    # Import global config file from given PATH
    with open("config/assembly_config.yaml") as yamlFile:
        config = yaml.safe_load(yamlFile)
        project = config["PROJECT"]
        species = config["SPECIES"]
        working_dir = config["WORKPATH"]
    # Start FYRAflow on the project
    print("Start FYRAflow assembly sub-module on project: " + project)
    if Runstep == "NGSfiltering":
        file_log_time = open(str(TimeLog) + "01_NGScontigFiltering_running_time.txt", "a+")
        file_log_time.write("\nProject name: " + project + "\n")
        file_log_time.write("Running step: NGS assembly filtering based on KAT filtered reads mapping." + "\n")
        file_log_time.write("Start time: " + time.ctime() + "\n")
        print("Let's start NGS assembly filtering based on KAT filtered reads mapping!")
        start_time = time.time()
        os.system("snakemake -s workflow/NGScontig_katReadsFiltering.smk --cores 48")
        end_time = time.time()
        file_log_time.write("Time of running NGS assembly filtering based on KAT filtered reads mapping: " + spend_time(start_time, end_time) + "\n")
        file_log_time.close()
        print("NGS assembly filtering based on KAT filtered reads mapping is done!")
    if Runstep == "QUAST_assessment":
        file_log_time = open(str(TimeLog) + "02_NGScontigQUASTassessment_running_time.txt", "a+")
        file_log_time.write("\nProject name: " + project + "\n")
        file_log_time.write("Running step: NGS assembly QUAST assessment." + "\n")
        file_log_time.write("Start time: " + time.ctime() + "\n")
        print("Let's start NGS assembly QUAST assessment!")
        start_time = time.time()
        os.system("snakemake -s workflow/NGScontig_assessment.smk --cores 48")
        os.system("cat {0}/02_assessment/*.txt >> {0}/02_assessment/{1}_NGS_contigs.QUAST_assessment".format(working_dir,species))
        os.system("rm {0}/02_assessment/*.txt ".format(working_dir))
        os.system("mv {0}/02_assessment/{1}_NGS_contigs.QUAST_assessment {0}/02_assessment/{1}_NGS_contigs.QUAST_assessment.txt".format(working_dir,species))
        end_time = time.time()
        file_log_time.write("Time of running NGS assembly QUAST assessment: " + spend_time(start_time, end_time) + "\n")
        file_log_time.close()
        print("NGS assembly QUAST assessment is done!")
    if Runstep == "TGScorrectness_assessment":
        file_log_time = open(str(TimeLog) + "TGScontigBaseLevelAssessment_running_time.txt","a+")
        file_log_time.write("\nProject name: " + project + "\n")
        file_log_time.write("Running step: TGS assembly base level correctness assessment." + "\n")
        file_log_time.write("Start time: " + time.ctime() + "\n")
        print("Let's start TGS assembly base level correctness assessment!")
        start_time = time.time()
        os.system("snakemake -s workflow/TGScontig_assessment.smk --cores 48")
        os.system("python scripts/run_deepvariant.py")
        end_time = time.time()
        file_log_time.write("Time of running TGS assembly base level correctness by NGS variant calling: " + spend_time(start_time, end_time) + "\n")
        file_log_time.close()
        print("TGS assembly base level correctness assessment is done!")
