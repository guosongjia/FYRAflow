# This is the main page of FYRAflow
import yaml
import os
import time
import argparse

# argparser setting
parser = argparse.ArgumentParser(description='This is the main python script for Fission Yeast Re-sequencing Analysis workflow.')
parser.add_argument("-s","-step",required=True,type=str,dest='step',help='Input the step the workflow: qualityControl, kmerFiltering, readsMapping, duplicatesMark')
parser.add_argument("-l","-log",required=True,type=str,dest='log',help='Input the path to save the running time log files.')
args = parser.parse_args()
if args.step:
    Runstep = str(args.step)
if args.log:
    TimeLog = str(args.log)

# Import global config file from given PATH
with open("config/main_config.yaml") as yamlFile:
    config = yaml.safe_load(yamlFile)

project = config["PROJECT"]
species = config["SPECIES"]

# Start FYRAflow on the project
print("Start FYRAflow on project: " + project)

def spend_time(start_time, end_time):
    seconds = end_time - start_time
    hours = seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    
    return "%d:%02d:%02d" % (hours, minutes, seconds)

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
    os.system("snakemake -s workflow/katFiltering.smk --cores 32")
    end_time = time.time()
    file_log_time.write("Time of running kmer statistic and filtering: " + spend_time(start_time, end_time) + "\n")
    file_log_time.close()
    print("KAT kmer statistic and filtering is done!")

# Do readsMapping when meeded
if Runstep == "readsMapping":
    file_log_time = open(str(TimeLog) + "03_readsMapping_running_time.txt", "a+")
    file_log_time.write("\nProject name: " + project + "\n")
    file_log_time.write("Running step: reads mapping and bam sorting." + "\n")
    file_log_time.write("Start time: " + time.ctime() + "\n")
    print("Let's start reads mapping and sorting!")
    start_time = time.time()
    os.system("snakemake -s workflow/readsMapping.smk --cores 32")
    end_time = time.time()
    file_log_time.write("Time of running reads mapping and bam sorting: " + spend_time(start_time, end_time) + "\n")
    file_log_time.close()
    print("BWA MEM reads mapping, samtools sorting and gatk markduplicates is done!")



