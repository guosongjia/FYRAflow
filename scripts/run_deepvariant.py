##############################################################
#  script: run_deepvariant.py
#  author: Guo-Song Jia
#  last edited: 2021.07.13
#  description: Script for in FYRAflow workflow. Run Deepvariant calling to TGS assembly for accuracy assessment.
##############################################################

import yaml
import os
import csv
with open("config/assembly_config.yaml") as yamlFile:
    yamlLoader = yaml.safe_load(yamlFile)
assembly_path = yamlLoader['RAWassemblyPATH']
working_dir = yamlLoader['WORKPATH']
with open(yamlLoader['SAMPLES']) as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = list()
    appendixList = list()
    for row in reader:
        sampleList.append(row[0])
        appendixList.append(row[1])
        sampleList = list(set(sampleList))


os.system('mkdir {0}/Assessment/01C_correctness/01.2_Deepvariant'.format(working_dir))
saving_PATH = working_dir + "/Assessment/01C_correctness/01.2_Deepvariant"
for sample in sampleList:
    for appendix in appendixList:
        sampleName=str(sample) + "_" + str(appendix)
        os.system('docker run -v "{0}/":"/input" -v "{1}/":"/output" \
                google/deepvariant:"1.1.0" /opt/deepvariant/bin/run_deepvariant --model_type=WGS \
                --ref=/input/{2}.fasta --reads=/input/{2}.kat_filtered.bwa.sorted.rmdup.RGadded.bam \
                --output_vcf=/output/{2}.deepvariant.vcf'.format(assembly_path,saving_PATH,sampleName))
        os.system('sudo chown jiaguosong:yeast {0}/*'.format(saving_PATH))
        os.system('bash scripts/reliable_variant_select.sh deepvariant {0}/{1}.deepvariant.vcf {0} {1}'.format(saving_PATH,sampleName))
