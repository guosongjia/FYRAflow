# FYRAflow: Fission Yeast Re-sequencing Analysis workflow
- - - 
## Description
### Guo-Song Jia@Li-Lin Du's lab
FYRAflow is a `snakemake` based workflow for fission yeast species re-sequencing data analysis. \
Over view of the workflow:\
![image](https://github.com/guosongjia/Private_scripts/blob/master/FYRAflow_flowchart.png)
## Full List of Tools used in this pipeline
- `Fastp` 
- `Kat`
- `BWA`
- `Samtools`
- `Picard Tools` (included in GATK4 packages)
- `GATK4`
- `Vcftools`
- `Whatshap`
- `R` (Dependency for some GATK steps)
## Quick Start
### Installation
#### Clone the repository:
`git clone https://github.com/guosongjia/FYRAflow.git`
#### Create the environment:
`conda env create -n fyraflow -f env.ymal`
#### Activate the environment:
`conda activate fyraflow`
### Set up configuration file
#### Customize the workflow based on your need in `config/main_config.yaml` 
```
# Please check the parameters, and adjust them according to your circumstance

# ================== Shared parameters for some or all of the sub-workflows ==================

## Project Name
PROJECT: test

## Fission yeast species of re-sequencing analysis
SPECIES: S_octosporus # "S_pombe" or "S_octosporus" or "S_japonicus" or "S_versatilis" or "S_osmophilus" or "S_cryophilus"

## Input the path of fastq files
READSPATH: example/raw_reads

## Input the working directory
WORKPATH: example/

## Input the sample list file
SAMPLE: example/sample.list

## Input the final output directory
OUTPUTPATH: example/
```
### Run FYRAflow:
```
python main.py -h
optional arguments:
  -h, --help           show this help message and exit
  -s STEP, -step STEP  Input the step the workflow: qualityControl, kmerFiltering, readsMapping_duplicatesMark, gvcf_calling
  -l LOG, -log LOG     Input the path to save the running time log files.

# Run "qualityControl" step as an example
python main.py -s qualityControl -l ./ 
```

