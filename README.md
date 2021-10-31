# FYRAflow: Fission Yeast Re-sequencing Analysis workflow
- - - 
## Description
### Guo-Song Jia@Li-Lin Du's lab
FYRAflow is a `snakemake` based workflow for fission yeast species re-sequencing data analysis. \
Overview of the `FYRAflow variant` submodule:\
<!-- <img scr="https://github.com/guosongjia/Private_scripts/blob/master/FYRAflow_flowchart_new.jpg" width=300> -->
<!-- ![image](https://github.com/guosongjia/Private_scripts/blob/master/FYRAflow_flowchart_new.jpg){:height="50%" width="50%"} -->
<div align=center><img src="https://github.com/guosongjia/Private_scripts/blob/master/FYRA_flow_flow-chart_21.10.31.png"/></div>
## Full List of Tools used in this pipeline
- `Fastp` 
- `Kat`
- `BWA`
- `Samtools`
- `Picard Tools` (included in GATK4 packages)
- `GATK4`
- `Vcftools`
- `QUAST`
- `Whatshap`
- `R` (Dependency for some GATK steps)
## Quick Start
### Installation
#### Clone the repository:
`git clone https://github.com/guosongjia/FYRAflow.git`
#### Create the environment:
`conda env create -n fyraflow -f env.yaml`
#### Activate the environment:
`conda activate fyraflow`
### Set up configuration file
#### Customize the workflow based on your need in `config/*.yaml`
#### Using `config/variant_config.yaml` as an example 
```
# Please check the parameters, and adjust them according to your circumstance

# ================== Parameters for fyraflow variant sub-modules ==================

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
This is the main python script for Fission Yeast Re-sequencing Analysis workflow.
usage: python main.py <mode> --help 
modes: variant, assembly

# Run "variant" submodule as an example
python main.py variant -h
usage: python main.py variant -s <running step> -l <PATH to store log files>

fyraflow variant mode options

positional arguments:
  variant

optional arguments:
  -h, --help            show this help message and exit
  -s STEP, --step STEP  Input the step the workflow: qualityControl,
                        kmerFiltering, readsMapping_duplicatesMark,
                        gvcf_calling, jointCallingFiltering
  -l LOG, --log LOG     Input the path to save the running time log files.
```

