# FYRAflow: Fission Yeast Re-sequencing Analysis workflow
- - - 
## Description
### Guo-Song Jia@Li-Lin Du's lab
FYRAflow is a `snakemake` based workflow for fission yeast species re-sequencing data analysis. \
The raw NGS reads data will be analyzed and trimmed using `fastp`. Then trimmed reads will be filtered based on kmer statistic and the low frequency kmer reads will be removed using `kat`.
## Quick start
### Installation
#### Clone the repository:
`git clone https://github.com/guosongjia/FYRAflow.git`
#### Create the environment:
`conda env create -n fyraflow -f env.ymal`
#### Activate the environment:
`conda activate fyraflow`
### Set up configuration file
#### Customize the workflow based on your need in `config/main_config.yaml` 

### Run FYRAflow:
```
python main.py -h 

```

