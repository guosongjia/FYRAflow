import csv
configfile: "config/main_config.yaml"

working_dir = config["WORKPATH"]
species = config["SPECIES"]
## Generate reference genome index files
species_index = "genome/" + species + "/"+ species +".fasta"

with open(str(config["SAMPLE"]),'r') as sampleFile:
    reader = csv.reader(sampleFile)
    sampleList = [row[0] for row in reader]

