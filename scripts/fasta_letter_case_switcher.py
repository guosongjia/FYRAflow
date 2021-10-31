import argparse
from Bio import SeqIO
##############################################################
#  script: fasta_letter_case_switcher.py
#  author: Guo-Song Jia
#  last edited: 2021.10.25
#  description: Script for in FYRAflow workflow. Switch the letter case of FASTA file between uppercase and lowercase.
##############################################################

# Argparse settings
parser = argparse.ArgumentParser(
    description='This is a build-in script in FYRAflow workflow. It can switch the letter case of FASTA file between uppercase and lowercase')
parser.add_argument("-i", "--input", required=True, dest='inputFasta',
                    help='Input fasta file for letter case transformation.')
parser.add_argument("-o", "--output", required=True, dest='outputFasta',
                    help='Output fasta file name after letter case transformation')
parser.add_argument("-c", "--case", required=True, dest='case',
                    help='Type the letter case for transformation: upper, lower')
args = parser.parse_args()
if args.case:
    caseTransfer = str(args.case)
if args.inputFasta:
    inputFile = SeqIO.parse(args.inputFasta, "fasta")
if args.outputFasta:
    outputFile = args.outputFasta


if caseTransfer == "upper":
    with open(outputFile, 'w') as Output:
        for seq_record in inputFile:
            Output.write(">" + str(seq_record.id) + "\n")
            Output.write(str(seq_record.seq.upper()) + "\n")
elif caseTransfer == "lower":
    with open(outputFile, 'w') as Output:
        for seq_record in inputFile:
            Output.write(">" + str(seq_record.id) + "\n")
            Output.write(str(seq_record.seq.lower()) + "\n")
