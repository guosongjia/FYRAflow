import argparse
import os
from Bio.Seq import Seq
##############################################################
#  script: extract_sequence_from_tidy_gff3.py
#  author: Guo-Song Jia
#  last edited: 2021.12.12
#  description: Script for in FYRAflow workflow. Extract different types of sequences from the gff3 file after tidying-up.
##############################################################

# Argparse settings
parser = argparse.ArgumentParser(
    description='This is a build-in script in FYRAflow workflow. It can extract different types of sequences (gene, CDS, protein) from the gff3 file after tidying-up.')
parser.add_argument("-g", "--genome", required=True, dest='genomeFasta',
                    help='Input the genome fasta file.')
parser.add_argument("-t", "--type", required=True, dest='featureType',
                    help='Input the sequence type for sequence extraction.')
parser.add_argument("-a", "--annotation", required=True, dest='annotationGffFile',
                    help='Input the genome annotation gff3 file.')
parser.add_argument("-o", "--output", required=True, dest='outputFastaFile',
                    help='Output the fasta file containing sequence corresponding to.')
args = parser.parse_args()

# Read the coding gene name as list
def extractGeneName(annotationGffFile):
    with open(annotationGffFile) as readin:
        geneNameList = ["".join(i.split("\t")[8].split(";")[0].split("=")[1]) for i in readin.readlines() if not i.startswith("#") and i.split("\t")[2] == "gene" ]
    return geneNameList

# Extract gene sequence
def extractGeneSequence(genomeFasta,annotationGffFile,geneName):
    with open(annotationGffFile) as readin:
        outputVal = list()
        coordinateDict = {}
        coordinateDict["name"] = str(geneName)
        outputVal.append(">" + str(geneName))
        for i in readin.readlines():
            if not i.startswith("#") and i.split("\t")[2] == "gene" and str(geneName) == "".join(i.split("\t")[8].split(";")[0].split("=")[1]):
                coordinateDict["contig"] = str(i.split("\t")[0])
                coordinateDict["start"] = int(i.split("\t")[3])
                coordinateDict["end"] = int(i.split("\t")[4])
                coordinateDict["chain"] = str(i.split("\t")[6])
    if coordinateDict["chain"] == "+":
        output = os.popen('samtools faidx {0} {1}:{2}-{3}|grep -v ">" |tr -d "\n" |tr "atgcn" "ATGCN"'.format(genomeFasta,coordinateDict["contig"],coordinateDict["start"],coordinateDict["end"]))
        valRead = "".join([i for i in output.readlines()])
        outputVal.append(valRead + "\n")
    elif coordinateDict["chain"] == "-":
        output = os.popen('samtools faidx {0} {1}:{2}-{3}|grep -v ">"|tr -d "\n"| tr "[ATGCatgcNn]" "[TACGTACGNN]" | rev | tr -d "\n"'.format(genomeFasta,coordinateDict["contig"],coordinateDict["start"],coordinateDict["end"]))
        valRead = "".join([i for i in output.readlines()])
        outputVal.append(valRead + "\n")
    return outputVal

# Extract CDS sequence
def extractCdsSequence(genomeFasta,annotationGffFile,geneName):
    with open(annotationGffFile) as readin:
        outputVal = list()
        coordinateDict = {}
        coordinateDict["name"] = str(geneName)
        cdsPosList = list()
        for i in readin.readlines():
            if not i.startswith("#") and i.split("\t")[2] == "gene" and str(geneName) == "".join(i.split("\t")[8].split(";")[0].split("=")[1]):
                coordinateDict["contig"] = str(i.split("\t")[0])
                coordinateDict["chain"] = str(i.split("\t")[6])
            elif not i.startswith("#") and i.split("\t")[2] == "CDS" and str(geneName) == "".join(i.split("\t")[8].split(";")[0].split(".")[1]):
                cdsPosList.append(int(i.split("\t")[3]))
                cdsPosList.append(int(i.split("\t")[4]))
        cdsPosList.sort()
        cdsNumber = int(len(cdsPosList)/2)
        for num in range(cdsNumber):
            coordinateDict["cds-" + str(num+1)] = str(cdsPosList[num*2]) + "-" + str(cdsPosList[num*2+1])
        if coordinateDict["chain"] == "+":
            for num in range(cdsNumber):
                startPos = coordinateDict["cds-"+str(num+1)].split("-")[0]
                endPos = coordinateDict["cds-"+str(num+1)].split("-")[1]
                cdsRange = str(startPos) + "-" + str(endPos)
                output = os.popen('samtools faidx {0} {1}:{2}|grep -v ">" |tr -d "\n" |tr "atgcn" "ATGCN"'.format(genomeFasta,coordinateDict["contig"],cdsRange))
                valRead = "".join([i for i in output.readlines()])
                outputVal.append(valRead)
        elif coordinateDict["chain"] == "-":
            for num in range(cdsNumber):
                startPos = coordinateDict["cds-"+str(cdsNumber-num)].split("-")[0]
                endPos = coordinateDict["cds-"+str(cdsNumber-num)].split("-")[1]
                cdsRange = str(startPos) + "-" + str(endPos)
                output = os.popen('samtools faidx {0} {1}:{2}|grep -v ">" |tr -d "\n" |tr "[ATGCatgcNn]" "[TACGTACGNN]" | rev | tr -d "\n"'.format(genomeFasta,coordinateDict["contig"],cdsRange))
                valRead = "".join([i for i in output.readlines()])
                outputVal.append(valRead)
        sequence = "".join(outputVal)
        outputSequence = list()
        outputSequence.append(">" + str(geneName))
        outputSequence.append(sequence + "\n")
    return outputSequence

# Extract protein sequence
def translateCDStoProtein(geneCDSsequence):
    geneticCodeDict = {'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TCN':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L','TAC':'Y',
                       'TAT':'Y','TAA':'*','TAG':'*','TGC':'C','TGT':'C','TGA':'*','TGG':'W','CTA':'L','CTC':'L','CTG':'L',
                       'CTT':'L','CTN':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P','CCN':'P','CAC':'H','CAT':'H','CAA':'Q',
                       'CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R','CGN':'R','ATA':'I','ATC':'I','ATT':'I','ATG':'M',
                       'ACA':'T','ACC':'T','ACG':'T','ACT':'T','ACN':'T','AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S',
                       'AGT':'S','AGA':'R','AGG':'R','GTA':'V','GTC':'V','GTG':'V','GTT':'V','GTN':'V','GCA':'A','GCC':'A',
                       'GCG':'A','GCT':'A','GCN':'A','GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G',
                       'GGT':'G','GGN':'G','MGR':'R','AAY':'N','GAY':'D','TGY':'C','CAR':'Q','GAR':'E','CAY':'H','ATH':'I',
                       'YTR':'L','AAR':'K','TTY':'F','TCN':'S','AGY':'S','ACN':'T','TAY':'Y','TAR':'*','TRA':'*','---':'-'}
    triples = [geneCDSsequence[i:i+3] for i in range(0, len(geneCDSsequence), 3)]
    aminoacids = [geneticCodeDict[triple] for triple in triples]
    protein = "".join(aminoacids)
    return protein

# Generate gene list
geneList = extractGeneName(args.annotationGffFile)
if str(args.featureType) == "gene":
    with open(str(args.outputFastaFile),'w') as outputFile:
        for gene in geneList:
            outputFile.write("\n".join(extractGeneSequence(args.genomeFasta,args.annotationGffFile,str(gene))))
elif str(args.featureType) == "CDS":
    with open(str(args.outputFastaFile),'w') as outputFile:
        for gene in geneList:
            outputFile.write("\n".join(extractCdsSequence(args.genomeFasta,args.annotationGffFile,str(gene))))
elif str(args.featureType) == "protein":
    with open(str(args.outputFastaFile),'w') as outputFile:
        for gene in geneList:
            cdsSequence = Seq(extractCdsSequence(args.genomeFasta,args.annotationGffFile,str(gene))[1].replace('\n',''))
            proteinSequence = cdsSequence.translate()
            outputFile.write(">" + str(gene) + "\n")
            outputFile.write(str(proteinSequence) + "\n")