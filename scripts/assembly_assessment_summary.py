import argparse
import sys

# Input main functions
if len(sys.argv) > 1 and sys.argv[1] == 'quast':
    mode = 'quast'
elif len(sys.argv) > 1 and sys.argv[1] == 'BUSCO':
    mode = 'BUSCO'
else:
    sys.stderr.write('This is a supplementary python script for summarizing NGS assembly statistic.\n')
    sys.stderr.write('usage: python main.py <mode> --help \n')
    sys.stderr.write('modes: quast, BUSCO\n')
    sys.exit(1)

if mode == 'quast':
    parser = argparse.ArgumentParser(description='summarizing quast outputs',usage='python assembly_assessment_summary.py quast -r <raw NGS contigs QUAST report.tsv> -f <filtered NGS contigs QUAST report.tsv> -o <outputfile>')
    parser.add_argument('quast')
    parser.add_argument("-r","--raw",required=True,type=str,dest='raw',help='Input the report.tsv file of raw NGS contigs obtained by QUAST.')
    parser.add_argument("-f","--filtered",required=True,type=str,dest='filtered',help='Input the report.tsv file of filtered NGS contigs obtained by QUAST.')
    parser.add_argument("-o","--output",required=True,type=str,dest='output',help='Output file name')
    args = parser.parse_args()
    if args.raw:
        raw_report = str(args.raw)
    if args.filtered:
        filtered_report = str(args.filtered)
    if args.output:
        outputFile = str(args.output)
# if mode == 'BUSCO':
#     parser = argparse.ArgumentParser(description='summarizing BUSCO outputs',usage='python assembly_assessment_summary.py quast -r <raw NGS contigs QUAST report.tsv> -f <filtered NGS contigs QUAST report.tsv> -o <outputfile>')

# Obtain assembly assessment results from QUAST report.tsv files
with open(outputFile,'a+') as outputSummary:
    with open(raw_report) as raw_assessment:
        for lines in raw_assessment.readlines():
            if lines.startswith("Assembly"):
                assembly_name = lines.split("\t")[1].split("_")[0].rstrip()
            elif lines.startswith("# contigs (>= 0 bp)"):
                raw_contig_num = lines.split("\t")[1].rstrip()
            elif lines.startswith("Total length (>= 0 bp)"):
                raw_contig_length = lines.split("\t")[1].rstrip()
            elif lines.startswith("GC (%)"):
                raw_contig_GC = lines.split("\t")[1].rstrip()
            elif lines.startswith("N50"):
                raw_contig_n50 = lines.split("\t")[1].rstrip()
            elif lines.startswith("NG50"):
                raw_contig_ng50 = lines.split("\t")[1].rstrip()
    with open(filtered_report) as filtered_assessment:
        for lines in filtered_assessment.readlines():
            if lines.startswith("# contigs (>= 0 bp)"):
                filtered_contig_num = lines.split("\t")[1].rstrip()
            elif lines.startswith("Total length (>= 0 bp)"):
                filtered_contig_length = lines.split("\t")[1].rstrip()
            elif lines.startswith("GC (%)"):
                filtered_contig_GC = lines.split("\t")[1].rstrip()
            elif lines.startswith("N50"):
                filtered_contig_n50 = lines.split("\t")[1].rstrip()
            elif lines.startswith("NG50"):
                filtered_contig_ng50 = lines.split("\t")[1].rstrip()
    outputSummary.write(str(assembly_name)+"\t" + 
    str(raw_contig_num)+"\t"+str(raw_contig_length)+"\t"+str(raw_contig_GC)+"\t"+str(raw_contig_n50)+"\t"+str(raw_contig_ng50)+"\t"+
    str(filtered_contig_num)+"\t"+str(filtered_contig_length)+"\t"+str(filtered_contig_GC)+"\t"+str(filtered_contig_n50)+"\t"+str(filtered_contig_ng50))