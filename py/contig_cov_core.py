import argparse
import util

#############################################################################################
# Purpose: input the contig fasta, and two unpaired, filtered map files mapped using        #
# bwa mem to return a file with each contig's coverage which is defined as the amount of    #
# reads that map to that contig divided by the length                                       #
##############################################################################################

parser = argparse.ArgumentParser(description='Get contig coverages from parsed bwa map file')
parser.add_argument('ifn_fasta', metavar='<renamed contig fasta>', type=str, help='renamed contig fasta')
parser.add_argument('ifn_bwa_R1_table', metavar='<bwa R1 parsed filtered table>', type=str, help='perl parsed bwa table for R1')
parser.add_argument('ifn_bwa_R2_table', metavar='<bwa R2 parsed filtered table>', type=str, help='perl parsed bwa table for R2')
parser.add_argument('read_length', metavar='avg read length', type=float, help='size of the reads')
parser.add_argument('ofn_coverages', metavar='<output path for file that gives contig coverages>', type=str, help='contig coverages output')

args = parser.parse_args()

read_length = args.read_length
bwa_head = util.get_header_dict(args.ifn_bwa_R1_table)

contig_to_readc = {}
contig_to_length = {}

with open(args.ifn_fasta) as fasta:
    contig = length = None
    for line in fasta:
        line = line.rstrip()
        if line[0] == ">":
            contig = line[1:]
        else:
            length = int(len(line))
            contig_to_readc[contig] = 0
            contig_to_length[contig] = length

print("Going over R1")

with open(args.ifn_bwa_R1_table) as R1:
    next(R1)
    for line in R1:
        line = line.rstrip().split("\t")
        contig = line[bwa_head['contig']]
        contig_to_readc[contig] += 1

print("Going over R2")

with open(args.ifn_bwa_R2_table) as R2:
    next(R2)
    for line in R2:
        line = line.rstrip().split("\t")
        contig = line[bwa_head['contig']]
        contig_to_readc[contig] += 1

ofile = open(args.ofn_coverages,"w+")
ofile.write("contig\tcov\tlength\n")
for contig in sorted(contig_to_readc.keys()):
    length = contig_to_length[contig]
    cov = round(float(contig_to_readc[contig]) * read_length / length, 2)
    ofile.write(contig + "\t" + str(cov) + "\t" + str(length) + "\n")