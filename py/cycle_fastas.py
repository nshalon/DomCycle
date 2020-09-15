import argparse
import util


###############################################################################################
# Purpose: input adjacency list, paired mapped reads, kmer size, contig table and             #
# outputs variable edge coverage                                                               #
###############################################################################################

parser = argparse.ArgumentParser(description='write the fasta sequences of putative cycles')
parser.add_argument('ifn_cyc_contig_table', metavar='<table of contigs in each cycle>', type=str, help='cycle_contig_table')
parser.add_argument('ifn_renamed_contig_fasta', metavar='<contig fasta>', type=str, help='contig fasta')
parser.add_argument('ifn_read_stats', metavar='<kmer size>', type=str, help='kmer size used for assembly')
parser.add_argument('ofn_cyc_fasta', metavar='<output cycle fastas>', type=str, help='output cycle fastas')

args = parser.parse_args()

k = int(util.get_k(args.ifn_read_stats, True))
print("K:", k)
contigname_to_seq = util.parse_fasta(args.ifn_renamed_contig_fasta)

ofile = open(args.ofn_cyc_fasta,"w+")
print("Writing fastas...")

with open(args.ifn_cyc_contig_table) as path:
    header = util.get_header_dict(args.ifn_cyc_contig_table)
    next(path)
    cycle = ""
    cycle_sequence = ""
    for line in path:
        line = util.split(line)
        curr_cycle, contig, orientation, fastg_edge = line[header["cycle"]], line[header["original_contig"]], line[header["orientation"]], bool(line[header["FASTG_edge"]])
        if cycle != curr_cycle:
            if cycle != "":
                ofile.write(cycle_sequence + "\n")
            cycle_sequence = ""
            cycle = curr_cycle
            ofile.write(">" + cycle + "\n")
        contig_sequence = contigname_to_seq[contig]
        if orientation == "-":
            contig_sequence = util.reverse_complement(contig_sequence)
        if fastg_edge:
            cycle_sequence += contig_sequence[:-k]
        else:
            cycle_sequence += contig_sequence
    ofile.write(cycle_sequence)





