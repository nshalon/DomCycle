import argparse
import util


###############################################################################################
# Purpose: input adjacency list, paired mapped reads, kmer size, contig table and             #
# outputs variable edge coverage                                                               #
###############################################################################################

parser = argparse.ArgumentParser(description='write the fasta sequences of putative cycles')
parser.add_argument('ifn_R1', metavar='<filtered R1.fastq>', type=str, help='min contig size')
parser.add_argument('ifn_R2', metavar='<filtered R2.fastq>', type=str, help='contig fasta input')
parser.add_argument('ofn_singletons', metavar='<unpaired reads>', type=str, help='contig fasta input')
parser.add_argument('ofn_paired', metavar='<contig fasta input>', type=str, help='contig fasta input')

args = parser.parse_args()

header = util.split(open(args.ifn_R1).readline())

print("Going over R1")

singletons = set()
id_to_stats = {}
with open(args.ifn_R1) as R1:
    next(R1)
    dup = {}
    for line in R1:
        id = util.split(line)[0]
        id_to_stats[id] = line.rstrip()
        singletons.add(id)

print("Going over R2")

opaired = open(args.ofn_paired,"w+")
opaired.write(util.write_line("id","contig1", "coord1", "back_coord1", "strand1", "edit_dist1", "score1", "match_length1", "cigar1", "substitute1", "insert1", "delete1", "clip1", "unique1","contig2", "coord2", "back_coord2", "strand2", "edit_dist2", "score2", "match_length2", "cigar2", "substitute2", "insert2", "delete2", "clip2", "unique2"))
with open(args.ifn_R2) as R2:
    next(R2)
    for line in R2:
        id = util.split(line)[0]
        if id_to_stats.get(id) != None:
            if id in singletons:
                singletons.remove(id)
            opaired.write(id_to_stats[id])
            for element in util.split(line)[1:]:
                opaired.write("\t" + element)
            opaired.write("\n")
        else:
            singletons.add(id)
            id_to_stats[id] = line.rstrip()

print("Writing singletons...")

osingletons = open(args.ofn_singletons, "w+")
osingletons.write(util.write_line("id", "contig", "coord", "back_coord", "strand", "edit_dist", "score", "match_length", "cigar", "substitute", "insert", "delete", "clip", "unique"))
for id in singletons:
    osingletons.write(util.write_line(id_to_stats[id]))

print("Percent reads with a pair", round(100*(len(id_to_stats.keys()) - len(singletons)) / len(id_to_stats.keys()),2))
print("Number of reads without a pair", len(singletons))