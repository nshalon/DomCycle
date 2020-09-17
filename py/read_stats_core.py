import argparse
import util
import numpy as np

parser = argparse.ArgumentParser(description='Get contig coverages from parsed bwa map file')
parser.add_argument('ifn_paired_table', metavar='paired read table', type=str, help='contig table w/ lengths')
parser.add_argument('k', metavar='assembly kmer', type=str, help='contig table w/ lengths',)
parser.add_argument('read_length', metavar='pre-trimmed read length', type=str, help='avg length of reads before any trimming')
parser.add_argument('read_stats', metavar='read stats file path', type=str, help='output read stats file path')

args = parser.parse_args()

insert_dist = util.get_insert(args.ifn_paired_table)
max_distance = np.percentile(sorted(insert_dist), 99.5) + 200
med_insert = np.median(insert_dist)
ofile = open(args.read_stats, "w+")
ofile.write(util.write_line("read_length", "k", "med_insert", "max_distance"))
ofile.write(util.write_line(args.read_length, args.k, med_insert, max_distance))