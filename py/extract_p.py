import argparse
import util
import itertools
from scipy.stats import binom

###############################################################################################
# Purpose: input adjacency list, paired mapped reads, kmer size, contig table and             #
# outputs variable edge coverage                                                               #
###############################################################################################

parser = argparse.ArgumentParser(description='write the fasta sequences of putative cycles')
parser.add_argument('ifn_cycle_sum', metavar='<minimum contig size to keep>', type=str, help='min contig size')
parser.add_argument('ifn_cycle_fasta', metavar='<contig fasta input>', type=str, help='contig fasta input')
parser.add_argument('ifn_cycle_contigs', metavar='<contig fasta input>', type=str, help='contig fasta input')
parser.add_argument('max_pval', metavar='<contig fasta input>', type=float, help='contig fasta input')
parser.add_argument('min_length', metavar='<contig fasta input>', type=int, help='contig fasta input')
parser.add_argument('min_lbc', metavar='<lower bound on coverage>', type=float, help='contig fasta input')
parser.add_argument('paired_weird_mult', metavar='<lower bound on coverage>', type=float, help='contig fasta input')
parser.add_argument('ofn_all_cycle_sum', metavar='<contig fasta input>', type=str, help='contig fasta input')
parser.add_argument('ofn_dominant_cycle_sum', metavar='<contig fasta input>', type=str, help='contig fasta input')
parser.add_argument('ofn_cycle_fasta', metavar='<minimum contig size to keep>', type=str, help='min contig size')
parser.add_argument('ofn_cycle_contigs', metavar='<contig fasta input>', type=str, help='contig fasta input')

args = parser.parse_args()

cyc_to_length, cyc_to_contig, contig_to_length = util.get_cyc_contig_stats(args.ifn_cycle_contigs)
cycle_scores = util.get_cycle_scores(args.ifn_cycle_sum)

collapsed_cycle = {}
cycles = list(cyc_to_length.keys())
for cyc1, cyc2 in itertools.combinations(cycles, 2):
    contigs1 = set(cyc_to_contig[cyc1])
    contigs2 = set(cyc_to_contig[cyc2])
    intersect_contigs = list(set(contigs1) & set(contigs2))
    union_contigs = list(set().union(contigs1, contigs2))
    intersect_sum = sum([contig_to_length[contig] for contig in intersect_contigs])
    union_sum = sum([contig_to_length[contig] for contig in union_contigs])
    if intersect_sum / float(union_sum) >= 0.95:
        sorted_cyc_scores = sorted([(cyc1, cycle_scores[cyc1]), (cyc2, cycle_scores[cyc2])], key=lambda tup: tup[1])
        collapsed_cycle[sorted_cyc_scores[0][0]] = True

ofn_sum = open(args.ofn_all_cycle_sum, "w+")
ofn_sum.write(util.write_line("sample", "cycle", "bottleneck", "median_support", "non_support_cov", "avg_external", "avg_paired_weird", "avg_singleton", "pval", "class", "score", "lower_bound_cov", "length"))
ofn_sum_dom = open(args.ofn_dominant_cycle_sum, "w+")
ofn_sum_dom.write(util.write_line("sample", "cycle", "bottleneck", "median_support", "non_support_cov", "avg_external", "avg_paired_weird", "avg_singleton", "pval", "class", "score", "lower_bound_cov", "length"))
good_cycles = []
good_cycles_to_stats = {}
with open(args.ifn_cycle_sum) as sum:
    header = util.get_header_dict(args.ifn_cycle_sum)
    next(sum)
    for line in sum:
        line_s = util.split(line)
        # cycle_length = line[header[""]]
        sample = line_s[header["sample"]]
        cycle = line_s[header["cycle"]]
        avg_out = (float(line_s[header["non_support_out_reads"]]) + float(line_s[header["non_support_in_reads"]])) / 2
        avg_paired_weird = args.paired_weird_mult * (float(line_s[header["paired_weird_out_reads"]]) + float(line_s[header["paired_weird_in_reads"]])) / 2
        total_non_support = avg_out + avg_paired_weird
        total_num_reads = float(line_s[header["paired_read_count"]])
        bottleneck_reads = float(line_s[header["bottleneck_reads"]])
        pval = binom.cdf(total_non_support, total_num_reads, bottleneck_reads / total_num_reads) # score = 2, yes filter paired weird added
        avg_external = (float(line_s[header["non_support_out"]]) + float(line_s[header["non_support_in"]])) / 2
        avg_paired_weird_cov = (float(line_s[header["paired_weird_out"]]) + float(line_s[header["paired_weird_in"]])) / 2
        avg_singleton = (float(line_s[header["singletons_out"]]) + float(line_s[header["singletons_in"]])) / 2
        bottleneck = float(line_s[header["bottleneck"]])
        median_support = float(line_s[header["avg_support_cov"]])
        non_support_cov = avg_external + args.paired_weird_mult * avg_paired_weird_cov
        low_coverage_est = bottleneck - non_support_cov
        if non_support_cov == 0 and bottleneck == 0:
            score = 0
        elif non_support_cov == 0 and bottleneck != 0:
            score = "inf"
        else:
            score = round(bottleneck / non_support_cov, 3)
        if pval < args.max_pval:
            cycle_class = "dominant"
            dominant = True
        else:
            dominant = False
            cycle_class = "not_dominant"
        ofn_sum.write(util.write_line(sample, cycle, bottleneck, median_support, round(non_support_cov,3), round(avg_external, 3), round(avg_paired_weird_cov,3), round(avg_singleton, 3), round(pval, 5), cycle_class, score, round(low_coverage_est,3), cyc_to_length[cycle]))
        if dominant:
            ofn_sum_dom.write(util.write_line(sample, cycle, bottleneck, median_support, round(non_support_cov,3), round(avg_external,3), round(avg_paired_weird_cov,3), round(avg_singleton, 3), round(pval, 5), cycle_class, score, round(low_coverage_est,3), cyc_to_length[cycle]))
            good_cycles.append(cycle)

ofn_contigs = open(args.ofn_cycle_contigs,"w+")
with open(args.ifn_cycle_contigs) as contigs:
    header = util.get_header_dict(args.ifn_cycle_contigs)
    head_line = contigs.readline()
    ofn_contigs.write(head_line)
    for line in contigs:
        line_s = util.split(line)
        cycle = line_s[header["cycle"]]
        if cycle in good_cycles:
            ofn_contigs.write(line)

ofn = open(args.ofn_cycle_fasta,"w+")
seqs = util.parse_fasta(args.ifn_cycle_fasta)
for cyc in good_cycles:
    ofn.write(">" + cyc + "\n" + seqs[cyc] + "\n")