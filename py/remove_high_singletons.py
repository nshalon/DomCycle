import argparse
import util
import sys
from scipy.stats import binom
import os

###############################################################################################
# Purpose: filter on local score                                                              #
###############################################################################################

parser = argparse.ArgumentParser(description='remove cycles with high singletons')
parser.add_argument('ifn_long_table', metavar='<long table of cycle coverage>', type=str, help='base pair level resolution of coverage')
parser.add_argument('ifn_cyc_stats_all', metavar='<cycle stats>', type=str, help='cycle stats')
parser.add_argument('ifn_cyc_contig_table', metavar='<cycle contig table>', type=str, help='cycle contig table')
parser.add_argument('ifn_cycle_cov_summary', metavar='<cycle coverage summary file>', type=str, help='cycle coverage summary file')
parser.add_argument('ifn_cyc_fa', metavar='<cycle fasta>', type=str, help='cycle fasta')
parser.add_argument('read_stats', metavar='<library read stats>', type=str, help='library read stats')
parser.add_argument('min_singleton_score', metavar='<minimum local score to threshold cycles>', type=float, help='minimum local score to threshold cycles')
parser.add_argument('max_singleton_score_pval', metavar='<max local score pval>', type=float, help='max local score pval')
parser.add_argument('odir', metavar='<out directory>', type=str, help='out directory')

args = parser.parse_args()

dominant_cycles = set()
with open(args.ifn_cycle_cov_summary) as cycle_cov_summary:
    header = util.get_header_dict(args.ifn_cycle_cov_summary)
    cycle_cov_summary.readline()
    line = util.split(cycle_cov_summary.readline())
    print(line)
    if len(line) > 1:
        total_num_reads = int(line[header["paired_read_count"]])
    else:
        total_num_reads = 0
        odominant_stats = open(os.path.join(args.odir, "dominant_cycle_stats"), "w+")
        ostats = open(os.path.join(args.odir, "cycle_stats_singletons"), "w+")
        odominant_stats.write(
            util.write_line("sample", "cycle", "bottleneck", "median_support", "avg_external", "avg_inter",
                            "avg_intra-nonsupport", "avg_singleton", "pval", "class", "score", "lower_bound_cov",
                            "length", "MDC", "singleton_score", "singleton_score_pval"))
        ostats.write(util.write_line("sample", "cycle", "bottleneck", "median_support", "avg_external", "avg_inter",
                                     "avg_intra-nonsupport", "avg_singleton", "pval", "class", "score",
                                     "lower_bound_cov", "length", "MDC", "singleton_score", "singleton_score_pval"))
        ofasta = open(os.path.join(args.odir, "cycles.fasta"), "w+")

factor = (2 * util.get_read_length(args.read_stats, True)) / util.get_read_insert(args.read_stats)

dom_cycles_not_passing = set()
cycle_to_singleton_stats = {}
with open(args.ifn_long_table) as long_table:
    header = util.get_header_dict(args.ifn_long_table)
    next(long_table)
    cycle = "CYC1"
    singleton_score = sys.maxsize
    min_score = sys.maxsize
    min_score_stats = (sys.maxsize, (sys.maxsize, 0))
    for line in long_table:
        line = util.split(line)
        line_cycle = line[header["cycle"]]
        if line_cycle != cycle:
            singleton_pval = binom.cdf(min_score_stats[1][1], total_num_reads, min_score_stats[1][0] / total_num_reads)
            cycle_to_singleton_stats[cycle] = (min_score_stats[0], singleton_pval)
            min_score = sys.maxsize
            cycle = line_cycle
        singleton_and_nonsupport = (1 / factor) * max(float(line[header["in_singleton"]]) + float(line[header["in_paired_weird"]]) + float(line[header["in_cov"]]), float(line[header["out_singleton"]]) + float(line[header["out_paired_weird"]]) + float(line[header["out_cov"]]))
        support = (1 / factor) * float(line[header["support"]])
        if singleton_and_nonsupport != 0:
            singleton_score = support / singleton_and_nonsupport
        else:
            singleton_score = sys.maxsize
        if singleton_score < min_score:
            min_score = singleton_score
            min_score_stats = (singleton_score, (support, singleton_and_nonsupport))
    singleton_pval = binom.cdf(min_score_stats[1][1], total_num_reads, min_score_stats[1][0] / total_num_reads)
    cycle_to_singleton_stats[cycle] = (min_score_stats[0], singleton_pval)

odominant_stats = open(os.path.join(args.odir, "dominant_cycle_stats_singletons"), "w+")
ostats = open(os.path.join(args.odir, "cycle_stats_singletons"), "w+")
odominant_stats.write(util.write_line("sample", "cycle", "bottleneck", "median_support", "avg_external", "avg_inter", "avg_intra-nonsupport", "avg_singleton", "pval", "class", "score", "lower_bound_cov", "length", "MDC", "singleton_score","singleton_score_pval"))
ostats.write(util.write_line("sample", "cycle", "bottleneck", "median_support", "avg_external", "avg_inter", "avg_intra-nonsupport", "avg_singleton", "pval", "class", "score", "lower_bound_cov", "length", "MDC", "singleton_score","singleton_score_pval"))
singleton_filtered_cycles = set()
with open(args.ifn_cyc_stats_all) as cycle_stats:
    header = util.get_header_dict(args.ifn_cyc_stats_all)
    next(cycle_stats)
    for line in cycle_stats:
        line_s = util.split(line)
        singleton_score, singleton_pval = cycle_to_singleton_stats[line_s[header["cycle"]]]
        if singleton_score == sys.maxsize:
            singleton_score_parse = "inf"
        else:
            singleton_score_parse = str(round(singleton_score, 2))
        external_pval = float(line_s[header["pval"]])
        mdc = str(round(float(line_s[header["median_support"]]) - float(line_s[header["avg_external"]]), 3))
        ostats.write(line[:-1] + "\t" + mdc + "\t" + singleton_score_parse + "\t" + str(round(singleton_pval, 5)) + "\n")
        if singleton_score > args.min_singleton_score and singleton_pval < args.max_singleton_score_pval and line_s[header["class"]] == "dominant":
            odominant_stats.write(line[:-1] + "\t" + mdc + "\t" + singleton_score_parse + "\t" + str(round(singleton_pval, 5)) + "\n")
            singleton_filtered_cycles.add(line_s[header["cycle"]])

ocontigs = open(os.path.join(args.odir, "cycle_contig_table"), "w+")
with open(args.ifn_cyc_contig_table) as icontigs:
    header = util.get_header_dict(args.ifn_cyc_contig_table)
    ocontigs.write(icontigs.readline())
    for line in icontigs:
        line_split = util.split(line)
        if line_split[header["cycle"]] in singleton_filtered_cycles:
            ocontigs.write(line)
ofasta = open(os.path.join(args.odir, "cycles.fasta"), "w+")
cycle_seqs = util.parse_fasta(args.ifn_cyc_fa)
for cycle, seq in cycle_seqs.items():
    if cycle in singleton_filtered_cycles:
        ofasta.write(">" + cycle + "\n" + seq + "\n")


