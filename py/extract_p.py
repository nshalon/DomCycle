import argparse
import util
import itertools
from scipy.stats import binom

###############################################################################################
# Purpose: filter on global score                                                             #
###############################################################################################

parser = argparse.ArgumentParser(description='global score filtering')
parser.add_argument('ifn_cycle_sum', metavar='<input cycle summary>', type=str, help='input cycle summary')
parser.add_argument('ifn_cycle_fasta', metavar='<fasta of input cycles>', type=str, help='fasta of input cycles')
parser.add_argument('ifn_cycle_contigs', metavar='<input cycle contig table>', type=str, help='input cycle contig table')
parser.add_argument('max_pval', metavar='<max pval for global score>', type=float, help='max pval for global score')
parser.add_argument('min_length', metavar='<min cycle length to report>', type=int, help='minimum cycle length to report')
parser.add_argument('min_lbc', metavar='<min AMC to report>', type=float, help='min AMC to report')
parser.add_argument('paired_weird_mult', metavar='<multiple on intra-nonsupport coverage>', type=float, help='multiple on intra-nonsupport coverage')
parser.add_argument('min_score', metavar='<minimum score to threshold on>', type=float, help='minimum score to threshold on')
parser.add_argument('ofn_all_cycle_sum', metavar='<output cycle summary>', type=str, help='output cycle summary')
parser.add_argument('ofn_dominant_cycle_sum', metavar='<output global score cycle table>', type=str, help='output global score cycle table')
parser.add_argument('ofn_cycle_fasta', metavar='<output global score cycle fasta>', type=str, help='output global score cycle fasta')
parser.add_argument('ofn_cycle_contigs', metavar='<output global score dominant cycle contig table>', type=str, help='output global score dominant cycle contig table')

args = parser.parse_args()

cyc_to_length, cyc_to_contig, contig_to_length = util.get_cyc_contig_stats(args.ifn_cycle_contigs)
cycle_scores = util.get_cycle_scores(args.ifn_cycle_sum)

ofn_sum = open(args.ofn_all_cycle_sum, "w+")
ofn_sum.write(util.write_line("sample", "cycle", "bottleneck", "median_support", "avg_external", "avg_inter", "avg_intra-nonsupport", "avg_singleton", "pval", "class", "score", "lower_bound_cov", "length"))
ofn_sum_dom = open(args.ofn_dominant_cycle_sum, "w+")
ofn_sum_dom.write(util.write_line("sample", "cycle", "bottleneck", "median_support", "avg_external", "avg_inter", "avg_intra-nonsupport", "avg_singleton", "pval", "class", "score", "lower_bound_cov", "length"))
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
        avg_paired_weird = (float(line_s[header["paired_weird_out_reads"]]) + float(line_s[header["paired_weird_in_reads"]])) / 2
        total_non_support = avg_out + args.paired_weird_mult * avg_paired_weird
        total_num_reads = float(line_s[header["paired_read_count"]])
        bottleneck_reads = float(line_s[header["bottleneck_reads"]])
        pval = binom.cdf(total_non_support, total_num_reads, bottleneck_reads / (args.min_score * total_num_reads)) # score = 1
        if bottleneck_reads / (args.min_score * total_num_reads) > 1: # when score is changed
            pval = 0
        non_support_cov = (float(line_s[header["non_support_out"]]) + float(line_s[header["non_support_in"]])) / 2
        avg_paired_weird_cov = (float(line_s[header["paired_weird_out"]]) + float(line_s[header["paired_weird_in"]])) / 2
        avg_singleton = (float(line_s[header["singletons_out"]]) + float(line_s[header["singletons_in"]])) / 2
        bottleneck = float(line_s[header["bottleneck"]])
        median_support = float(line_s[header["avg_support_cov"]])
        avg_external = non_support_cov + args.paired_weird_mult * avg_paired_weird_cov
        low_coverage_est = bottleneck - non_support_cov
        if non_support_cov == 0 and bottleneck == 0:
            score = 0
        elif avg_external == 0:
            score = "inf"
        else:
            score = round(bottleneck / avg_external, 3)
        if pval < args.max_pval:
            cycle_class = "dominant"
            dominant = True
        else:
            dominant = False
            cycle_class = "not_dominant"
        ofn_sum.write(util.write_line(sample, cycle, bottleneck, median_support, round(avg_external, 3), round(non_support_cov,3), round(avg_paired_weird_cov,3), round(avg_singleton, 3), round(pval, 5), cycle_class, score, round(low_coverage_est,3), cyc_to_length[cycle]))
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