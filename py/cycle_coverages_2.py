import argparse
import util
import numpy as np

parser = argparse.ArgumentParser(description='Determine basepair cycle coverage values...')
parser.add_argument('sample', metavar='<sample name>', type=str, help='sample name')
parser.add_argument('ifn_cyc_table', metavar='<cycle table>', type=str, help='cycle table')
parser.add_argument('ifn_paired_map', metavar='<bwa paired read table>', type=str, help='bwa paired read table')
parser.add_argument('ifn_singleton_map', metavar='<bwa singleton read table>', type=str, help='bwa singleton read table')
parser.add_argument('ifn_contig_table', metavar='<cycle contig table>', type=str, help='cycle contig table')
parser.add_argument('read_stats', metavar='<read stats>', type=str, help='read stats file')
parser.add_argument('ofn_cyc_table', metavar='<output cycle table>', type=str, help='base-breakdown of coverge stats for each cycle')
parser.add_argument('ofn_cyc_coverage_summary', metavar='<cycle coverage summary>', type=str, help='cycle coverage summary')
parser.add_argument('ofn_out_contigs', metavar='<contigs leading out of cycles>', type=str, help='contigs leading out of cycles')

args = parser.parse_args()

k = util.get_k(args.read_stats, True)
max_distance = util.get_max_distance(args.read_stats)
read_len = int(util.get_read_length(args.read_stats, True))
insert = int(util.get_read_insert(args.read_stats))
print("Insert", insert)

print("Max insert distance", max_distance)
print("Fetching cycle and contig information...")

contig_to_cyc, cyc_lengths, cycle_to_contig = util.get_cycle_specs(args.ifn_cyc_table)
cycle_to_covs = {}
cycle_nonsupport_in = {}
cycle_nonsupport_out = {}
cycle_pairedweird_in = {}
cycle_pairedweird_out = {}
cycle_singletons_in = {}
cycle_singletons_out = {}
for cycle, length in cyc_lengths:
    cycle_to_covs[cycle] = np.zeros((7,length), dtype=int) # 1: support 2: out non-support 3: in non-support 4: out paired-weird 5: in paired-weird 6: out single 7: in single
    cycle_nonsupport_in[cycle] = 0
    cycle_nonsupport_out[cycle] = 0
    cycle_pairedweird_in[cycle] = 0
    cycle_pairedweird_out[cycle] = 0
    cycle_singletons_in[cycle] = 0
    cycle_singletons_out[cycle] = 0

contig_to_cov, contig_to_length = util.get_contig_stats(args.ifn_contig_table)
print("Iterating through bwa table...")

trans_contigs = {} # k: cycle v: [{(out_contig, #reads)}]

line_count = 0
with open(args.ifn_paired_map) as map:
    header = util.get_header_dict(args.ifn_paired_map)
    next(map)
    for line in map:
        line_count += 1
        if line_count % 1e6 == 0:
            print("mapping line count:", line_count)
        line = util.split(line)
        contig1, contig2 = line[header['contig1']], line[header['contig2']]
        if contig_to_cyc.get(contig1) is None and contig_to_cyc.get(contig2) is None: # both reads don't map to a cycle
            continue
        contig1_length = contig_to_length[contig1]
        contig2_length = contig_to_length[contig2]
        read1_start, read2_start = int(line[header['coord1']]), int(line[header['coord2']])
        read1_stop, read2_stop = int(line[header['back_coord1']]), int(line[header['back_coord2']])
        dist1, dist2 = int(line[header['edit_dist1']]), int(line[header['edit_dist2']])
        read1_usable = util.read_not_in_overlap(contig1_length, read1_start, read1_stop, k, dist1) # read_not_in_overlap() returns True if read doesn't belong to contig_cyc
        read2_usable = util.read_not_in_overlap(contig2_length, read2_start, read2_stop, k, dist2)
        read_strand1, read_strand2 = int(line[header['strand1']]), int(line[header['strand2']])
        if not (read1_usable and read2_usable) or not (line[header["unique1"]] == "T" and line[header["unique2"]] == "T"): # if either read in k-overlap regions for reads that map to cycle contigs, then skip
            continue
        cycles = util.contig_to_cycs(contig1, contig2, contig_to_cyc, read_strand1, read_strand2)
        for (cycle1,cycle2), (cumsum1, cumsum2), (contig1_len, contig2_len), (ori1, ori2), (cyc1_len,cyc2_len) in cycles:
            if cycle1 == None and cycle2 == None: continue
            if cycle1 != None: # contig1 belongs to a cycle
                cycle_coord1_start, cycle_strand1 = util.get_cycle_coords(cumsum1, contig1_len, ori1, read1_start, read1_stop)
            else:
                cycle_coord1_start, cycle_strand1 = None, None # contig1 doesn't belong to a cycle
            if cycle2 != None: # contig2 belongs to a cycle
                cycle_coord2_start, cycle_strand2 = util.get_cycle_coords(cumsum2, contig2_len, ori2, read2_start, read2_stop)
            else:
                cycle_coord2_start, cycle_strand2 = None, None
            if cycle1 == cycle2: # True = cycle supporting or paired weird
                sorted_coords = sorted([(cycle_coord1_start, cycle_strand1),(cycle_coord2_start, cycle_strand2)], key=lambda tup: tup[1], reverse=True) # arrange by which one is cycle strand 1 and -1
                start, stop = sorted_coords[0][0], sorted_coords[1][0]
                if start > stop: # read goes around the cycle?
                    total_distance = cyc1_len - start + stop
                    if total_distance > max_distance or cycle_strand1 == cycle_strand2: # weird paired reads either on same strand or too long away from each other
                        for coord,strand in sorted_coords:
                            if strand == 1:
                                # stop = min(cyc1_len - 1, coord + read_len)
                                stop = min(cyc1_len - 1, coord + insert)
                                cycle_pairedweird_out[cycle1] += 1
                                cycle_to_covs[cycle1][4][(coord - 1):stop] += 1
                            else:
                                start = max(0, coord)
                                cycle_pairedweird_in[cycle1] += 1
                                cycle_to_covs[cycle1][3][start:coord] += 1
                    else: # supporting read!
                        cycle_to_covs[cycle1][0][(start - 1):] += 1
                        cycle_to_covs[cycle1][0][:stop] += 1
                else:
                    total_distance = stop - start
                    if total_distance > max_distance or cycle_strand1 == cycle_strand2: # weird paired non-supporting
                        for coord,strand in sorted_coords:
                            if strand == 1:
                                stop = min(cyc1_len - 1, (coord + insert))
                                cycle_pairedweird_out[cycle1] += 1
                                cycle_to_covs[cycle1][4][(coord - 1):stop] += 1
                            else:
                                start = max(0, coord)
                                cycle_pairedweird_in[cycle1] += 1
                                cycle_to_covs[cycle1][3][(start):coord] += 1
                    else: # supporting read!
                        cycle_to_covs[cycle1][0][(start-1):stop] += 1
            else: # cycle non-supporting
                cycle_stats = [(cycle1, cycle_coord1_start, cycle_strand1) , (cycle2, cycle_coord2_start, cycle_strand2)]
                for cycle_ind, (cycle, coord_start, cyc_strand) in enumerate(cycle_stats):
                    cycle_contigs = cycle_to_contig.get(cycle,[None]) # [None] for if the cycle == None which is handled in next step
                    if (contig1 in cycle_contigs and contig2 in cycle_contigs) or (cycle == None): # read doesn't map to a cycle or both contigs belong to cycle, so don't count as non-support
                        continue
                    if cycle_ind == 0:
                        out_contig = contig2
                    else:
                        out_contig = contig1
                    trans_contigs.setdefault(cycle, {})  # get the out contigs and reads connected to other contigs
                    out_contig_reads = trans_contigs[cycle].get(out_contig, 0)
                    trans_contigs[cycle][out_contig] = out_contig_reads + 1
                    if cyc_strand == 1: # out non-supporting cov
                        cycle_nonsupport_out[cycle] += 1
                        cyc_length = len(cycle_to_covs[cycle][0])
                        # cycle_to_covs[cycle][1][(coord_start - 1):(coord_start + read_len)] += 1
                        cycle_to_covs[cycle][1][(coord_start - 1):(coord_start + insert)] += 1
                    else: # in non-supporting cov
                        cycle_nonsupport_in[cycle] += 1
                        # stop = coord_start - read_len
                        stop = coord_start - insert
                        if stop < 0:
                            stop = 0
                        cycle_to_covs[cycle][2][stop:coord_start] += 1

with open(args.ifn_singleton_map) as singletons:
    header = util.get_header_dict(args.ifn_singleton_map)
    next(singletons)
    for line in singletons:
        line = util.split(line)
        contig = line[header["contig"]]
        if contig_to_cyc.get(contig) is not None:
            contig_length = contig_to_length[contig]
            read_start, read_stop = int(line[header['coord']]), int(line[header['back_coord']])
            read_usable = util.read_not_in_overlap(contig_length, read_start, read_stop, k)  # read_not_in_overlap() returns True if read doesn't belong to contig_cyc
            read_strand = int(line[header['strand']])
            if not read_usable:  # if either read in k-overlap regions for reads that map to cycle contigs, then skip
                continue
            cycles = util.contig_to_cycs(contig, "singleton", contig_to_cyc)
            for cycle, cumsum, contig_len, ori, cyc_len in cycles:
                cycle_coord_start, cycle_strand = util.get_cycle_coords(cumsum, contig_len, ori, read_start, read_stop)
                if cycle_strand == 1:
                    # stop = min(cyc_len, cycle_coord_start + read_len)
                    stop = min(cyc_len, cycle_coord_start + insert)
                    cycle_singletons_out[cycle] += 1
                    cycle_to_covs[cycle][6][cycle_coord_start-1:stop] += 1
                else:
                    # start = max(0, cycle_coord_start - read_len)
                    start = max(0, cycle_coord_start - insert)
                    cycle_singletons_in[cycle] += 1
                    cycle_to_covs[cycle][5][start:cycle_coord_start] += 1

print("Writing results...")

ostats = open(args.ofn_cyc_coverage_summary,"w+")
ostats.write(util.write_line("sample","cycle", "avg_support_cov","mean_support", "bottleneck", "bottleneck_reads", "bottleneck_pos", "max_nonsupport", "max_nonsupport_reads","non_support_out", "non_support_out_reads", "non_support_in",
                             "non_support_in_reads", "paired_weird_out", "paired_weird_out_reads", "paired_weird_in","paired_weird_in_reads", "singletons_out", "singletons_in","cyc_length", "paired_read_count"))
ofile = open(args.ofn_cyc_table,"w+")
ofile.write(util.write_line("cycle", "base", "support", "out_cov", "in_cov", "out_paired_weird", "in_paired_weird", "out_singleton", "in_singleton"))

cov_factor = (2 * read_len) / float(insert) # factor to get coverage bins into x-coverage units
cycle_count = 0
for cycle, length in cyc_lengths:
    if cycle_count % 100 == 0:
        print("Print cycle",cycle_count)
    cycle_count += 1
    cycle_covs = cycle_to_covs[cycle]
    bottleneck_reads = np.amin(cycle_covs[0])
    bottleneck = round(cov_factor * bottleneck_reads, 2)
    bottleneck_pos = np.argmin(cycle_covs[0]) + 1
    avg_support = round(cov_factor * int(np.median(cycle_covs[0])), 2)
    mean_support = round(cov_factor * int(np.mean(cycle_covs[0])), 2)
    out_nonsupport = cycle_nonsupport_out[cycle]
    in_nonsupport = cycle_nonsupport_in[cycle]
    max_nonsupport = max(out_nonsupport, in_nonsupport)
    pairedweird_out = cycle_pairedweird_out[cycle]
    pairedweird_in = cycle_pairedweird_in[cycle]
    singletons_out = cycle_singletons_out[cycle]
    singletons_in = cycle_singletons_in[cycle]
    ostats.write(util.write_line(args.sample, cycle, avg_support, mean_support, bottleneck, bottleneck_reads, bottleneck_pos, round(max_nonsupport * cov_factor, 2), max_nonsupport, round(out_nonsupport * cov_factor, 2), out_nonsupport, round(in_nonsupport * cov_factor, 2), in_nonsupport,
                                 round(pairedweird_out * cov_factor, 2), pairedweird_out, round(pairedweird_in * cov_factor, 2), pairedweird_in, singletons_out, singletons_in, length, line_count))
    for base in range(length):
        ofile.write(util.write_line(cycle, base + 1, cov_factor * cycle_covs[0][base], cov_factor * cycle_covs[1][base], cov_factor * cycle_covs[2][base], cov_factor * cycle_covs[3][base], cov_factor * cycle_covs[4][base], cov_factor * cycle_covs[5][base], cov_factor * cycle_covs[6][base]))

print("Writing out contigs...")
ofn_out_contigs = open(args.ofn_out_contigs, "w+")
ofn_out_contigs.write(util.write_line("cycle","out_contig","num_reads"))
for cycle in trans_contigs.keys():
    for out_contig in trans_contigs[cycle].keys():
        ofn_out_contigs.write(util.write_line(cycle, out_contig, trans_contigs[cycle][out_contig]))




