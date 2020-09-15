import argparse
import util
import numpy as np
import sys
#############################################################################
# Purpose: find shortest cycle back to (-) node of each (+) node. Start point is (+) node for each (a) strand
# Input: all nodes parsed and renamed from fastg, edges, and weights
# Output: table of all shortest cycles, a separate adjacency matrix for each, and a weight list for r graphs
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_edge_summary', metavar='<summary stats for each edge>', type=str, help='input fastq files')
parser.add_argument('ifn_read_stats', metavar='<kmer size>', type=str, help='kmer size used for assembly')
parser.add_argument('alpha', metavar='<out coverage dampening factor>', type=float, help='out coverage dampening factor')
parser.add_argument('ofn_table', metavar='<adjacency matrix for edge>', type=str, help='input fastq files')

args = parser.parse_args()

print("Collecting all edge statistics now...")
parent_to_children, children_to_parents, edge_to_cov, start_edges, edge_type, edge_to_length = util.get_edge_structures(args.ifn_edge_summary)

alpha = args.alpha # out coverage dampening factor (a=1 means score > 1 to report further)

print("Finding dominant cycles over", len(start_edges), "edges")
# print(start_edges)

ofile = open(args.ofn_table, "w+")
ofile.write(util.write_line("cycle","original_contig","FASTG_edge","renamed_contig","contig_len","cum_sum","orientation","cycle_length","contig_index"))
k = util.get_k(args.ifn_read_stats, True)

out = lambda in_weight, out_weight: alpha * ((out_weight + in_weight) / 2)
orientations = {"+": "-", "-": "+"}
previously_found_cycles = set()
cycle_num = 0
for start_edge in start_edges:
    bottleneck = edge_to_cov[start_edge]
    out_weight = 0
    in_weight = 0
    start_head, start_tail = start_edge
    bottleneck = min(edge_to_cov[start_edge], max([edge_to_cov[(start_tail,v)] for v in parent_to_children.get(start_tail,[])], default=sys.maxsize)) # bottleneck is the min(w(start_edge),max(w(neighbors)))
    cycle_path = [start_edge]
    cycle_contigs = {start_tail, start_head}
    is_open = True
    fringe = start_tail
    while out(in_weight, out_weight) < bottleneck and is_open:
        all_successors = parent_to_children.get(fringe, [])[:]
        successor_covs = [edge_to_cov[(fringe, s)] for s in all_successors]
        if len([cov for cov in successor_covs if cov >= bottleneck]) < 1:
            break
        successor = all_successors.pop(np.argmax(successor_covs))
        is_open = successor not in cycle_contigs
        in_weight += sum([edge_to_cov[(in_neighbor, successor)] for in_neighbor in children_to_parents[successor] if in_neighbor != fringe])
        out_weight += sum([edge_to_cov[(fringe, out_neighbor)] for out_neighbor in parent_to_children[fringe] if out_neighbor != successor])
        cycle_path.append((fringe, successor))
        cycle_contigs.add(successor)
        fringe = successor
    if fringe == start_head and out(in_weight, out_weight) < bottleneck: # dominant cycle exists
        sorted_path = tuple(sorted(cycle_contigs))
        complement_sorted_path = util.get_sorted_comp_path(cycle_contigs)
        if sorted_path in previously_found_cycles or complement_sorted_path in previously_found_cycles:
            continue
        previously_found_cycles.add(sorted_path)
        previously_found_cycles.add(complement_sorted_path)
        cycle_num += 1
        cycle = "CYC" + str(cycle_num)
        cycle_length = int(sum([edge_to_length[e] for i, e in enumerate(cycle_path) if i % 2 == 0]) - (k * (len(cycle_path) / 2)))
        print("Found", cycle, "with out weight", round(out(in_weight, out_weight), 2), "and length", cycle_length)
        cum_sum = 1 # starting cycle position of contig
        ori = "+"
        contig_strand = start_tail[-2]
        contig_index = 0
        # print("cycle",cycle_path)
        for edge_num, edge in enumerate(cycle_path): # just for writing out appropriate output
            if edge_num % 2 == 0:
                contig_len = edge_to_length[edge]
            else:
                contig_index += 1
                contig = edge[0][:-2]
                is_fastg = edge_type[edge]
                ofile.write(util.write_line(cycle, contig, is_fastg, contig, contig_len, cum_sum, ori, cycle_length, contig_index))
                cum_sum += contig_len - k
                if contig_strand != edge[1][-2]: # switch orientation of contig in cycle path if necessary
                    contig_strand = edge[1][-2]
                    ori = orientations[ori]

