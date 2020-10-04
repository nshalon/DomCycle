import argparse
import util

###############################################################################################
# Purpose: input adjacency list, paired mapped reads, kmer size, contig table and             #
# outputs variable edge coverage                                                               #
###############################################################################################

parser = argparse.ArgumentParser(description='get variable edge coverages')
parser.add_argument('ifn_paired_reads', metavar='<table of bwa mapped paired reads>', type=str, help='bwa mapped paired reads')
parser.add_argument('ifn_adj_list', metavar='<adjacency list from fastg>', type=str, help='adjacency list for fastg')
parser.add_argument('ifn_contigs', metavar='<contig table input>', type=str, help='all contig coverage and length')
parser.add_argument('ifn_read_stats', metavar='<read length .txt>', type=str, help='size of the reads')
parser.add_argument('ofn_covs', metavar='<output coverages for variable edges>', type=str, help='output coverages of variable edges')
args = parser.parse_args()

read_length = util.get_read_length(args.ifn_read_stats, True)
insert = util.get_read_insert(args.ifn_read_stats)

k = util.get_k(args.ifn_read_stats, True)
max_distance = util.get_max_distance(args.ifn_read_stats)

factor = ( 2 * read_length) / ( insert )
correction = 1 / (1 - ( k / ( insert + k ) ) )

bwa_head = util.get_header_dict(args.ifn_paired_reads) # get the dict that returns field to index

parent_to_child, connected  = util.get_parent_to_child(args.ifn_adj_list) # parse the adj list

contig_to_cov, contig_to_length = util.get_contig_stats(args.ifn_contigs) # get contig length and coverages

edge_to_coverage = {}
nongraph_edge_to_coverage = {}
nongraph_edge_parent_to_child = {}
with open(args.ifn_paired_reads) as paired_reads:
    next(paired_reads)
    for line in paired_reads:
        line = util.split(line)
        contig1, contig2 = line[bwa_head['contig1']],line[bwa_head['contig2']]
        contig1_len, contig2_len = contig_to_length[contig1],contig_to_length[contig2]
        read1_start, read2_start = int(line[bwa_head['coord1']]),int(line[bwa_head['coord2']])
        read1_stop, read2_stop = int(line[bwa_head['back_coord1']]), int(line[bwa_head['back_coord2']])
        dist1, dist2 = int(line[bwa_head['edit_dist1']]), int(line[bwa_head['edit_dist2']])
        strand1, strand2 = int(line[bwa_head['strand1']]), int(line[bwa_head['strand2']])
        if contig1 == contig2:
            sorted_starts_coords = sorted([(read1_start,strand1), (read2_start,strand2)], key=lambda start_strand: start_strand[0])
            if sorted_starts_coords[0][1] == 1 or (sorted_starts_coords[1][0] - sorted_starts_coords[0][0] < contig1_len/2): # if the positive facing read falls before the negative facing read (as it should for normal intra-contig reads), skip
                continue
        read1_usable = util.read_not_in_overlap(contig1_len, read1_start, read1_stop, k, dist1)
        read2_usable = util.read_not_in_overlap(contig2_len, read2_start, read2_stop, k, dist2)
        read_far = util.get_reads_far(read1_stop, contig1_len, strand1 , read2_stop, contig2_len, strand2, max_distance)
        if (read1_usable and read2_usable) and not read_far:
            suffix1, suffix2 = util.get_node_suffix(read1_start,read1_stop,read2_start,read2_stop)
            node1 = contig1 + suffix1
            node2 = contig2 + suffix2
            edge = (node1, node2)
            coverage = edge_to_coverage.get(edge, 0)
            edge_to_coverage[edge] = coverage + 1
            comp_edge = util.get_comp_edge(edge)
            coverage = edge_to_coverage.get(comp_edge, 0)
            edge_to_coverage[comp_edge] = coverage + 1
            if connected.get(edge) == None: # if nodes not connected in graph
                if node2 not in nongraph_edge_parent_to_child.get(node1, []):
                    nongraph_edge_parent_to_child.setdefault(node1, []).append(node2)
                    nongraph_edge_parent_to_child.setdefault(comp_edge[0], []).append(comp_edge[1])

ofile = open(args.ofn_covs,"w+")
ofile.write(util.write_line("parent","child","coverage","length","weight","graph_edge"))
with open(args.ifn_adj_list) as adj_list:
    for line in adj_list:
        parent,g_children = util.get_parent_children(line)
        children = [(child,True) for child in g_children]
        if parent[-1] == "-":
            invariable_edge = True
        else: invariable_edge = False
        if nongraph_edge_parent_to_child.get(parent) != None: # var edge with children not in graph
            non_graph_children = [(child,False) for child in nongraph_edge_parent_to_child[parent] ]
            children += non_graph_children
        children.sort(key=lambda child: child[0])
        for child in children:
            if child == '':
                continue
            if invariable_edge:
                coverage = contig_to_cov[parent[:-2]]
                if coverage == 0: coverage = 1e-10  # avoid division by 0
                length = contig_to_length[parent[:-2]]
                weight = round(float(length)/float(coverage),5)
                in_graph = True
            else:
                in_graph = child[1]
                coverage = round(edge_to_coverage.get((parent,child[0]),0) * factor * correction, 3)
                length = 200
                if coverage == 0: coverage = 1e-10 # avoid division by 0
                weight = round(float(length)/coverage,5)
            child = child[0]
            ofile.write(util.write_line(parent,child,coverage,length,weight,in_graph))













