import numpy as np

####################################
# Generic functions                #
####################################

def get_read_length(file, paired=False):
    if paired:
        with open(file) as file:
            next(file)
            line = split(file.readline().rstrip())
            return float(line[0])
    else:
        with open(file) as file:
            return float(file.readline().rstrip())

def get_read_insert(file):
    with open(file) as file:
        next(file)
        line = split(file.readline().rstrip())
        return float(line[2])

def get_max_distance(file):
    with open(file) as file:
        next(file)
        line = split(file.readline().rstrip())
        return float(line[3])

def get_k(file, paired=False):
    if paired:
        with open(file) as file:
            next(file)
            line = split(file.readline().rstrip())
            return float(line[1])
    else:
        with open(file) as file:
            next(file)
            return float(file.readline().rstrip())

def split(line):
    return line.rstrip().split("\t")

def write_line(*argv):
    line = ''
    for arg in argv:
        line += str(arg) + "\t"
    return line[:-1] + "\n"

def get_header_dict(file):
    with open(file) as file:
        header = split(file.readline())
        header_dict = {}
        for index,field in enumerate(header):
            header_dict[field] = index
    return header_dict

def parse_adj_list(file):
    parent_to_child = {}
    child_to_parent = {}
    start_nodes = []
    with open(file) as nodes:
        for edge_list in nodes:
            edges = edge_list.rstrip().split(':')
            parent = edges[0]
            children = edges[1].split(',')
            for child in children:
                if child == '':
                    continue
                parent_to_child.setdefault(parent, []).append(child)
                child_to_parent.setdefault(child, []).append(parent)
            if parent[-2:] == 'a+':
                start_nodes.append(parent)
    return parent_to_child,child_to_parent,start_nodes

def get_parent_to_child(adj_list_file):
    parent_to_child = {}
    connected = {}
    with open(adj_list_file) as nodes:
        for edge_list in nodes:
            edges = edge_list.rstrip().split(':')
            parent = edges[0]
            children = edges[1].split(',')
            for child in children:
                if child == '':
                    continue
                parent_to_child.setdefault(parent, []).append(child)
                # contig_edge = tuple(sorted([parent[:-2],child[:-2]]))
                connected[(parent,child)] = True
    return parent_to_child,connected

def edge_summary_parse(edge_summary_file):
    parent_to_child = {}
    nodes = []
    with open(edge_summary_file) as file:
        header = get_header_dict(edge_summary_file)
        next(file)
        for line in file:
            line = split(line)
            parent,child = line[header["parent"]],line[header["child"]]
            nodes.append(parent)
            nodes.append(child)
            parent_to_child.setdefault(parent,[]).append(child)
            parent_to_child.setdefault(child, [])
    nodes = np.unique(nodes)
    return sorted(nodes),parent_to_child

def get_parent_children(line):
    edges = line.rstrip().split(":")
    parent,children = edges
    children = children.split(",")
    if "" in children: children.remove("")
    return (parent,children)

def get_connected(file):
    with open(file) as adj_list:
        connected = {}
        for line in adj_list:
            edges = line.rstrip().split(':')
            parent, children = edges
            children = children.split(',')
            for child in children:
                if child == '':
                    continue
                parent_contig = parent[:-2]
                child_contig = child[:-2]
                connected[(parent_contig,child_contig)] = True
    return connected

def get_contig_stats(file):
    contig_to_length = {}
    contig_to_cov = {}
    with open(file) as contig_table:
        next(contig_table)
        for line in contig_table:
            line = line.rstrip().split("\t")
            contig,cov,length = line
            contig_to_cov[contig] = float(cov)
            contig_to_length[contig] = int(length)

    return (contig_to_cov,contig_to_length)

def get_edge_stats(file):
    edge_to_weight = {}
    edge_to_cov = {}
    edge_to_length = {}
    parent_to_children = {}
    child_to_parents = {}
    start_nodes = []
    with open(file) as edge_summary:
        next(edge_summary)
        for edge in edge_summary:
            edgeStats = edge.rstrip().replace(" ", "\t").split("\t")
            parent, child, coverage, length, weight = edgeStats[:5]
            if parent[-2:] == "a+" and parent_to_children.get(parent) == None:
                start_nodes.append(parent)
            edge = (parent, child)
            edge_to_weight[edge] = float(weight)
            edge_to_cov[edge] = float(coverage)
            edge_to_length[edge] = float(length)
            parent_to_children.setdefault(parent,[]).append(child)
            child_to_parents.setdefault(child, []).append(parent)

    return edge_to_weight, edge_to_cov, edge_to_length, parent_to_children, child_to_parents, start_nodes

def get_edge_structures(file):
    edge_to_cov = {}
    parent_to_children = {}
    child_to_parents = {}
    start_edges = []
    edge_to_graph_type = {}
    edge_to_length = {}
    with open(file) as table:
        header = get_header_dict(file)
        next(table)
        for line in table:
            line = split(line)
            parent = line[header["parent"]]
            child = line[header["child"]]
            graph_edge = line[header["graph_edge"]]
            edge = (parent, child)
            edge_to_cov[edge] = float(line[header["coverage"]])
            parent_to_children.setdefault(parent,[]).append(child)
            child_to_parents.setdefault(child, []).append(parent)
            edge_to_graph_type[edge] = str(graph_edge)
            if edge_to_length.get(edge) == None:
                edge_to_length[edge] = int(line[header["length"]])
            if parent[-1:] == "-":
                start_edges.append(edge)
    return parent_to_children, child_to_parents, edge_to_cov, start_edges, edge_to_graph_type, edge_to_length


def reverse_complement(seq):
    comp = {"A":"T","T":"A","G":"C","C":"G"}
    comp_seq = ""
    for base in reversed(seq): comp_seq += comp[base]
    return comp_seq

def parse_fasta(fasta):
    with open(fasta) as fa:
        name_to_seq = {}
        seq_name = "start"
        sequence = ""
        for line in fa:
            line = line.rstrip()
            if len(line) < 1: continue
            if line[0] == ">":
                if seq_name != "start":
                    name_to_seq[seq_name] = sequence
                seq_name = line[1:]
                sequence = ""
            else:
                sequence += line
        if seq_name != "start":
            name_to_seq[seq_name] = sequence
    return name_to_seq

def get_popped_tandem(tandem_edge, edge_to_coverage, edge_to_length, max_size, min_coverage, max_fraction):
    tandem_edge_coverage = edge_to_coverage[tandem_edge]
    contig_edge = tuple(reversed(tandem_edge))
    contig_coverage = edge_to_coverage[contig_edge]
    contig_size = edge_to_length[contig_edge]
    if contig_coverage > min_coverage and contig_size < max_size and tandem_edge_coverage > min_coverage and tandem_edge_coverage < max_fraction * contig_coverage:
        return True
    else: return False

def get_simple_repeat(edge, edge_to_coverage, max_error, min_coverage, parent_to_children, children_to_parents):
    repeat_cov = edge_to_coverage[edge]
    if repeat_cov < min_coverage:
        return False
    looped_back_contigs = [] # list of tuples of insert contig and strand and bottleneck coverage of the loop
    expected_insert_cov = 0.5 * repeat_cov
    goal_node = edge[0]
    for child in parent_to_children.get(edge[1], []):
        putative_insert_contig_edge = (child, (child[:-1] + "+"))
        insert_children = parent_to_children.get(putative_insert_contig_edge[1], [])
        insert_parents = children_to_parents[child]
        if goal_node in insert_children and len(insert_children) == 1 and len(insert_parents) == 1:
            insert_cov = edge_to_coverage[putative_insert_contig_edge]
            if abs((insert_cov - expected_insert_cov)/expected_insert_cov) < max_error:
                looped_back_contigs.append((child, edge_to_coverage[(edge[1], child)]))
    if len(looped_back_contigs) == 0:
        return False
    sorted_bn_simple_repeat = sorted(looped_back_contigs, key=lambda tup: tup[1])
    return edge, sorted_bn_simple_repeat[0]

def get_inverted_repeat(edge, edge_to_coverage, max_error, min_coverage, parent_to_children, children_to_parents):
    repeat_cov = edge_to_coverage[edge]
    if repeat_cov < min_coverage:
        return False
    expected_insert_cov = 0.5 * repeat_cov
    goal_node = get_comp_edge(edge)[0]
    edge_children = parent_to_children.get(edge[1], [])
    if len(edge_children) == 2:
        if edge_children[0][:-2] == edge_children[1][:-2]:
            insert_children = ( (edge_children[0][:-1] + "+"), (edge_children[1][:-1] + "+"))
            insert_forw_children = parent_to_children.get(insert_children[0],[])
            insert_rev_children = parent_to_children.get(insert_children[1],[])
            if goal_node in insert_forw_children and goal_node in insert_rev_children and len(insert_forw_children) == 1 and len(insert_rev_children) == 1 :
                insert_cov = edge_to_coverage[(edge_children[0], insert_children[0])]
                if abs((insert_cov - expected_insert_cov) / expected_insert_cov) < max_error:
                    return edge, (edge_children[1], insert_cov)
    return False

def get_hypers(hyper_file):
    hyper_parents = set()
    hyper_parents_to_nodes = {}
    hyper_parents_to_edges = {}
    hyper_parents_to_cov = {}
    with open(hyper_file) as edges:
        header = get_header_dict(hyper_file)
        next(edges)
        for line in edges:
            line = split(line)
            parent = line[header["parent"]]
            cov = float(line[header["coverage"]])
            path = line[header["child"]].split(":")
            edges = [(first, second) for first, second in zip(path, path[1:])]
            hyper_parents.add(parent)
            hyper_parents_to_nodes[parent] = path
            hyper_parents_to_edges[parent] = edges
            hyper_parents_to_cov[parent] = cov
    return hyper_parents, hyper_parents_to_nodes, hyper_parents_to_edges, hyper_parents_to_cov

####################################
# Paired read table functions      #
####################################

def read_not_in_overlap(contig_len,read_start,read_stop, k, edit_distance=0, buffer=3):
    k_buffer = k + buffer
    if (read_start > k_buffer or read_stop > k_buffer) and (read_start < contig_len - k_buffer or read_stop < contig_len - k_buffer):
        if read_stop < k_buffer or read_start < k_buffer or read_start > contig_len - k_buffer or read_stop > contig_len - k_buffer:
            if edit_distance == 0:
                return True
            else:
                return False
        else: return True
    else:
        return False

def get_reads_far(read1_stop, contig1_length, read_strand1 , read2_stop, contig2_length, read_strand2, max_distance):
    distance = 0
    contig_stats = [(read1_stop, contig1_length, read_strand1), (read2_stop, contig2_length, read_strand2)]
    for stop,length,strand in contig_stats:
        if strand == 1:
            distance += length - stop
        else:
            distance += stop
    return distance >= max_distance

def get_node_suffix(start1,stop1,start2,stop2):
    strand1 = np.sign(start1 - stop1)
    strand2 = np.sign(start2 - stop2)
    if strand1 == strand2:
        if strand1 == 1:
            return ("a+", "b-")
        else: return ("b+","a-")
    else:
        if strand1 == 1:
            return ("a+","a-")
        else: return("b+","b-")

def get_comp_edge(edge):
    comp_orientation = {'a+': 'b-', 'b-': 'a+', 'a-': 'b+', 'b+': 'a-'}
    comp_child,comp_parent = edge
    comp_child_node = comp_child[:-2] + comp_orientation[comp_child[-2:]]
    comp_parent_node = comp_parent[:-2] + comp_orientation[comp_parent[-2:]]
    return (comp_parent_node,comp_child_node)

def get_insert(file):
    head = get_header_dict(file)
    inserts = []
    count = 0
    with open(file) as map:
        next(map)
        for line in map:
            line = split(line)
            if line[head['contig1']] == line[head['contig2']]:
                sorted_starts_coords = sorted([(int(line[head['back_coord1']]), int(line[head['strand1']])), (int(line[head['back_coord2']]),int(line[head['strand2']]))] , key=lambda start_strand: start_strand[0] )
                if sorted_starts_coords[0][1] == 1:
                    inserts.append(sorted_starts_coords[1][0] - sorted_starts_coords[0][0])
                    count += 1
            if count > 10000:
                break
    return inserts


####################################
# Cycle finder algorithm           #
####################################

def get_sorted_comp_path(path):
    comp_orientation = {'a+':'b-','b-':'a+','a-':'b+','b+':'a-'}
    comp_path = []
    for node in path:
        contig = node[:-2]
        comp = comp_orientation[node[-2:]]
        comp_path.append(contig+comp)
    return tuple(sorted(comp_path))

def get_renamed_contigs(contig_table, cycle_summary, ofile):
    header = get_header_dict(cycle_summary)
    contig_count = {} # count how many times each contig appears in cycles
    contig_to_name = {} # store renamed values
    line_count = 0
    repeat_num = 1
    with open(cycle_summary) as cyc_sum:
        next(cyc_sum)
        for line in cyc_sum:
            line = split(line)
            line_count += 1
            if line_count % 2 != 0:
                contig = line[header["edge"]].split(">")[0][:-2]
                count = contig_count.setdefault(contig,0) + 1
                contig_count[contig] = count
                if count != 1 and contig_to_name[contig][0] != "R": # if this contig is a repeat and hasn't been labeled as such yet
                    contig_to_name[contig] = "R" + str(repeat_num)
                    repeat_num += 1
                elif count == 1: # if this contig is, at this point in the contig cycle table, unique to the cycle
                    cycle = line[header["cycle"]][3:]
                    edge_count = (int(line[header["edge_count"]]) + 1) / 2
                    contig_to_name[contig] = cycle + "_" + str(int(edge_count))
    contigs = list(get_contig_stats(contig_table)[1].keys())
    ofile = open(ofile, "w+")
    ofile.write(write_line("original","cycle_name"))
    for contig in contigs:
        rename = contig_to_name.get(contig)
        if rename == None:
            ofile.write(write_line(contig,contig))
        else:
            ofile.write(write_line(contig,rename))

    return contig_to_name

def get_cycle_specs(cyc_table):
    con_to_cyc = {}
    putative_cycles = [] # list of tuples of (cyc, cyc_length)
    cycle_to_contigs = {} # key: cycles value: list of contigs that belong to that cycle
    with open(cyc_table) as cycs:
        header = get_header_dict(cyc_table)
        next(cycs)
        for line in cycs:
            line = split(line)
            contig = line[header["original_contig"]]
            cycle_stats = con_to_cyc.setdefault(contig ,[])
            cycle = line[header["cycle"]]
            cumsum = int(float(line[header["cum_sum"]]))
            contig_len = int(float(line[header["contig_len"]]))
            orientation = line[header["orientation"]]
            cycle_length = int(float(line[header["cycle_length"]]))
            cycle_stats.append((cycle, cumsum, contig_len, orientation, cycle_length))
            con_to_cyc[contig] = cycle_stats
            cycle_to_contigs.setdefault(cycle,[]).append(contig)
            if (cycle, cycle_length) not in putative_cycles:
                putative_cycles.append((cycle,cycle_length))

    return con_to_cyc, putative_cycles, cycle_to_contigs

""" 
returns: -tuple of cycles that read1, read2 belong to (None if a read doesn't belong to a cycle)
"""
def contig_to_cycs(contig1, contig2, cycle_members, strand1=None, strand2=None):
    orientation_to_strand = {("+", 1) : 1, ("+", -1) : -1, ("-", 1) : -1, ("-", -1) : 1}
    same_strand_cycle_stats = []
    cycles_to_cycle_stats = {}
    if contig2 != "singleton":
        cycles1 = cycle_members.get(contig1,[ tuple( [None]*5 ) ]) # positions of contig1 in any cycles
        cycles2 = cycle_members.get(contig2,[ tuple( [None]*5 ) ]) # positions of contig2 in any cycles
        if contig1 == contig2: # reads are intracontig
            return [tuple(zip(cycles1[i],cycles2[i])) for i in range(len(cycles1))]
        else: # both contigs belong to different cycles, different number of cycles return in list 1 and 2, ...
            read_cycles = []
            for cycle1, cumsum1, c1_len, ori1, cyc1_len in cycles1: # the loops exist in the case of a repeat cycle contig
                for cycle2, cumsum2, c2_len, ori2,cyc2_len in cycles2:
                    cycle_stats = tuple([(cycle1,cycle2) , (cumsum1, cumsum2), (c1_len, c2_len), (ori1, ori2), (cyc1_len,cyc2_len)])
                    # if ori1 != None and ori2 != None:
                    #     if orientation_to_strand[(ori1, strand1)] == orientation_to_strand[(ori2, strand2)]:
                    #         same_strand_cycle_stats.append(cycle_stats)
                    cycles_to_cycle_stats.setdefault((cycle1, cycle2), []).append(cycle_stats)
                    read_cycles.append(cycle_stats)
            # for cycle_stat in same_strand_cycle_stats: # take out flipped repeats (a repeat that is one direction then the other direction)
            #     ((cycle1, cycle2), (cumsum1, cumsum2), (c1_len, c2_len), (ori1, ori2), (cyc1_len, cyc2_len)) = cycle_stat
            #     if cycle1 == cycle2:
            #         better_cycles = [1 for ((cycle1, cycle2), (cumsum1, cumsum2), (c1_len, c2_len), (ori1, ori2), (cyc1_len, cyc2_len)) in cycles_to_cycle_stats[(cycle1, cycle2)] if orientation_to_strand[(ori1, strand1)] != orientation_to_strand[(ori2, strand2)]]
            #         if len(better_cycles) == 1:
            #             read_cycles.remove(cycle_stat)
            return read_cycles
    else:
        return cycle_members[contig1]


"""
returns: -the cycle relative coordinates of a read that maps to a cycle in a paired read that supports a cycle 
"""
def get_cycle_coords(cycle_contig_start, contig_length, contig_orientation, read_start, read_stop):
    if cycle_contig_start == None: return None # if read doesn't map to a cycle
    if contig_orientation == "+": # if the cycle and contig are pointing the same way
        start = cycle_contig_start + read_start
        stop = cycle_contig_start + read_stop
    else:
        start = cycle_contig_start + (contig_length - read_start)
        stop = cycle_contig_start + (contig_length - read_stop)
    strand = np.sign(start - stop)
    return stop, strand


def get_edge_type(file):
    edge_to_type = {}
    with open(file) as edge_summary:
        next(edge_summary)
        for edge in edge_summary:
            edgeStats = edge.rstrip().replace(" ", "\t").split("\t")
            parent, child, type = edgeStats[0], edgeStats[1], edgeStats[-1]
            edge = (parent, child)
            edge_to_type[edge] = type
    return edge_to_type


def get_cyc_contig_stats(file):
    cyc_to_length = {}
    cyc_to_contig = {}
    contig_to_length = {}
    with open(file) as contigs:
        header = get_header_dict(file)
        next(contigs)
        for line in contigs:
            line_s = split(line)
            cycle = line_s[header["cycle"]]
            contig = line_s[header["original_contig"]]
            contig_to_length[contig] = int(line_s[header["contig_len"]])
            cyc_to_contig.setdefault(cycle, []).append(contig)
            cyc_to_length[cycle] = float(line_s[header["cycle_length"]])
    return cyc_to_length, cyc_to_contig, contig_to_length

def get_cycle_scores(file):
    cycle_scores = {}
    header = get_header_dict(file)
    with open(file) as cycles:
        next(cycles)
        for line in cycles:
            line = split(line)
            if float(line[header["max_nonsupport"]]) != 0:
                score = float(line[header["bottleneck"]]) / float(line[header["max_nonsupport"]])
            else:
                score = 100
            cycle = line[header["cycle"]]
            cycle_scores[cycle] = score
    return cycle_scores

####################################
# FUNCTIONAL CLASSIFICATION        #
####################################

def is_int(string):
    try:
        int(string)
        return True
    except ValueError:
        return False

####################################
# TIMEPOINT                        #
####################################

def read_num(file):
    with open(file) as read_stats:
        next(read_stats)
        return int(split(read_stats.readline())[1])

def isnt_cyc(seq1, seq2):
    return not (seq1[:3] == "CYC" or seq2[:3] == "CYC")

def is_support(seq1, seq2):
    return seq1 == seq2

def get_cycle_specs_tp(cyc_table):
    cycle_junction_coords = {}
    putative_cycles = [] # list of tuples of (cyc, cyc_length)
    cycle_to_contigs = {} # key: cycles value: list of contigs that belong to that cycle
    with open(cyc_table) as cycs:
        header = get_header_dict(cyc_table)
        next(cycs)
        for line in cycs:
            line = split(line)
            contig = line[header["original_contig"]]
            cycle = line[header["cycle"]]
            cycle_junction_coords.setdefault(cycle, []).append(float(line[header["cum_sum"]]))
            cycle_length = int(float(line[header["cycle_length"]]))
            cycle_to_contigs.setdefault(cycle,[]).append(contig)
            putative_cycles.append((cycle,cycle_length))
    return cycle_junction_coords, [(cyc, int(length)) for cyc, length in np.unique(putative_cycles, axis = 0) ], cycle_to_contigs

def cyc_read_usable(start, stop, cyc_coords, k):
    if cyc_coords == None:
        return True
    for coord in cyc_coords:
        if start > stop:
            if start < coord + k + 2 and stop > coord - 2:
                return False
        else:
            if stop < coord + k + 2 and start > coord - 2:
                return False
    return True