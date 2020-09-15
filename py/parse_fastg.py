import argparse
import numpy as np
import re
from operator import itemgetter

#############################################################################
# Purpose: input a fastg and
# 1) Prepare adjacency list format, in which header of each list is a strand (a or b) 
# and direction (+ or -). In this parsing, any contig that contains a ' means that if it's
# in the first part of fastg on the "b" strand (b+) and if on latter part of fastg header it's
# on the "a" strand (a-)
# 2) That means that for each contig, it has four nodes: a+,a-,b+,b-
# 3) that means that the adjacency list for "C1:C2',C3,C5';" would look like this:
# C1a+: C2b-,C3a-,C5b-
# 4) The adjacency list will always have two fixed edges, with the directionality of a- -> a+ ; a+ -> neighbor
# and b- -> b+ ; b+ 
#############################################################################

parser = argparse.ArgumentParser(description='rename reads to have FW and RV')
parser.add_argument('ifn_fastg', metavar='<file containing just the headers of each connection in fastg>', type=str, help='headers means everything after ">" and no bases')
parser.add_argument('ofn_adjacency_list', metavar='<output text representation of adjacency list>', type=str, help='adjacency list with format PARENT:CHILD,CHILD')
parser.add_argument('ofn_rename', metavar='<file that shows map of original name to rename>', type=str, help='output full file NODE C1')
parser.add_argument('ofn_renamed_fasta', metavar='<fasta that has the renamed CXXX headers>', type=str, help='output renamed fasta')
parser.add_argument('ofn_adjacency_matrix', metavar='<output adjacency matrix>', type=str, help='output adjacency matrix, with rows as parent, 1 as connection, 0 as no connection')

args = parser.parse_args()

# parse the headers and create a list of all contigs in fastg
# then map from each contig to its length and rename contigs
# according to increasing lengths

contig_to_length = {}
original_contig_to_rename = {}
adjacency_list = {}
contig_to_sequence = {}

with open(args.ifn_fastg) as fastg:
    line_count = 0
    for line in fastg:
        if line_count % 4 == 0:
            line = line.replace('>','')
            line = line.replace(';',':')
            line = line.split(':')
            contig = line[0]
            contig = contig.replace("'","")
            contig = contig.replace(" ","")
            length_contig = int(contig.split("_")[3])
            contig_to_length[contig] = length_contig
        elif line_count % 4 == 1:
            sequence = line.rstrip()
            contig_to_sequence[contig] = sequence
        
        if line_count % 10000 == 0:
            print("fastg line count:",line_count)
        line_count += 1
    #now convert hash table to list of tuples and sort by value
    contig_to_length_list = [(k,v) for k,v in contig_to_length.items()]
    sorted_contig_to_length_list = sorted(contig_to_length_list,key=itemgetter(1))
    length_order = 1
    rename_to_original_contig = {}
    #creates mapping of original fastg contig name to renamed order
    #and vice versa
    ofile = open(args.ofn_rename,'w+')
    sorted_contig_to_length_list.reverse()
    ofasta = open(args.ofn_renamed_fasta,'w+')
    for contig,length in sorted_contig_to_length_list:
        rename = "C" + str(length_order)
        original_contig_to_rename[contig] = rename
        rename_to_original_contig[rename] = contig
        sequence = contig_to_sequence[contig]
        
        ofasta.write(">" + rename + "\n" + sequence + "\n")
        
        ofile.write(contig + " " + rename + " " + str(contig_to_length[contig]) + "\n")
        length_order += 1

    for rename in rename_to_original_contig.keys():
        contig_nodes = [rename + 'a+', rename + 'a-', rename + 'b+', rename + 'b-']
        for node in contig_nodes:
            adjacency_list[node] = []
            # always the '-' will be connected to '+' and create an edge
            if '-' in node:
                adjacency_list[node].append(contig_nodes[contig_nodes.index(node)-1])

print("Now writing adjacency lists and matrix")

with open(args.ifn_fastg) as fastg:
    line_count = -1
    for line in fastg:
        line_count += 1
        if line_count % 2 != 0: #make sure only right headers are read
            continue
        line = line.replace(">","").replace(",",":").replace(";","").rstrip().split(":")
        parent = line[0]
        children = line[1:]
        contig_node = ""
        if "'" == parent[-1]:
            contig_rename = original_contig_to_rename[parent[:-1]]
            contig_node = contig_rename + 'b+'
        else:
            contig_rename = original_contig_to_rename[parent]
            contig_node = contig_rename + 'a+'
        for child in children:
            if "'" == child[-1]:
                contig_rename = original_contig_to_rename[child[:-1]]
                child_node = contig_rename + 'b-'
                adjacency_list[contig_node].append(child_node)
            else:
                contig_rename = original_contig_to_rename[child]
                child_node = contig_rename + 'a-'
                adjacency_list[contig_node].append(child_node)
                
    adjacency_list_tuple = [(k,v) for k,v in adjacency_list.items()]
    sorted_adjacency_list_tuple = sorted(adjacency_list_tuple,key=itemgetter(0))
    ofile = open(args.ofn_adjacency_list,'w+')
    contig_list = []
    item_count = 0
    for item in sorted_adjacency_list_tuple:
        item_count += 1
        children = ""
        for child in item[1]:
            children += child + ","
        children = children[:-1]
        ofile.write(item[0] + ':' + children + '\n')
        contig_list.append(item[0])
        
    # omatrix = open(args.ofn_adjacency_matrix,"w+")
    # omatrix.write("\t")
    # for contig in contig_list:
    #     omatrix.write(contig + "\t")
    # omatrix.write("\n")
    # item_count = 0
    # for contig in contig_list:
    #     if item_count % 10000 == 0:
    #         print("contig number",item_count)
    #     item_count += 1
    #     omatrix.write(contig + "\t")
    #     adj_list_index = contig_list.index(contig)
    #     contig_children = sorted_adjacency_list_tuple[adj_list_index][1]
    #     for contig in contig_list:
    #         if contig in contig_children:
    #             omatrix.write("1")
    #         else:
    #             omatrix.write("0")
    #         omatrix.write("\t")
    #     omatrix.write("\n")
            
            
    
        
        



        
