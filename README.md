# DomCycle

DomCycle is a tool for *de novo* recognition of circular mobile genetic elements, such as plasmids and phages, in microbial metagenomic samples.

To access this pipeline:

`$ git clone https://github.com/nshalon/DomCycle.git`

DomCycle requires the following inputs:

- -g: assembly graph in FASTG format
- -1: R1 mapped to the assembly in .sam format
- -2: R2 mapped to the assembly in .sam format
- -k: kmer used for assembly

To run DomCycle:

`$ cd DomCycle`

`$ python3 DomCycle.py -g `*`FASTG_PATH`*` -1 `*`R1.sam`*` -2 `*`R2.sam`*` -k `*`kmer_size`*

To run an example:

`$ python3 DomCycle.py -g ex/k77.fastg -1 ex/R1_map.sam -2 ex/R2_map.sam -k 77`

Other configurable parameters:

- -r: average read length in library [DEFAULT: 150]
- -o: output directory [DEFAULT: ./output]
- -a: alpha, out-coverage dampening factor [DEFAULT: 0]
- -p: pval threshold for identifying dominant cycles [DEFAULT: 0.01]
- --minscore: minimum score for classifying dominant cycles [DEFAULT: 1]
- --sample: custom sample name [DEFAULT: "metagenome"]
- --min_quality: minimum mapping quality score when filtering reads [DEFAULT: 0]
- --min_match_len: minimum read mapping match length when filtering reads [DEFAULT: 50]
- --max_edit_distance: maximum edit distance when filtering reads [DEFAULT: 1]

DomCycle was tested with the following tools and specifications:

Assembly: MEGAHIT v1.1.3
Mapping: BWA V0.7.12
Python v3.5.1
Read Trim Length: 50bp

Python dependencies:
- (python version 3.x)
- numpy
- argparse
- pathlib
- scipy

Main output files and main column descriptions:
1. cycle_contig_table: a table where each row represents one contig in a dominant cycle. Each dominant cycle will have an ordered list of contigs that compose it, where two successive contigs in the table belonging to the same cycle have an edge between them and there's an additional edge between the last contig belonging to the cycle and the first contig belonging to the cycle in the table.
    - cycle: the dominant cycle containing the contig
    - FASTG_edge: whether the edge originated from the FASTG
    - cum_sum: the contig's start position in the cycle sequence space
    - orientation: denotes the relative orientation of the two contigs composing the cycle edge. (+) indicates that the two contigs are connected via the same reference strand and (-) indicates that they're connected via opposite reference strands
    - contig_index: indicates the contig is the ith contig in the cycle
2. cycle_stats: displays the prominent cycle coverage statistics that follow the manuscript's methods.
    - bottleneck: the minimum support x-coverage across all of the cycle bases
    - median_support: the cycle's median support x-coverage
    - avg external: the sum of avg_inter and avg_intra-nonsupport 
    - avg_inter: for each cycle strand, the total number of reads where one side maps to the cycle and the other side maps elsewhere, then averaged over both cycle strands
    - avg_intra-nonsupport: for each cycle strand, total number of reads where both sides map to the same cycle either far away from each other or to the same cycle strand, then averaged over both cycle strands
    - avg_singleton: the total number of reads where one side maps to the cycle and the other side doesn't map anywhere after filtering out reads according to the mapping quality thresholds set
    - pval: given a null hypothesis of bottleneck = avg_external, the p-value associated with the associated statistics under D ~ binomial(n = # reads in sample, p = (avg_external) / (# reads in sample))
    - score: the bottleneck / avg_external
    - lower_bound_cov: the bottleneck - avg_external
3. cycle_cov_long: basepair resolution of coverage statistics for each cycle, where each row in the table represents a base in the cycle and the coverages follow the above definitions (but displayed separately for either cycle strand)
4. cycles.fasta: a .fasta file of the dominant cycles after combining the contigs that compose each cycle

Other files produced from mapping, graph building, cycle identification (dominant and otherwise) can be found in the "raw" subdirectory under the specified output directory.
