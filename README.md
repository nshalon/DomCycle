# DomCycle: *de novo* identification of circular genetic elements in microbial and metagenomic samples

DomCycle relies on coverage profiles of contigs in an assembly graph to confidently infer circular genetic elements in a sequenced sample. First, through an inputted FASTG file and SAM mapped read files, DomCycle constructs an assembly graph. In the assembly graph, nodes represent contig ends and edges represent either (1) connections between contig ends supported by read pairs or (2) a 3'-5' connection of a single contig. Each edge has a coverage value. Using the final assembly graph and contig-level coverage values, the dominant cycles algorithm (Algorithm 1 in the manuscript) is used to identify putative circular elements. Each element is further vetted on a nucleotide-level coverage profile for final cycle output.

To read more about DomCycle and its performance, please refer to our [preprint manuscript](https://www.biorxiv.org/content/10.1101/2021.05.25.445656v1). 

## Requirements and Tested Software

Python 3.x, perl, and bash are required to run DomCycle. Please see requirements.txt to install python dependencies.

DomCycle was tested with the following tools and specifications:

Assembly: MEGAHIT (k=77), Read Mapping: BWA

Python dependencies:
- (python version 3.x)
- numpy
- argparse
- pathlib
- scipy

## Installation

To access this pipeline:

`$ git clone https://github.com/nshalon/DomCycle.git`

To install python dependencies:

`$ python3 -m pip install -r DomCycle/requirements.txt`

## Input

DomCycle requires the following inputs:

- -g: assembly graph in FASTG format
- -1: R1 mapped to the assembly FASTA in SAM format
- -2: R2 mapped to the assembly FASTA in SAM format
- -k: kmer used for assembly

Contig names in the FASTG file must match the contig names in the SAM files. When using MEGAHIT an issue occurs if the FASTG was generated using the megahit_toolkit (see [here](https://github.com/voutcn/megahit/wiki/Visualizing-MEGAHIT's-contig-graph)), and the reads were mapped directly to the FASTA file. In that case there is a mismatch between contig names in the FASTG and SAM files. To correct the contig names in the FASTG file we supply a script (pl/rename_fastg.pl) that recieves as input a FASTA file and an associated FASTG file, and produces a corrected FASTG file that is suitable as input for DomCycle, alongside a naming mapping file. Usage example:
`$ perl pl/rename_fastg.pl k77.fasta original_k77.fastg k77.fastg mapping_file`

To create the input SAM files, please map each read side separately to the FASTA e.g. `bwa mem contigs_filename.fa R1.fastq > R1_map.sam`. Additionally, please ensure that the SAM files contain the full, untruncated contig names.

Configurable parameters:

- -r: average read length in library [DEFAULT: 150]
- -o: output directory [DEFAULT: ./output]
- -a: alpha, out-coverage dampening factor [DEFAULT: 0]
- -p: pval threshold for identifying dominant cycles [DEFAULT: 0.01]
- --minscore: minimum score for classifying dominant cycles [DEFAULT: 1]
- --sample: custom sample name [DEFAULT: "metagenome"]
- --min_quality: minimum mapping quality score when filtering reads [DEFAULT: 0]
- --min_match_len: minimum read mapping match length when filtering reads [DEFAULT: 50]
- --max_edit_distance: maximum edit distance when filtering reads [DEFAULT: 1]

## Running DomCycle

To run DomCycle:

`$ cd DomCycle`

`$ python3 DomCycle.py -g `*`FASTG_PATH`*` -1 `*`R1.sam`*` -2 `*`R2.sam`*` -k `*`kmer_size`*

To run an example:

`$ python3 DomCycle.py -g ex/k77.fastg -1 ex/R1_map.sam -2 ex/R2_map.sam -k 77`

## Output

Primary output files and associate fields:
1. cycle_contig_table: a table where each row represents one contig in a dominant cycle. Each dominant cycle will have an ordered list of contigs that compose it, where two successive contigs in the table belonging to the same cycle have an edge between them and there's an additional edge between the last contig belonging to the cycle and the first contig belonging to the cycle in the table.
    - cycle: the dominant cycle containing the contig
    - FASTG_edge: whether the edge originated from the FASTG
    - cum_sum: the contig's start position in the cycle sequence space
    - orientation: denotes the relative orientation of the two contigs composing the cycle edge. (+) indicates that the two contigs are connected via the same reference strand and (-) indicates that they're connected via opposite reference strands
    - contig_index: indicates the contig is the ith contig in the cycle
2. dominant_cycle_stats: displays the prominent cycle coverage statistics that follow the manuscript's methods. "dominant_cycle_stats" refers to only the calculated dominant cycles.
    - bottleneck: the minimum support x-coverage across all of the cycle bases
    - median_support: the cycle's median support x-coverage
    - avg external: the sum of avg_inter and avg_intra-nonsupport 
    - avg_inter: for each cycle strand, the total number of reads where one side maps to the cycle and the other side maps elsewhere, then averaged over both cycle strands
    - avg_intra-nonsupport: for each cycle strand, total number of reads where both sides map to the same cycle either far away from each other or to the same cycle strand, then averaged over both cycle strands
    - avg_singleton: the total number of reads where one side maps to the cycle and the other side doesn't map anywhere after filtering out reads according to the mapping quality thresholds set
    - pval: given a null hypothesis of bottleneck = avg_external, the p-value associated with the associated statistics under D ~ binomial(n = # reads in sample, p = (avg_external) / (# reads in sample))
    - score: the bottleneck / avg_external
3. cycles.fasta: a .fasta file of the dominant cycles after combining the contigs that compose each cycle
4. cycle_covs_long: basepair resolution of coverage statistics for each cycle, where each row in the table represents a base in the cycle and the coverages follow the above definitions (but displayed separately for either cycle strand). You can remove this file to save disk space.

Other intermediate files produced during the run from mapping, graph building, cycle identification (dominant and otherwise) can be found in the folders "raw" and "dominant_cycles_pre" under the output directory.
