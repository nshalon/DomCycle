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

`$ python DomCycle.py -g `*`FASTG_PATH`*` -1 `*`R1.sam`*` -2 `*`R2.sam`*` -k `*`kmer_size`*

To run an example:

`$ python DomCycle.py -g ex/k77.fastg -1 ex/R1_map.sam -2 ex/R2_map.sam -k 77`

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
- numpy
- argparse
- pathlib
