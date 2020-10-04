import argparse
import shutil
import os

###############################################################################################
# Purpose: clean directories                                                                  #
###############################################################################################

parser = argparse.ArgumentParser(description='write the fasta sequences of putative cycles')
parser.add_argument('idir', metavar='<contig fasta input>', type=str, help='contig fasta input')

args = parser.parse_args()

raw_dir = os.path.join(args.idir, "raw")
shutil.rmtree(raw_dir)
if not os.path.exists(raw_dir): os.mkdir(raw_dir)
files = [f for f in os.listdir(args.idir) if os.path.isfile(os.path.join(args.idir, f))]
for f in files: shutil.move(os.path.join(args.idir, f), raw_dir)
shutil.copy(os.path.join(args.idir, "dominant_cycles", "cycle_contig_table"), args.idir)
shutil.copy(os.path.join(args.idir, "dominant_cycles", "cycle_stats"), args.idir)
shutil.copy(os.path.join(args.idir, "dominant_cycles", "cycles.fasta"), args.idir)
shutil.copy(os.path.join(raw_dir, "cycle_covs_long"), args.idir)
shutil.move(os.path.join(args.idir, "dominant_cycles"), raw_dir)