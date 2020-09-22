import argparse
import os
import subprocess
import datetime
import pathlib


def run_steps():
    parser = argparse.ArgumentParser(description='Run DomCycle')
    parser.add_argument('-g', "--graph", type=str, help='Path to FASTG file created from assembly')
    parser.add_argument('-1', "--r1", metavar='R1 MAPPED TO ASSEMBLY', type=str, help='R1 mapped to assembly in .sam format using BWA output')
    parser.add_argument('-2', "--r2", metavar='R2 MAPPED TO ASSEMBLY', type=str, help='R2 mapped to assembly in .sam format using BWA output')
    parser.add_argument('-k', "--k", metavar='KMER SIZE', type=int, help='kmer size used for assembly')
    parser.add_argument('-r', "--read_len", metavar='AVG (PRE-TRIMMED) READ LENGTH', type=float, help='input the average read length in the input library', default=150)
    parser.add_argument('-o', "--outdir", metavar='OUTPUT_DIR', type=str, help='Parent output directory', default="/".join([os.getcwd(),"output"]))
    parser.add_argument('-a', "--alpha", metavar='OUT COVERAGE DAMPENENING FACTOR', type=float, help='Out coverage dampening factor for finding cycles in graph', default=0)
    parser.add_argument('-p', "--maxpval", metavar='P-VALUE', type=float, help='Minimum p-value used to classify dominant cycles', default=0.01)
    parser.add_argument("--minscore", metavar='MIN SCORE', type=float, help='Minimum score for a dominant cycle', default=1)
    parser.add_argument('--sample', metavar='SAMPLE NAME', type=str, help='your custom sample name', default="metagenome")
    parser.add_argument('--min_quality', metavar='MIN QUALITY THRESHOLD', type=int, help='min map quality threshold', default=0)
    parser.add_argument('--min_match_len', metavar='MIN MATCH READ LENGTH', type=int, help='minimum match read length for filtering mapped reads', default=50)
    parser.add_argument('--max_edit_distance', metavar='MAX MISMATCH', type=int, help='maximum mismatches tolerated for filtering mapped reads', default=1)
    parser.add_argument('--steps', metavar='STEPS TO RUN', type=str, help='maximum mismatches tolerated for filtering mapped reads', default="1-10")

    args = parser.parse_args()

    script_dir = str(pathlib.Path(__file__).parent.absolute())
    steps = parse_steps(args.steps)
    odir = args.outdir
    py_dir = os.path.join(script_dir, "py")
    pl_dir = os.path.join(script_dir, "pl")

    if not os.path.exists(odir):
        os.mkdir(odir)
    else:
        print("WARNING!\nOutput directory already exists - consider changing the directory to avoid overwriting any previous results.")

    if args.outdir[0] != "/":
        args.outdir = os.path.join(os.getcwd(), args.outdir)
    if args.graph[0] != "/":
        args.graph = os.path.join(os.getcwd(), args.graph)
    if args.r1[0] != "/":
        args.r1 = os.path.join(os.getcwd(), args.r1)
    if args.r2[0] != "/":
        args.r2 = os.path.join(os.getcwd(), args.r2)

    if 1 in steps:
        print("\n\nParsing assembly graph...")
        subprocess.check_output(["python", "parse_fastg.py", args.graph, os.path.join(odir, "adjacency_list"),
                         os.path.join(odir, "contig_rename_map"), os.path.join(odir, "renamed_final_contigs.fa"),
                         os.path.join(odir, "adjacency_matrix")], cwd=py_dir)
        # except output.CalledProcessError as e:
        #     print("Error occurred!")
        #     print(e)
        #     return True

    if 2 in steps:
        print("\n\nParsing R1...")
        subprocess.check_output(["perl", "parse_bwa_sam.pl", args.r1, os.path.join(odir, "R1table"),
                         os.path.join(odir, ".R1table_stats")], cwd=pl_dir)
        print("Parsing R2...")
        subprocess.check_output(["perl", "parse_bwa_sam.pl", args.r2, os.path.join(odir, "R2table"),
                         os.path.join(odir, ".R2table_stats")], cwd=pl_dir)
    if 3 in steps:
        print("\n\nFiltering reads: min quality", args.min_quality, "; min match length", args.min_match_len, "; max edit distance", args.max_edit_distance)
        print("Filtering R1...")
        subprocess.check_output(["perl", "filter_map.pl", os.path.join(odir, "R1table"), str(args.min_quality),
                         str(args.min_match_len), str(args.max_edit_distance), os.path.join(odir, "filter_R1table"),
                         os.path.join(odir, ".filter_R1table_stats")], cwd=pl_dir)
        print("Filtering R2...")
        subprocess.check_output(["perl", "filter_map.pl", os.path.join(odir, "R2table"), str(args.min_quality),
                         str(args.min_match_len), str(args.max_edit_distance), os.path.join(odir, "filter_R2table"),
                         os.path.join(odir, ".filter_R2table_stats")], cwd=pl_dir)

    print("\n\nBUILDING GRAPH...")

    if 4 in steps:
        print("\n\nCalculating contig (internal edge) coverages...")
        subprocess.check_output(["python", "contig_cov_core.py", os.path.join(odir, "renamed_final_contigs.fa"),
                         os.path.join(odir, "filter_R1table"), os.path.join(odir, "filter_R2table"),
                         str(args.read_len), os.path.join(odir, "contig_table")], cwd=py_dir)

    if 5 in steps:
        print("\n\nPairing reads and calculating read statistics...")
        subprocess.check_output(["python", "pair_reads.py", os.path.join(odir, "filter_R1table"),
                         os.path.join(odir, "filter_R2table"), os.path.join(odir, "filter_singleton_table"),
                         os.path.join(odir, "filter_paired_table")], cwd=py_dir)
        subprocess.check_output(["python", "read_stats_core.py", os.path.join(odir, "filter_paired_table"),
                         str(args.k), str(args.read_len), os.path.join(odir, "read_stats.txt")], cwd=py_dir)

    if 6 in steps:
        print("\n\nFinding all external edges...")
        subprocess.check_output(["python", "find_new_variable_edges.py", os.path.join(odir, "filter_paired_table"),
                         os.path.join(odir, "adjacency_list"), os.path.join(odir, "contig_table"),
                         os.path.join(odir, "read_stats.txt"), os.path.join(odir, "edge_summary")], cwd=py_dir)

    print("\n\nGRAPH BUILT!")

    if 7 in steps:
        print("\n\nFinding cycles in the graph...")
        subprocess.check_output(["python", "dominant_cycles.py", os.path.join(odir, "edge_summary"),
                         os.path.join(odir, "read_stats.txt"), str(args.alpha),
                         os.path.join(odir, "cycle_contig_table")], cwd=py_dir)

    if 8 in steps:
        print("\n\nCalculating cycle coverage statistics in the cycle space...")
        subprocess.check_output(["python", "cycle_coverages_2.py", args.sample, os.path.join(odir, "cycle_contig_table"),
                         os.path.join(odir, "filter_paired_table"), os.path.join(odir, "filter_singleton_table"),
                         os.path.join(odir, "contig_table"), os.path.join(odir, "read_stats.txt"),
                         os.path.join(odir, "cycle_covs_long"), os.path.join(odir, "cycle_cov_summary"),
                         os.path.join(odir, "out_contigs")], cwd=py_dir)

    if 9 in steps:
        print("\n\nCreating cycle fastas...")
        subprocess.check_output(["python", "cycle_fastas.py", os.path.join(odir, "cycle_contig_table"),
                         os.path.join(odir, "renamed_final_contigs.fa"), os.path.join(odir, "read_stats.txt"),
                         os.path.join(odir, "cycles.fasta")], cwd=py_dir)

    if 10 in steps:
        if not os.path.exists(os.path.join(odir, "dominant_cycles")):
            os.mkdir(os.path.join(odir, "dominant_cycles"))
        print("\n\nIdentifying dominant cycles...\n")
        subprocess.check_output(["python", "extract_p.py", os.path.join(odir, "cycle_cov_summary"),
                         os.path.join(odir, "cycles.fasta"), os.path.join(odir, "cycle_contig_table"),
                         str(args.maxpval), "0", "0", str(args.minscore), os.path.join(odir, "cycle_stats"), os.path.join(odir, "dominant_cycles", "cycle_stats"),
                         os.path.join(odir, "dominant_cycles", "cycles.fasta"), os.path.join(odir, "dominant_cycles", "cycle_contig_table")],
                         cwd=py_dir)

        print(len(open(os.path.join(odir, "dominant_cycles", "cycle_stats")).readlines()) -1, "dominant cycles found!")

def parse_steps(steps):
    if "-" in steps: # sequential list of steps
        start_step = int(steps.split("-")[0])
        stop_step = int(steps.split("-")[1])
        return [(step + start_step) for step in range(stop_step)]
    elif "," in steps: # list of steps
        return [int(step) for step in steps.split(",")]
    else: # singular step
        return int(steps)

if __name__ == "__main__":
    print("Running DomCycle!")
    print("Start time:", datetime.datetime.now())
    error = run_steps()
    if not error:
        print("\n\nFinished running DomCycle!")
    print("Stop time:", datetime.datetime.now())