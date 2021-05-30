import argparse
import os
import subprocess
import datetime
import pathlib

def run_steps() -> None:
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
    parser.add_argument('--steps', metavar='STEPS TO RUN', type=str, help='maximum mismatches tolerated for filtering mapped reads', default="1-13")

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

    if not os.path.exists(os.path.join(odir, "dominant_cycles_pre")):
        os.mkdir(os.path.join(odir, "dominant_cycles_pre"))
    if not os.path.exists(os.path.join(odir, "dominant_cycles")):
        os.mkdir(os.path.join(odir, "dominant_cycles"))

    if args.outdir[0] != "/":
        args.outdir = os.path.join(os.getcwd(), args.outdir)
    if args.graph[0] != "/":
        args.graph = os.path.join(os.getcwd(), args.graph)
    if args.r1[0] != "/":
        args.r1 = os.path.join(os.getcwd(), args.r1)
    if args.r2[0] != "/":
        args.r2 = os.path.join(os.getcwd(), args.r2)

    cmds = get_step_dicts()
    for step in steps:
        cmd, messages = cmds[step]
        script_dir = py_dir if cmd[0] == "python3" else pl_dir
        for message in messages:
            print(message)
        subprocess.check_output(cmd, cwd=script_dir)

def get_step_dicts() -> dict[int, tuple[list[str], list[str]]]:
    """
    Retrieve a dict with pipeline steps and messages
    """
    cmds = {
        1 : (["python3", "parse_fastg.py", args.graph, os.path.join(odir, "adjacency_list"),
                         os.path.join(odir, "contig_rename_map"), os.path.join(odir, "renamed_final_contigs.fa"),
                         os.path.join(odir, "adjacency_matrix")], ["\n\nParsing assembly graph..."]),
        2 : (["perl", "parse_bwa_sam.pl", args.r1, os.path.join(odir, "R1table"),
                         os.path.join(odir, ".R1table_stats")], ["\n\nParsing R1..."]),
        3 : (["perl", "parse_bwa_sam.pl", args.r2, os.path.join(odir, "R2table"),
                         os.path.join(odir, ".R2table_stats")], ["\n\nParsing R2..."]),
        4 : (["python3", "contig_cov_core.py", os.path.join(odir, "renamed_final_contigs.fa"),
                         os.path.join(odir, "filter_R1table"), os.path.join(odir, "filter_R2table"),
                         str(args.read_len), os.path.join(odir, "contig_table")], ["\n\nBUILDING GRAPH...", "\n\nCalculating contig (internal edge) coverages..."]),
        5 : (["python3", "pair_reads.py", os.path.join(odir, "filter_R1table"),
                         os.path.join(odir, "filter_R2table"), os.path.join(odir, "filter_singleton_table"),
                         os.path.join(odir, "filter_paired_table")], ["\n\nPairing reads..."]),
        6 : (["python", "read_stats_core.py", os.path.join(odir, "filter_paired_table"),
                         str(args.k), str(args.read_len), os.path.join(odir, "read_stats.txt")], ["Calculating read statistics..."]),
        7 : (["python3", "find_new_variable_edges.py", os.path.join(odir, "filter_paired_table"),
                         os.path.join(odir, "adjacency_list"), os.path.join(odir, "contig_table"),
                         os.path.join(odir, "read_stats.txt"), os.path.join(odir, "edge_summary")], ["\n\nFinding all external edges..."]),
        8 : (["python3", "dominant_cycles.py", os.path.join(odir, "edge_summary"),
                         os.path.join(odir, "read_stats.txt"), str(args.alpha),
                         os.path.join(odir, "cycle_contig_table")], ["\n\nFinding cycles in the graph..."]),
        9 : (["python3", "cycle_coverages_2.py", args.sample, os.path.join(odir, "cycle_contig_table"),
                         os.path.join(odir, "filter_paired_table"), os.path.join(odir, "filter_singleton_table"),
                         os.path.join(odir, "contig_table"), os.path.join(odir, "read_stats.txt"),
                         os.path.join(odir, "cycle_covs_long"), os.path.join(odir, "cycle_cov_summary"),
                         os.path.join(odir, "out_contigs")], ["\n\nCalculating cycle coverage statistics in the cycle space..."]),
        10 : (["python3", "cycle_fastas.py", os.path.join(odir, "cycle_contig_table"),
                         os.path.join(odir, "renamed_final_contigs.fa"), os.path.join(odir, "read_stats.txt"),
                         os.path.join(odir, "cycles.fasta")], ["\n\nCreating cycle fastas..."]),
        11 : (["python3", "extract_p.py", os.path.join(odir, "cycle_cov_summary"),
                         os.path.join(odir, "cycles.fasta"), os.path.join(odir, "cycle_contig_table"),
                         str(args.maxpval), "0", "0", str(args.minscore), os.path.join(odir, "cycle_stats"), os.path.join(odir, "dominant_cycles_pre", "cycle_stats"),
                         os.path.join(odir, "dominant_cycles_pre", "cycles.fasta"), os.path.join(odir, "dominant_cycles_pre", "cycle_contig_table")], ["\n\nIdentifying cycles with low out reads...\n"]),
        12 : (["python3", "remove_high_singletons.py", os.path.join(odir, "cycle_covs_long"),
                         os.path.join(odir, "cycle_stats"), os.path.join(odir, "cycle_contig_table"),
                         os.path.join(odir, "cycle_cov_summary"), os.path.join(odir, "cycles.fasta"), str(args.minscore), str(args.maxpval),
                         os.path.join(odir, "dominant_cycles")], ["\n\nIdentifying dominant cycles...\n"]),
        13 : (["python3", "clean.py", os.path.join(odir)], [len(open(os.path.join(odir, "dominant_cycles", "dominant_cycle_stats")).readlines()) -1, "dominant cycles found!"])
    }

    return cmds

def parse_steps(steps: str):
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
    run_steps()
    print("\n\nFinished running DomCycle!")
    print("Stop time:", datetime.datetime.now())