#!/usr/bin/env python3
"""
Robust Phylogenetic Analysis Pipeline for ASFV with CLI arguments
"""

import os
import argparse
import tempfile
import subprocess
import multiprocessing
from Bio import SeqIO
from datetime import datetime

def parse_arguments():
    parser = argparse.ArgumentParser(description='ASFV Phylogenetic Analysis Pipeline')
    parser.add_argument('-r', '--reference', required=True, help='Path to reference genome (GenBank or FASTA)')
    parser.add_argument('-i', '--input', required=True, help='Path to text file listing sample assemblies')
    parser.add_argument('-o', '--output', default=None, help='Output directory (default: timestamped)')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of CPU threads (default: 8)')
    parser.add_argument('-b', '--bootstrap', type=int, default=1000, help='Bootstrap replicates (default: 1000)')
    return parser.parse_args()

def setup_environment(args):
    # Update the output directory format to include date and time
    output_dir = args.output or f"results_phylo_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    os.makedirs(output_dir, exist_ok=True)
    
    ref_path = args.reference
    if ref_path.lower().endswith(('.gb', '.gbk')):
        fasta_ref = os.path.join(output_dir, 'reference.fasta')
        with open(fasta_ref, 'w') as f_out:
            record = SeqIO.read(ref_path, 'genbank')
            SeqIO.write(record, f_out, 'fasta')
        return output_dir, fasta_ref
    return output_dir, ref_path

def validate_samples(sample_file):
    with open(sample_file) as f:
        samples = [line.strip() for line in f if line.strip()]
    valid_samples = []
    for path in samples:
        if not os.path.exists(path):
            raise SystemExit(f"ERROR: Sample file {path} not found")
        records = list(SeqIO.parse(path, 'fasta'))
        if len(records) != 1:
            raise SystemExit(f"ERROR: {path} contains multiple sequences")
        seq_len = len(records[0].seq)
        if seq_len < 170000:
            print(f"WARNING: {path} is unusually short ({seq_len} bp)")
        valid_samples.append(path)
    return valid_samples

def run_analysis(args, samples, output_dir, reference):
    print("\n=== Running MAFFT Alignment ===")
    all_inputs = [reference] + samples
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as combined_fasta:
        seen_ids = set()
        for fasta_file in all_inputs:
            for record in SeqIO.parse(fasta_file, "fasta"):
                orig_id = record.id
                new_id = orig_id
                counter = 1
                while new_id in seen_ids:
                    new_id = f"{orig_id}_{counter}"
                    counter += 1
                record.id = record.name = record.description = new_id
                seen_ids.add(new_id)
                SeqIO.write(record, combined_fasta, "fasta")
        tmp_path = combined_fasta.name

    mafft_cmd = [
        "mafft", "--thread", str(args.threads), "--auto", "--reorder", "--adjustdirection",
        "--anysymbol", "--namelength", "1000", tmp_path
    ]
    try:
        with open(f"{output_dir}/alignment.fasta", "w") as outfile:
            result = subprocess.run(mafft_cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
            if result.stderr:
                print("MAFFT Alignment Warnings:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"\nMAFFT Alignment Failed!\nError: {e.stderr}")
        raise SystemExit(1)
    finally:
        os.remove(tmp_path)

    print("\n=== Trimming Alignment ===")
    trimal_cmd = [
        "trimal", "-in", f"{output_dir}/alignment.fasta", "-out", f"{output_dir}/trimmed_alignment.fasta",
        "-gt", "0.9", "-cons", "60"
    ]
    subprocess.run(trimal_cmd, check=True)

    print("\n=== Building Phylogenetic Tree (with ModelFinder) ===")
    for f in os.listdir(output_dir):
        if f.startswith("asfv_tree"):
            os.remove(os.path.join(output_dir, f))

    iqtree_cmd = [
        "iqtree", "-s", f"{output_dir}/trimmed_alignment.fasta", "-m", "MFP",
        "-bb", str(args.bootstrap), "-alrt", "1000", "-nt", str(min(args.threads, multiprocessing.cpu_count())),
        "--redo", "-czb", "-seed", "12345", "-pre", f"{output_dir}/asfv_tree"
    ]
    try:
        subprocess.run(iqtree_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\nIQ-TREE Failed! Error: {e.stderr}")
        print("Attempting automatic thread configuration...")
        if "-nt" in iqtree_cmd:
            nt_index = iqtree_cmd.index("-nt")
            iqtree_cmd[nt_index + 1] = "AUTO"
        subprocess.run(iqtree_cmd, check=True)

    print("\n=== Generating Visualization ===")
    r_script = f"""    library(ggplot2)
    library(ggtree)
    library(ape)

    tree <- read.tree("{output_dir}/asfv_tree.treefile")
    p <- ggtree(tree, layout = "rectangular") +
        geom_tiplab(size = 3) +
        geom_nodelab(aes(label=label, color=as.numeric(label)), size=3) +
        scale_color_gradient2(low="gray", mid="blue", high="red", midpoint=70, name="Bootstrap") +
        theme_tree2() +
        ggtitle("ASFV Phylogeny (IQ-TREE, Bootstrap Support)")

    ggsave("{output_dir}/asfv_phylogeny.pdf", plot = p, width = 12, height = 8)
    ggsave("{output_dir}/asfv_phylogeny.svg", plot = p, width = 12, height = 8)
    """
    with open(f"{output_dir}/visualize_tree.R", "w") as f:
        f.write(r_script)
    subprocess.run(["Rscript", f"{output_dir}/visualize_tree.R"], check=True)

    print("\n=== Writing Summary Report ===")
    alignment_file = f"{output_dir}/trimmed_alignment.fasta"
    summary_file = f"{output_dir}/summary.txt"
    iqtree_log = f"{output_dir}/asfv_tree.iqtree"

    num_seqs, aln_len = 0, 0
    with open(alignment_file) as f:
        for line in f:
            if line.startswith(">"):
                num_seqs += 1
            else:
                aln_len += len(line.strip())

    model = "Unknown"
    if os.path.exists(iqtree_log):
        with open(iqtree_log) as f:
            for line in f:
                if "Best-fit model:" in line:
                    model = line.strip().split(":")[-1].strip()
                    break

    with open(summary_file, "w") as f:
        f.write("ASFV Phylogenetic Pipeline Summary\n")
        f.write("==================================\n\n")
        f.write(f"Alignment file: {alignment_file}\n")
        f.write(f"Number of sequences: {num_seqs}\n")
        f.write(f"Alignment length: {aln_len} bp\n")
        f.write(f"IQ-TREE best-fit model: {model}\n\n")
        f.write("Commands used:\n")
        f.write("MAFFT:\n")
        f.write("  mafft --thread {args.threads} --auto --reorder --adjustdirection \
--anysymbol --namelength 1000 input.fasta > alignment.fasta\n")
        f.write("IQ-TREE:\n")
        f.write(f"  {' '.join(iqtree_cmd)}\n")

    print(f"✅ Summary written to {summary_file}")

def main():
    args = parse_arguments()
    output_dir, reference = setup_environment(args)
    samples = validate_samples(args.input)
    print(f"\nStarting ASFV Phylogenetic Analysis with:")
    print(f"- Reference genome: {args.reference}")
    print(f"- Samples: {len(samples)} assemblies")
    print(f"- Threads: {args.threads}")
    print(f"- Bootstrap replicates: {args.bootstrap}")
    print(f"- Output directory: {output_dir}")
    run_analysis(args, samples, output_dir, reference)
    print("\nAnalysis successfully completed!")
    print(f"Results available in: {output_dir}")

if __name__ == "__main__":
    main()
