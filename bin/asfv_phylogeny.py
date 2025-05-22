#!/usr/bin/env python3
"""
ASFV Phylogenetic Pipeline with ETE3 Visualization
"""

import os
import argparse
import tempfile
import subprocess
import multiprocessing
from Bio import SeqIO
from datetime import datetime
from ete3 import Tree, TreeStyle, NodeStyle, TextFace

def parse_arguments():
    parser = argparse.ArgumentParser(description='ASFV Phylogenetic Analysis Pipeline')
    parser.add_argument('-r', '--reference', required=True, help='Path to reference genome (GenBank or FASTA)')
    parser.add_argument('-i', '--input', required=True, help='Path to text file listing sample assemblies')
    parser.add_argument('-o', '--output', default=None, help='Output directory (default: timestamped)')
    parser.add_argument('-t', '--threads', type=int, default=8, help='Number of CPU threads (default: 8)')
    parser.add_argument('-b', '--bootstrap', type=int, default=1000, help='Bootstrap replicates (default: 1000)')
    return parser.parse_args()

def setup_environment(args):
    output_dir = args.output or f"results_phylo_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    os.makedirs(output_dir, exist_ok=True)
    ref_path = args.reference
    if ref_path.lower().endswith(('.gb', '.gbk')):
        fasta_ref = os.path.join(output_dir, 'reference.fasta')
        with open(fasta_ref, 'w') as f_out:
            record = SeqIO.read(ref_path, 'genbank')
            record.id = record.name = record.description = record.id.replace("_R_", "")
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
        # Clean record ID
        for r in records:
            r.id = r.name = r.description = r.id.replace("_R_", "")
        valid_samples.append(path)
    return valid_samples

def clean_iqtree_newick(infile, outfile):
    with open(infile) as fin, open(outfile, "w") as fout:
        tree_str = fin.read().strip()
        cleaned = tree_str.replace("'", "")
        fout.write(cleaned + "\n")

def visualize_tree_with_ete3(treefile, output_dir):
    try:
        tree = Tree(treefile, format=1)
    except Exception as e:
        print(f"❌ ETE3 tree parsing error: {e}")
        return

    for node in tree.traverse():
        nstyle = NodeStyle()
        nstyle["fgcolor"] = "black"
        nstyle["size"] = 5
        nstyle["vt_line_color"] = "#555555"
        nstyle["hz_line_color"] = "#555555"
        nstyle["vt_line_width"] = 2
        nstyle["hz_line_width"] = 2
        nstyle["draw_descendants"] = True

        if not node.is_leaf():
            try:
                support = float(node.name.split("/")[0]) if "/" in node.name else float(node.name)
                node.name = f"{support:.1f}"
                if support >= 90:
                    nstyle["fgcolor"] = "green"
                elif support >= 70:
                    nstyle["fgcolor"] = "orange"
                else:
                    nstyle["fgcolor"] = "red"
                nstyle["size"] = 6
            except ValueError:
                pass

        node.set_style(nstyle)

        if node.is_leaf():
            label = node.name.replace("_R_", "")  # Final cleanup safeguard
            name_face = TextFace(label, fsize=12, fgcolor="black")
            name_face.margin_right = 10
            node.add_face(name_face, column=0, position="branch-right")

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.show_branch_support = True
    ts.branch_vertical_margin = 30
    ts.min_leaf_separation = 20
    ts.scale = 150
    ts.show_scale = True
    ts.title.add_face(TextFace("ASFV Phylogenetic Tree (ETE3)", fsize=16, bold=True), column=0)

    tree.render(f"{output_dir}/asfv_ete3_tree.png", w=2200, dpi=300, tree_style=ts)
    print(f"✅ ETE3 tree saved to {output_dir}/asfv_ete3_tree.png")

def run_analysis(args, samples, output_dir, reference):
    print("\n=== Running MAFFT Alignment ===")
    all_inputs = [reference] + samples
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as combined_fasta:
        seen_ids = set()
        for fasta_file in all_inputs:
            for record in SeqIO.parse(fasta_file, "fasta"):
                clean_id = record.id.replace("_R_", "").replace(" ", "_")
                while clean_id in seen_ids:
                    clean_id += "_dup"
                record.id = record.name = record.description = clean_id
                seen_ids.add(clean_id)
                SeqIO.write(record, combined_fasta, "fasta")
        tmp_path = combined_fasta.name

    subprocess.run([
        "mafft", "--thread", str(args.threads), "--auto", "--reorder",
        "--adjustdirection", "--anysymbol", "--namelength", "1000", tmp_path
    ], stdout=open(f"{output_dir}/alignment.fasta", "w"), check=True)
    os.remove(tmp_path)

    print("\n=== Trimming Alignment ===")
    subprocess.run([
        "trimal", "-in", f"{output_dir}/alignment.fasta",
        "-out", f"{output_dir}/trimmed_alignment.fasta", "-gt", "0.9", "-cons", "60"
    ], check=True)

    print("\n=== Running IQ-TREE ===")
    tree_prefix = f"{output_dir}/asfv_tree"
    subprocess.run([
        "iqtree", "-s", f"{output_dir}/trimmed_alignment.fasta", "-m", "MFP",
        "-bb", str(args.bootstrap), "-nt", str(args.threads),
        "--redo", "-czb", "-seed", "12345", "-pre", tree_prefix
    ], check=True)

    print("\n=== Generating ETE3 Visualization ===")
    clean_iqtree_newick(f"{tree_prefix}.treefile", f"{output_dir}/cleaned_tree.nwk")
    visualize_tree_with_ete3(f"{output_dir}/cleaned_tree.nwk", output_dir)

    print("\n=== Writing Summary ===")
    alignment_file = f"{output_dir}/trimmed_alignment.fasta"
    summary_file = f"{output_dir}/summary.txt"
    iqtree_log = f"{tree_prefix}.iqtree"
    model = "Unknown"
    if os.path.exists(iqtree_log):
        with open(iqtree_log) as f:
            for line in f:
                if "Best-fit model:" in line:
                    model = line.split(":")[-1].strip()
                    break

    num_seqs, aln_len = 0, 0
    with open(alignment_file) as f:
        for line in f:
            if line.startswith(">"):
                num_seqs += 1
            else:
                aln_len += len(line.strip())

    with open(summary_file, "w") as f:
        f.write("ASFV Phylogenetic Pipeline Summary\n")
        f.write("==================================\n")
        f.write(f"Alignment file: {alignment_file}\n")
        f.write(f"Sequences: {num_seqs}\n")
        f.write(f"Alignment length: {aln_len} bp\n")
        f.write(f"IQ-TREE Model: {model}\n")
        f.write("Command: MAFFT, TRIMAL, IQ-TREE, ETE3\n")

    print(f"✅ Summary saved to {summary_file}")

def main():
    args = parse_arguments()
    output_dir, reference = setup_environment(args)
    samples = validate_samples(args.input)
    run_analysis(args, samples, output_dir, reference)

if __name__ == "__main__":
    main()

