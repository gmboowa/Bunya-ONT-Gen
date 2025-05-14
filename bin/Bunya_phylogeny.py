#!/usr/bin/env python3
"""
Robust Phylogenetic Analysis Pipeline for Bunyavirus with CLI arguments
"""

import os
import argparse
import tempfile
import subprocess
import multiprocessing
import shutil
from Bio import SeqIO
from datetime import datetime

def parse_arguments():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(description='Bunya Phylogenetic Analysis Pipeline')
    parser.add_argument('-r', '--reference', required=True,
                        help='Path to reference genome (GenBank or FASTA)')
    parser.add_argument('-i', '--input', required=True,
                        help='Path to text file listing sample assemblies')
    parser.add_argument('-o', '--output', default=None,
                        help='Output directory (default: timestamped)')
    parser.add_argument('-t', '--threads', type=int, default=8,
                        help='Number of CPU threads (default: 8)')
    parser.add_argument('-b', '--bootstrap', type=int, default=1000,
                        help='Bootstrap replicates (default: 1000)')
    return parser.parse_args()

def setup_environment(args):
    """Create output directory and process reference genome"""
    output_dir = args.output or f"phylogeny_{datetime.now().strftime('%Y%m%d')}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert reference to FASTA if GenBank
    ref_path = args.reference
    if ref_path.lower().endswith(('.gb', '.gbk')):
        fasta_ref = os.path.join(output_dir, 'reference.fasta')
        with open(fasta_ref, 'w') as f_out:
            record = SeqIO.read(ref_path, 'genbank')
            SeqIO.write(record, f_out, 'fasta')
        return output_dir, fasta_ref
    return output_dir, ref_path

def validate_samples(sample_file):
    """Validate input samples and return list of paths"""
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
        if seq_len < 11000:
            print(f"WARNING: {path} is unusually short ({seq_len} bp)")
        
        valid_samples.append(path)
    
    return valid_samples

def run_analysis(args, samples, output_dir, reference):
    """Execute the analysis pipeline"""
    # 1. Multiple Sequence Alignment
    print("\n=== Running MAFFT Alignment ===")
    
    # Combine all sequences into single temporary FASTA
    all_inputs = [reference] + samples
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as combined_fasta:
        for fasta_file in all_inputs:
            for record in SeqIO.parse(fasta_file, "fasta"):
                SeqIO.write(record, combined_fasta, "fasta")
        tmp_path = combined_fasta.name

    mafft_cmd = [
        "mafft",
        "--thread", str(args.threads),
        "--auto",
        "--reorder",
        "--adjustdirection",
        "--anysymbol",
        "--namelength", "1000",
        tmp_path
    ]

    try:
        with open(f"{output_dir}/alignment.fasta", "w") as outfile:
            result = subprocess.run(
                mafft_cmd,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
            if result.stderr:
                print("MAFFT Alignment Warnings:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"\nMAFFT Alignment Failed!\nError: {e.stderr}")
        raise SystemExit(1)
    finally:
        os.remove(tmp_path)

    # 2. Alignment Trimming
    print("\n=== Trimming Alignment ===")
    trimal_cmd = [
        "trimal",
        "-in", f"{output_dir}/alignment.fasta",
        "-out", f"{output_dir}/trimmed_alignment.fasta",
        "-gt", "0.9",
        "-cons", "60"
    ]
    subprocess.run(trimal_cmd, check=True)

    # 3. Phylogenetic Inference with Smart Thread Handling
    print("\n=== Building Phylogenetic Tree (with ModelFinder) ===")
    
    # Clean up previous runs
    for f in os.listdir(output_dir):
        if f.startswith("asfv_tree"):
            os.remove(os.path.join(output_dir, f))
    
    # Build IQ-TREE command
    iqtree_cmd = [
        "iqtree",
        "-s", f"{output_dir}/trimmed_alignment.fasta",
        "-m", "MFP",
        "-bb", str(args.bootstrap),
        "-alrt", "1000",
        "-nt", str(min(args.threads, multiprocessing.cpu_count())),
        "--redo",
        "-czb",
        "-seed", "12345",
        "-pre", f"{output_dir}/Bunya_tree"
    ]

    try:
        subprocess.run(iqtree_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\nIQ-TREE Failed! Error: {e.stderr}")
        print("Attempting automatic thread configuration...")
        iqtree_cmd[9] = "AUTO"
        subprocess.run(iqtree_cmd, check=True)

    # 4. Visualization
    print("\n=== Generating Visualization ===")
    r_script = f"""
    library(ggplot2)
    library(ggtree)
    tree <- read.tree("{output_dir}/Bunya_tree.treefile")
    p <- ggtree(tree) + 
      geom_tiplab(size=3) +
      geom_nodelab(aes(label=ifelse(as.numeric(label) > 70, label, '')), 
                   size=3, color="red") +
      theme_tree2() + 
      ggtitle("Bunya Phylogeny (IQ-TREE, {args.bootstrap} Bootstraps)")
    ggsave("{output_dir}/Bunya_phylogeny.pdf", plot=p, width=12, height=8)
    """
    with open(f"{output_dir}/visualize_tree.R", "w") as f:
        f.write(r_script)
    
    subprocess.run(["Rscript", f"{output_dir}/visualize_tree.R"], check=True)

def main():
    args = parse_arguments()
    
    # Setup environment and validate inputs
    output_dir, reference = setup_environment(args)
    samples = validate_samples(args.input)
    
    print(f"\nStarting ASFV Phylogenetic Analysis with:")
    print(f"- Reference genome: {args.reference}")
    print(f"- Samples: {len(samples)} assemblies")
    print(f"- Threads: {args.threads}")
    print(f"- Bootstrap replicates: {args.bootstrap}")
    print(f"- Output directory: {output_dir}")
    
    # Execute analysis pipeline
    run_analysis(args, samples, output_dir, reference)
    
    print("\nAnalysis successfully completed!")
    print(f"Results available in: {output_dir}")

if __name__ == "__main__":
    main()