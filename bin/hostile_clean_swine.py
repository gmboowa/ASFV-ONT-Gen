#!/usr/bin/env python3
"""
hostile_clean_swine.py - Enhanced version with better file handling
"""

import argparse
import subprocess
import re
import zipfile
import shutil
import sys
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='ONT read cleaning pipeline with smart reference handling',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   
    parser.add_argument('--fasta', required=True,
                        help='NCBI accession (GCF_/GCA_) or path to FASTA')
    parser.add_argument('-i', '--input', required=True,
                        help='File list containing FASTQ paths')
    parser.add_argument('-o', '--output', required=True,
                        help='Main output directory')
    parser.add_argument('--index-dir', default='minimap2_index',
                        help='Directory for reference indexes')
    parser.add_argument('--reference-dir', default='references',
                        help='Directory for downloaded references')
   
    return parser.parse_args()

def validate_dependencies():
    """Check for required executables"""
    for cmd in ['datasets', 'minimap2', 'hostile']:
        if not shutil.which(cmd):
            raise SystemExit(f"Missing dependency: {cmd}. Install with 'conda install -c bioconda ncbi-datasets-cli minimap2 hostile'")

def extract_accession(input_str):
    """Extract valid NCBI accession from input string"""
    match = re.search(r'(GC[AF]_\d+)', input_str)
    if not match:
        raise SystemExit(f"Invalid accession format: {input_str}")
    return match.group(1).split('.')[0]

def find_existing_reference(accession, reference_dir):
    """Check for existing reference files"""
    ref_dir = Path(reference_dir)
    patterns = [
        f"{accession}*_genomic.fna",
        f"{accession}*.fna",
        f"{accession}*.fa",
        f"{accession}*.fasta"
    ]
   
    for pattern in patterns:
        matches = list(ref_dir.rglob(pattern))
        if matches:
            return matches[0]
    return None

def get_reference(fasta_input, reference_dir):
    """Get reference path with smart checking"""
    local_path = Path(fasta_input)
    if local_path.exists():
        print(f"Using local reference: {local_path.resolve()}")
        return local_path

    base_accession = extract_accession(fasta_input)
    existing_ref = find_existing_reference(base_accession, reference_dir)
    if existing_ref:
        print(f"Using existing reference: {existing_ref.resolve()}")
        return existing_ref

    print(f"Downloading reference: {base_accession}")
    ref_dir = Path(reference_dir)
    ref_dir.mkdir(parents=True, exist_ok=True)
    zip_path = ref_dir / f"{base_accession}.zip"
   
    try:
        subprocess.run(
            ['datasets', 'download', 'genome', 'accession', base_accession,
             '--include', 'genome', '--filename', str(zip_path)],
            check=True
        )
       
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(ref_dir)
       
        new_ref = find_existing_reference(base_accession, ref_dir)
        if not new_ref:
            raise SystemExit(f"Download failed for {base_accession}")
        print(f"Successfully downloaded reference: {new_ref.resolve()}")
        return new_ref
   
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"Failed to download {base_accession}: {e}")
    finally:
        if zip_path.exists():
            zip_path.unlink()

def build_index(fasta_path, index_dir):
    """Build minimap2 index if missing"""
    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)
    index_file = index_dir / "reference.mmi"
   
    if not index_file.exists():
        print(f"Building minimap2 index at: {index_file.resolve()}")
        subprocess.run(
            ['minimap2', '-x', 'map-ont', '-d', str(index_file), str(fasta_path)],
            check=True
        )
    else:
        print(f"Using existing index: {index_file.resolve()}")
    return index_file

def validate_fastq_file(fastq_path):
    """Validate FASTQ file exists and is accessible"""
    path = Path(fastq_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTQ file not found: {path.resolve()}")
    if not path.is_file():
        raise ValueError(f"Path is not a file: {path.resolve()}")
    if path.stat().st_size == 0:
        raise ValueError(f"Empty FASTQ file: {path.resolve()}")
    return path.resolve()

def process_samples(input_list, index_file, output_root):
    """Clean samples using hostile with enhanced error handling"""
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    processed = 0
    failed_samples = []
   
    with open(input_list) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            try:
                fastq_path = validate_fastq_file(line)
                sample_id = fastq_path.stem.split('.')[0]
                sample_out = output_root / f"{sample_id}_cleaned"
                sample_out.mkdir(parents=True, exist_ok=True)
               
                print(f"\nProcessing sample: {sample_id}")
                print(f"Input FASTQ: {fastq_path}")
                print(f"Output directory: {sample_out.resolve()}")
               
                cmd = [
                    'hostile', 'clean',
                    '--fastq1', str(fastq_path),
                    '--index', str(index_file),
                    '--aligner', 'minimap2',
                    '-o', str(sample_out)
                ]
               
                subprocess.run(cmd, check=True)
                processed += 1
                print(f"Successfully processed: {sample_id}")
               
            except Exception as e:
                print(f"❌ Failed to process {line}: {str(e)}")
                failed_samples.append((line, str(e)))
                continue
   
    if failed_samples:
        print("\nFailed samples:")
        for sample, error in failed_samples:
            print(f"- {sample}: {error}")
   
    return processed

def main():
    args = parse_arguments()
    validate_dependencies()
   
    try:
        # Handle reference
        ref_path = get_reference(args.fasta, args.reference_dir)
        index_file = build_index(ref_path, args.index_dir)
       
        # Process samples
        print(f"\n{' PROCESSING STARTED ':=^60}")
        processed = process_samples(args.input, index_file, args.output)
        print(f"\nCompleted {processed} samples successfully")
       
    except Exception as e:
        print(f"\n❌ Pipeline failed: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
