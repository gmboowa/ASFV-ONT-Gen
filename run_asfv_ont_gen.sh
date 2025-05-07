#!/bin/bash

# Example pipeline runner
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i) INPUT="$2"; shift ;;
        -r) REF="$2"; shift ;;
        -o) OUT="$2"; shift ;;
        --min_length) MINLEN="$2"; shift ;;
        --threads) THREADS="$2"; shift ;;
    esac
    shift
done

mkdir -p $OUT

echo "Running ASFV-ONT-Gen with:"
echo "Input: $INPUT"
echo "Reference: $REF"
echo "Output: $OUT"

# Placeholder steps
fastp -i $INPUT -o $OUT/filtered.fastq
NanoStat --fastq $OUT/filtered.fastq -o $OUT/nanostat.txt
NanoPlot --fastq $OUT/filtered.fastq -o $OUT/nanoplot/
filtlong --min_length $MINLEN $OUT/filtered.fastq > $OUT/longreads.fastq
flye --nano-raw $OUT/longreads.fastq --out-dir $OUT/flye --genome-size 200k --threads $THREADS
medaka_consensus -i $OUT/longreads.fastq -d $REF -o $OUT/medaka -t $THREADS -m r941_min_sup_g507
racon $OUT/longreads.fastq $OUT/medaka/consensus.fasta $REF > $OUT/racon_polished.fasta
prokka $OUT/racon_polished.fasta --outdir $OUT/annotation --prefix asfv
mafft $OUT/racon_polished.fasta > $OUT/aligned.fasta
iqtree -s $OUT/aligned.fasta -nt AUTO -bb 1000 -m GTR+G
