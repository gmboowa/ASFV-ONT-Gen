
# ASFV-ONT-Gen: A Pipeline for African Swine Fever Virus ONT Whole-Genome Analysis

**ASFV-ONT-Gen** is a robust and modular pipeline designed for analyzing whole-genome sequencing data of African Swine Fever Virus (ASFV) using Oxford Nanopore Technologies (ONT). It covers preprocessing, assembly, annotation, and phylogenetic analysis using a reference genome such as `NC_044959.2.fasta`.

---

## Features

- Preconfigured for ASFV long-read sequencing
- Quality control using Fastp and NanoPlot
- Assembly with Flye, Raven, or Canu
- Reference-based polishing with Medaka and Racon
- Functional annotation and gene prediction
- Phylogenetic analysis using MAFFT and IQ-TREE

---

## Workflow Overview

![ASFV-ONT-Gen Workflow](ASFV-ONT-Gen_Workflow.png)

> *A visual representation of the ASFV-ONT-Gen pipeline, including both de novo and reference-based assembly, variant calling, indel analysis, and phylogenetic inference.*

---

## Input

- Raw FASTQ files from ONT sequencing
- ASFV reference genome (e.g., `NC_044959.2.fasta`)

---

## Installation

```bash
conda env create -f installer.yml
conda activate asfv-ont-gen
```

---

## Usage

```bash
bash run_asfv_ont_gen.sh \
  -i data/reads.fastq \
  -r references/NC_044959.2.fasta \
  -o results/ \
  --min_length 1000 \
  --threads 8
```

Edit the `config.yaml` for tool-specific parameters.

---

## Output

- Filtered reads
- Assembled genome and polished sequences
- Functional annotations
- Multiple sequence alignment (`.aln`)
- Phylogenetic tree (`.nwk`)

---

## Citation

Please cite this repository and each software tool used within the workflow.




