# ASFV-ONT-Gen: A Pipeline for African Swine Fever Virus ONT Whole-Genome Analysis

**ASFV-ONT-Gen** is a robust and modular pipeline designed for analyzing whole-genome sequencing data of African Swine Fever Virus (ASFV) using Oxford Nanopore Technologies (ONT). It covers preprocessing, assembly, annotation, and phylogenetic analysis using a reference genome such as `NC_044959.2.fasta`.

---

## Table of Contents
1. [Overview](#overview)
2. [Workflow](#workflow)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Inputs/Outputs](#inputs-and-outputs)
6. [Dependencies](#dependencies)
7. [Examples](#examples)
8. [Troubleshooting](#troubleshooting)
9. [License](#license)


![Pipeline Workflow](ASFV-ONT-Gen_Workflow.png)

## Installation <a name="installation"></a>

```bash
# Clone repository
git clone https://github.com/gmboowa/ASFV-ONT-Gen.git
cd ASFV-ONT-Gen

# Create conda environment
conda env create -f AFSV_ont.yml
conda activate AFSV_ont

# Install Kraken2 database (optional)
scripts/setup_kraken_db.sh
```

## Usage <a name="usage"></a>

```bash

python AFSV_ont_pipeline.py -inputs samples.txt -reference reference.fasta -threads 8

```

## Inputs and Outputs <a name="inputs-and-outputs"></a>

### Inputs

| File Type      | Description            |
|----------------|------------------------|
| *.fastq        | ONT sequence files     |
| reference.fasta| ASFV reference genome  |

### Output Structure

```bash
results/
├── qc/                  # Quality reports
├── assemblies/De no vo & reference-based          # Final genomes
├── variants_called/            # SNP/indel calls
├── phylogeny/           # Tree files
└── reports/             # Summary reports
```

## Dependencies <a name="dependencies"></a>

| Tool       | Version | Role        |
|------------|---------|-------------|
| Flye       | 2.9+    | Assembly    |
| Medaka     | 1.7+    | Polishing   |
| Minimap2   | 2.24+   | Alignment   |
| FastQC     | 0.11+   | QC          |
| MultiQC    | 1.13+   | Reporting   |

## Examples <a name="examples"></a>

```bash
# Run test dataset
python AFSV_ont_pipeline.py \
  -inputs test_data/samples.txt \
  -reference test_data/NC_044959.2.fasta
```

## Troubleshooting <a name="troubleshooting"></a>

**Common Issues:**

- **Memory errors:** Reduce thread count with `-threads 4`
- **Assembly failures:** Check `results/qc/nanoplot/` for read quality
- **Dependency issues:** Update conda with `conda env update -f AFSV_ont.yml`

## License <a name="license"></a>
MIT License 
