# ASFV-ONT-Gen: A Pipeline for African Swine Fever Virus ONT Whole-Genome Analysis

**ASFV-ONT-Gen** is a robust and modular pipeline designed for analyzing whole-genome sequencing data of African Swine Fever Virus (ASFV) using Oxford Nanopore Technologies (ONT). It covers preprocessing, assembly, annotation, and phylogenetic analysis using a reference genome such as `NC_044959.2.fasta`.

---

# ASFV Oxford Nanopore Pipeline

A bioinformatics pipeline for processing Oxford Nanopore sequencing data, performing quality control, assembly, variant calling, annotation and taxonomic classification.

## Features

- **Quality Control**: FastQC and NanoPlot for read quality assessment.
- **Read Mapping**: Minimap2 for alignment to a reference genome.
- **Taxonomic Classification**: Kraken2 with a custom viral database.
- **Assembly**: Both de novo (Flye) and reference-based (Medaka) assembly.
- **Variant Calling**: Bcftools for SNP/indel detection and snpEff for annotation.
- **Multi-threading**: Utilizes concurrent processing for efficiency.
- **Summary Reports**: MultiQC, mapped reads summary, and assembly statistics.

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installed.
- Java 21+ (required for snpEff).

### Setup

1. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/AFSV_ONT_Pipeline.git
   cd AFSV_ONT_Pipeline


![Pipeline Workflow](ASFV-ONT-Gen_Workflow.png)

## Installation <a name="installation"></a>

Create and activate the Conda environment:

```bash
# Clone repository
git clone https://github.com/gmboowa/ASFV-ONT-Gen.git
cd ASFV-ONT-Gen

# Create conda environment
conda env create -f AFSV_ont.yml
conda activate AFSV_ont
conda install -c bioconda -c conda-forge fastqc nanoplot minimap2 samtools bcftools medaka multiqc spades kraken2 mafft fasttree seqtk flye krona snpeff -y

# Install Kraken2 database (optional)
scripts/setup_kraken_db.sh
```

## Usage <a name="usage"></a>

```bash

Command

python AFSV_ont_pipeline.py -inputs samples.txt -reference reference.fasta -threads 8

```

## Inputs and Outputs <a name="inputs-and-outputs"></a>

### Inputs

| File Type      | Description            |
|----------------|------------------------|
| *.fastq        | ONT sequence files     |
| reference.fasta| ASFV reference genome  |

### Output Structure

Results are saved in ./results:

```bash
results/
├── qc/                  # FastQC and NanoPlot reports
├── nanoplot/            # NanoPlot visualizations
├── multiqc/             # MultiQC summary
├── <summary>/
|   |-species_abundance_summary.csv    # Kraken2 species abundance summaries
├── krona/               # Krona interactive reports
├── <sample_name>/       # Per-sample outputs
│   ├── denovo_assembly/ # Flye assembly files
│   ├── reference_assembly/ # Medaka consensus
│   ├── variants/        # VCF and annotated variants --- SNPS and Indels
│   └── kraken2_report.txt
├── mapped_reads_summary.tsv       # Mapped reads statistics
└── assembly_stats_summary.tsv     # Assembly metrics (N50, L50, etc.)
```

## Dependencies <a name="dependencies"></a>

Workflow Overview

Quality Control: Generate QC reports for raw reads.

Read Mapping: Align reads to the reference genome.

Taxonomic Classification: Identify species composition with Kraken2.

Assembly:

De novo assembly using Flye.

Reference-based consensus with Medaka.

Variant Calling: Identify variants and annotate using snpEff.

Summarization: Aggregate results into TSV files and MultiQC reports.

Notes
The Kraken2 database (kraken2_viral_db) will be auto-built if missing.

Ensure SNPEFF_HOME is set (automatically configured in the script).

Edit Entrez.email in the script if NCBI API access is required.

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
