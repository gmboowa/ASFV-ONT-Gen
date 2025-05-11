# ASFV-ONT-Gen: A Pipeline for African Swine Fever Virus ONT Whole-Genome Analysis

**ASFV-ONT-Gen** is a modular pipeline designed for analyzing whole-genome sequencing data of African Swine Fever Virus (ASFV) using Oxford Nanopore Technologies (ONT). It covers preprocessing, assembly, annotation, and phylogenetic analysis using a selected reference genome such as `NC_044959.2.fasta`.

---

## Features

- **Quality control**: FastQC + NanoPlot + MultiQC for read quality assessment.
- **Read mapping**: Minimap2 for alignment to a selected reference genome.
- **Taxonomic classification**: Kraken2 with viral database.
- **Assembly**: Both de novo (Flye + Medaka polish) & reference-based assembly.
- **Variant calling & annotation**: Bcftools for SNP/INDEL detection & snpEff for annotation.
- **Multi-threading**: Utilizes concurrent processing for efficiency.
- **Summary reports**: MultiQC, mapped reads summary, & assembly statistics.
- **Phylogenetics**:
  - MAFFT for multiple sequence alignment
  - IQ-TREE for robust phylogeny (ModelFinder + UFBoot)
  - ggtree for tree visualization

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installed.
- Java 21+ (required for snpEff).
  
### Install remaining dependencies
   conda install -c bioconda -c conda-forge fastqc nanoplot minimap2 samtools bcftools \
   medaka multiqc spades kraken2 mafft fasttree seqtk flye krona snpeff iqtree trimal r-ggplot2 r-ggtree

### Pipeline Workflow


![Pipeline Workflow](ASFV-ONT-Gen_Workflow.png)

## Installation <a name="installation"></a>

Make sure you have Miniconda or Anaconda on your system

Create and activate the Conda environment:

```

# Create conda environment for assembly, variant calling and taxonomic profiling

conda env create -f AFSV_ont.yml
conda activate AFSV_ont
conda install -c bioconda -c conda-forge fastqc nanoplot minimap2 samtools bcftools medaka
multiqc spades kraken2 mafft fasttree seqtk flye krona snpeff -y

# Install Kraken2 database 
scripts/setup_kraken_db.sh

## Usage

# Run whole-genome analysis

**`python AFSV_ont_pipeline.py -inputs fastq_sample.txt -reference reference.fasta -threads 8`**



# Create conda environment for phylogenetic analysis

conda env create -f asfv_phylogeny.yml

conda activate asfv_phylogeny

## Usage 

# Run phylogenetic inference

**`asfv_phylogeny.py  -r NC_044959.2.gb -i fasta_sample.txt -o results -t 8 -b 1000`**


```

## Inputs and Outputs <a name="inputs-and-outputs"></a>

### Inputs for genomic analysis

| File Type      | Description            |
|----------------|------------------------|
| *.fastq        | ONT sequence files     |
| reference.fasta| ASFV reference genome  |

### Inputs for phylogenetic analysis 

Configuration 


| Parameter         | Description                          | Default       |
|-------------------|--------------------------------------|---------------|
| `-inputs`         | File with FASTQ paths                | Required      |
| `-reference`      | ASFV reference genome                | Required      |
| `-threads`        | CPU threads for parallel steps       | `8`           |
| `-bootstrap`      | Phylogenetic bootstrap replicates    | `1000`        |
| `-min_coverage`   | Consensus calling threshold          | `20x`         |
| `-tree_model`     | IQ-TREE substitution model           | `MFP (auto)`  |
| `-aln_consensus`  | trimAl conservation threshold        | `60%`         |

### Output structure

Results are saved in ./results:

```
results/

├── 00_RawDataQC/          # FastQC/NanoPlot reports
├── 01_Assemblies/         # Flye & Medaka outputs
├── 02_Variants/           # VCF files & annotations
├── 03_Phylogeny/
│   ├── alignment.fasta    # MAFFT multiple alignment
│   ├── trimmed_alignment/ # trimAl filtered sequences
│   ├── iqtree_results/    # Tree files + support values
│   └── phylogeny.pdf      # Final ggtree visualization
├── 04_Taxonomy/           # Kraken2/Krona reports
└── reports/               # MultiQC + summary stats


```

## Troubleshooting

**Common Issues:**

- **Memory errors:** Reduce thread count with `-threads 4`
- **Assembly failures:** Check `results/qc/nanoplot/` for read quality
- **Dependency issues:** Update conda with `conda env update -f AFSV_ont.yml`

## License
MIT License 
