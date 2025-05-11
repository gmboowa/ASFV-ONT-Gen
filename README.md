
# ASFV-ONT-Gen: A Pipeline for African Swine Fever Virus ONT Whole-Genome Analysis

**ASFV-ONT-Gen** is a modular pipeline for analyzing whole-genome sequencing data of African Swine Fever Virus (ASFV) using Oxford Nanopore Technologies (ONT). It covers preprocessing, assembly, annotation, and phylogenetic analysis using a reference genome such as `NC_044959.2.fasta`.

---

## 🚀 Features

- **Quality control**: FastQC, NanoPlot, MultiQC
- **Read mapping**: Minimap2
- **Taxonomic classification**: Kraken2 with viral database
- **Assembly**: De novo (Flye + Medaka polish) & reference-based
- **Variant calling & annotation**: BCFtools and SnpEff
- **Multi-threading**: Efficient parallel processing
- **Summary reports**: MultiQC, assembly stats, mapping summaries
- **Phylogenetics**: MAFFT, IQ-TREE, and ggtree

---
### Pipeline workflow


![Pipeline Workflow](ASFV-ONT-Gen_Workflow.png)

## 🛠️ Installation

### Prerequisites

- [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)

- Java 21+ (for SnpEff)

### Environment setup for genomic analysis

```bash
conda env create -f AFSV_ont.yml
conda activate AFSV_ont
conda install -c bioconda -c conda-forge fastqc nanoplot minimap2 samtools bcftools medaka \
multiqc spades kraken2 mafft fasttree seqtk flye krona snpeff -y
```

### Setup Kraken2 database

```bash
bash scripts/setup_kraken_db.sh
```

---

## 📦 Usage

### Run whole-genome analysis

```bash
python AFSV_ont_pipeline.py -inputs fastq_sample.txt -reference reference.fasta -threads 8
```

---

### Environment setup for phylogenetic analysis

```bash
conda env create -f asfv_phylogeny.yml
conda activate asfv_phylogeny
```

### Run phylogenetic inference

```bash
asfv_phylogeny.py -r NC_044959.2.gb -i fasta_sample.txt -o results -t 8 -b 1000
```

---

## Inputs and outputs

### Inputs for genomic analysis

| File Type        | Description                |
|------------------|----------------------------|
| `*.fastq`        | ONT sequence files         |
| `reference.fasta`| ASFV reference genome      |

### Inputs for Phylogenetic Analysis

| Parameter         | Description                          | Default       |
|-------------------|--------------------------------------|---------------|
| `-inputs`         | File with FASTQ paths                | Required      |
| `-reference`      | ASFV reference genome                | Required      |
| `-threads`        | CPU threads for parallel steps       | `8`           |
| `-bootstrap`      | Phylogenetic bootstrap replicates    | `1000`        |
| `-min_coverage`   | Consensus calling threshold          | `20x`         |
| `-tree_model`     | IQ-TREE substitution model           | `MFP (auto)`  |
| `-aln_consensus`  | trimAl conservation threshold        | `60%`         |

---

### Output Structure

```text
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

---

## Troubleshooting

| Issue               | Recommendation                                                  |
|---------------------|-----------------------------------------------------------------|
| Memory errors       | Reduce thread count (e.g., `-threads 4`)                        |
| Assembly failures   | Inspect quality reports in `results/qc/nanoplot/`              |
| Dependency problems | Update Conda with `conda env update -f AFSV_ont.yml`           |

---

## License

MIT License
