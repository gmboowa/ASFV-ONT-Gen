import ssl
from Bio import Entrez
import logging
import subprocess
import shutil
import os
import re
import time
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd

# Fix SSL issue for NCBI access
ssl._create_default_https_context = ssl._create_unverified_context
Entrez.email = "gmboowa@gmail.com"

BIOCONDA_TOOLS = [
    "fastqc", "nanoplot", "minimap2", "samtools", "bcftools", "medaka",
    "multiqc", "spades", "kraken2", "mafft", "fasttree", "seqtk", "flye", "krona", "snpEff"
]

TOOL_MAPPING = {
    "fastqc": "fastqc",
    "nanoplot": "nanoplot",
    "minimap2": "minimap2",
    "samtools": "samtools",
    "bcftools": "bcftools",
    "medaka": "medaka_consensus",
    "multiqc": "multiqc",
    "spades": "spades.py",
    "kraken2": "kraken2",
    "mafft": "mafft",
    "fasttree": "FastTree",
    "seqtk": "seqtk",
    "flye": "flye",
    "krona": "ktImportText",
    "snpEff": "snpEff"
}

mapped_reads_counts = []
assembly_stats = []

def setup_directories(results_dir):
    SUMMARY_DIR = results_dir / "summary"
    KRONA_DIR = results_dir / "krona"
    for subdir in ["qc", "nanoplot", "summary", "krona", "multiqc"]:
        (results_dir / subdir).mkdir(parents=True, exist_ok=True)
    return SUMMARY_DIR, KRONA_DIR

def create_multiqc_config(results_dir):
    config_content = """\
table_columns_visible:
  fastqc_top_overrepresented_sequences_table: False

table_columns_hidden:
  fastqc_top_overrepresented_sequences_table: True

suppress_errors: true
"""
    config_path = results_dir / "multiqc" / "multiqc_config.yaml"
    with open(config_path, 'w') as f:
        f.write(config_content)
    return config_path

def quality_control(fastq_list, results_dir):
    for fq in fastq_list:
        run_cmd(f"fastqc {fq} -o {results_dir}/qc")
        run_cmd(
            f"nanoplot --fastq {fq} -o {results_dir}/nanoplot "
            f"--minlength 100 --downsample 10000 --plots dot "
            f"--title {Path(fq).stem} --loglength"
        )

def run_cmd(cmd):
    logging.info(f"Running command: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}\n{e}")
        
def ensure_java_version():
    print("\n🔍 Java Check:")
    try:
        result = subprocess.run(["java", "-version"],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                text=True)
        version_match = re.search(r'version "(\d+)\.', result.stdout)
        if version_match:
            java_version = int(version_match.group(1))
            if java_version >= 21:
                print(f"✅ Java {java_version} verified")
                return
            print(f"❌ Java {java_version} detected. Requires Java 21+")
            raise RuntimeError("Insufficient Java version")
        raise RuntimeError("Could not determine Java version")
    except FileNotFoundError:
        print("❌ Java not found in PATH")
        raise

def configure_snpeff(ref_fasta):
    genome_id = Path(ref_fasta).stem

    # Locate active conda environment
    conda_prefix = Path(os.environ.get("CONDA_PREFIX", ""))
    if not conda_prefix or not conda_prefix.exists():
        raise EnvironmentError("Conda environment not detected or CONDA_PREFIX not set.")

    # Search for snpEff.jar in conda environment share directory
    jar_path = next(conda_prefix.rglob("snpEff.jar"), None)

    # If not found, install snpEff
    if not jar_path or not jar_path.exists():
        print("⚠️ snpEff.jar not found in current conda environment. Attempting to install via conda...")
        try:
            run_cmd("conda install -y -c bioconda snpeff")
            jar_path = next(conda_prefix.rglob("snpEff.jar"), None)
        except Exception as e:
            logging.error(f"Failed to install snpEff: {e}")
            raise FileNotFoundError("Failed to install snpEff via conda.")

    if not jar_path or not jar_path.exists():
        raise FileNotFoundError("❌ snpEff.jar still not found after installation.")

    snpeff_home = jar_path.parent.parent
    config_file = snpeff_home / "snpEff.config"
    data_dir = snpeff_home / "data" / genome_id
    data_dir.mkdir(parents=True, exist_ok=True)

    gbk_file = data_dir / "genes.gbk"
    if not gbk_file.exists():
        logging.info(f"Downloading GenBank file for {genome_id}")
        try:
            handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
            with open(gbk_file, "w") as file:
                file.write(handle.read())
            handle.close()
        except Exception as e:
            logging.error(f"Failed to download GenBank file: {e}")
            return genome_id, None, None

    genome_line = f"{genome_id}.genome : {genome_id}"
    if config_file.exists():
        lines = config_file.read_text().splitlines()
        if not any(l.strip().startswith(f"{genome_id}.genome") for l in lines):
            with open(config_file, "a") as f:
                f.write(f"\n{genome_line}\n")
            logging.info(f"Added genome entry to config: {genome_line}")
    else:
        config_file.write_text(genome_line + "\n")
        logging.info(f"Created config file with entry: {genome_line}")

    try:
        run_cmd(f"java -Xmx4g -jar {jar_path} build -genbank -v {genome_id} -c {config_file}")
        snpeff_db_path = snpeff_home / "data" / genome_id
        if not any(snpeff_db_path.glob("*.bin")):
            raise RuntimeError(f"snpEff database for {genome_id} not created properly in {snpeff_db_path}")
    except Exception as e:
        logging.warning(f"snpEff build failed: {e}")
        return genome_id, None, None
    return genome_id, str(jar_path), str(config_file)

def verify_tools():
    print("\n🔍 Tool Verification:")
    missing = []
    for pkg, cmd in TOOL_MAPPING.items():
        if not shutil.which(cmd):
            missing.append(pkg)
    if missing:
        print(f"⚠ Missing tools: {', '.join(missing)}")
        print("⏳ Attempting conda installation...")
        run_cmd(f"conda install -y -c bioconda -c conda-forge {' '.join(missing)}")
        missing_post = [pkg for pkg, cmd in TOOL_MAPPING.items() if not shutil.which(cmd)]
        if missing_post:
            raise RuntimeError(f"Failed to install: {', '.join(missing_post)}")
    print("✅ All tools verified")

def check_or_build_kraken2_db():
    print("\n🔍 Kraken2 Database:")
    kraken_db_path = Path("./kraken2_viral_db")
    essential_files = ["hash.k2d", "opts.k2d", "taxo.k2d"]

    if all((kraken_db_path / f).exists() for f in essential_files):
        print("✅ Pre-built database verified")
        return kraken_db_path

    print("⚠ Database not found. Building Kraken2 database...")
    kraken_db_path.mkdir(parents=True, exist_ok=True)
    try:
        run_cmd(f"kraken2-build --download-taxonomy --db {kraken_db_path} --use-ftp")
        run_cmd(f"kraken2-build --download-library viral --db {kraken_db_path}")
        run_cmd(f"kraken2-build --build --db {kraken_db_path}")
        run_cmd(f"kraken2-build --clean --db {kraken_db_path} --skip-maps")
        if all((kraken_db_path / f).exists() for f in essential_files):
            print("✅ Database built successfully")
            return kraken_db_path
        raise RuntimeError("Database build failed - missing essential files")
    except Exception as e:
        print("❌ Database build failed. Consider using a pre-built database")
        raise

def extract_sample_name(fastq_path):
    return Path(fastq_path).stem.replace(".fastq", "")

def should_call_variants(fq_path):
    path = Path(fq_path)
    return path.exists() and path.stat().st_size > 0

def summarize_kraken2_report(report_path, sample_name):
    try:
        df = pd.read_csv(report_path, sep="\t", header=None, names=[
            "Percentage", "Reads_covered", "Reads_assigned",
            "Taxonomic_rank", "NCBI_ID", "Taxon_name"])
        df = df[df["Taxonomic_rank"] == "S"]
        df["Taxon_name"] = df["Taxon_name"].str.strip()
        df["Sample"] = sample_name
        return df.groupby(["Sample", "Taxon_name"])[["Percentage", "Reads_covered"]].sum().reset_index()
    except Exception as e:
        logging.error(f"Kraken2 summary failed for {sample_name}: {e}")
        return None

def get_assembly_stats(assembly_fasta, sample, assembly_type):
    try:
        run_cmd(f"samtools faidx {assembly_fasta}")
        fai_file = f"{assembly_fasta}.fai"
        if not Path(fai_file).exists():
            raise FileNotFoundError(f"Index file {fai_file} not found")

        contig_lengths = [int(line.strip().split('\t')[1]) for line in open(fai_file)]
        contig_lengths.sort(reverse=True)
        total_length = sum(contig_lengths)
        num_contigs = len(contig_lengths)

        half_length = total_length / 2
        cumulative_length = 0
        n50, l50 = 0, 0
        for i, length in enumerate(contig_lengths):
            cumulative_length += length
            if cumulative_length >= half_length:
                n50 = length
                l50 = i + 1
                break

        avg_length = total_length / num_contigs
        largest_contig = contig_lengths[0]

        return {
            "Sample": sample,
            "Assembly_Type": assembly_type,
            "Total_Length": total_length,
            "Number_of_Contigs": num_contigs,
            "Largest_Contig": largest_contig,
            "N50": n50,
            "L50": l50,
            "Average_Contig_Length": avg_length
        }
    except Exception as e:
        logging.error(f"Failed to calculate {assembly_type} assembly stats for {sample}: {e}")
        return None

def extract_sample(fq, ref_fasta, sample_dir, sample, genome_id, jar_path, config_file, threads, KRONA_DIR):
    kraken_db = check_or_build_kraken2_db()
    bam = sample_dir / f"{sample}_aligned.bam"
    sorted_bam = sample_dir / f"{sample}_aligned.sorted.bam"
    mapped_fq = sample_dir / f"{sample}_mapped_reads.fastq"
    dedup_fq = sample_dir / f"{sample}_mapped_reads_dedup.fastq"
    ann_vcf = sample_dir / f"{sample}_variants_annotated.vcf"
    vcf = sample_dir / f"{sample}_variants.vcf.gz"
    kraken_report = sample_dir / "kraken2_report.txt"
    kraken_output = sample_dir / "kraken2_output.txt"
    krona_output = KRONA_DIR / f"{sample}_krona.html"
    denovo_dir = sample_dir / "denovo_assembly"
    reference_dir = sample_dir / "reference_assembly"
    variant_txt = sample_dir / f"{sample}_variant_summary.txt"

    try:
        run_cmd(f"minimap2 -ax map-ont -t {threads} {ref_fasta} {fq} | samtools view -Sb - > {bam}")
        run_cmd(f"samtools sort -@ {threads} -o {sorted_bam} {bam}")
        run_cmd(f"samtools index {sorted_bam}")
        run_cmd(f"samtools fastq -F 4 {sorted_bam} > {mapped_fq}")
        run_cmd(f"awk '{{if(NR%4==1){{$0=sprintf(\"@%s\", NR/4)}} print}}' {mapped_fq} > {dedup_fq}")

        read_count = sum(1 for _ in open(dedup_fq)) // 4
        mapped_reads_counts.append({
            "Sample": sample,
            "Extracted_Mapped_Reads": dedup_fq.name,
            "Read_Count": read_count
        })

        run_cmd(f"kraken2 --db {kraken_db} --threads {threads} --report {kraken_report} --output {kraken_output} {dedup_fq}")
        run_cmd(f"cut -f2,3 {kraken_output} | ktImportTaxonomy -m 3 -o {krona_output}")
        sample_report = summarize_kraken2_report(kraken_report, sample)

        if should_call_variants(dedup_fq):
            denovo_dir.mkdir(exist_ok=True)
            run_cmd(f"flye --nano-raw {dedup_fq} --out-dir {denovo_dir} --threads {threads}")
            denovo_assembly_fasta = denovo_dir / "assembly.fasta"
            if denovo_assembly_fasta.exists():
                stats = get_assembly_stats(denovo_assembly_fasta, sample, "de_novo")
                if stats:
                    assembly_stats.append(stats)

            reference_dir.mkdir(exist_ok=True)
            run_cmd(f"medaka_consensus -i {dedup_fq} -d {ref_fasta} -o {reference_dir} -t {threads} -m r941_min_high_g360")
            ref_assembly_fasta = reference_dir / "consensus.fasta"
            if ref_assembly_fasta.exists():
                stats = get_assembly_stats(ref_assembly_fasta, sample, "reference_based")
                if stats:
                    assembly_stats.append(stats)

            run_cmd(f"bcftools mpileup -f {ref_fasta} --threads {threads} {sorted_bam} | bcftools call -mv -Oz -o {vcf} --threads {threads}")
            run_cmd(f"bcftools index {vcf}")
            if jar_path and config_file:
                run_cmd(f"java -Xmx4g -jar {jar_path} ann -noStats -v {genome_id} -c {config_file} {vcf} > {ann_vcf}")
                run_cmd(f"bcftools stats {vcf} > {variant_txt}")

        return sample_report, str(kraken_output)
    except Exception as e:
        logging.error(f"Processing sample {sample} failed: {e}")
        return None, None

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-inputs", required=True, help="Path to file with list of input FASTQ files")
    parser.add_argument("-reference", required=True, help="Reference genome FASTA file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads to use")
    args = parser.parse_args()

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    results_dir = Path(f"results_{timestamp}")
    summary_dir, krona_dir = setup_directories(results_dir)
    multiqc_config = create_multiqc_config(results_dir)

    print("\n🚀 Starting pre-run verification checks")
    ensure_java_version()
    verify_tools()
    check_or_build_kraken2_db()
    genome_id, jar_path, config_file = configure_snpeff(args.reference)
    print("\n✅ All pre-run checks completed successfully")

    fastq_list = [line.strip() for line in open(args.inputs) if line.strip()]
    quality_control(fastq_list, results_dir)

    all_kraken_reports = []
    kraken_output_files = []

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = []
        for fq in fastq_list:
            sample = extract_sample_name(fq)
            sample_dir = results_dir / sample
            sample_dir.mkdir(exist_ok=True)
            futures.append(executor.submit(
                extract_sample,
                fq, args.reference,
                sample_dir, sample,
                genome_id, jar_path, config_file,
                args.threads,
                krona_dir
            ))
        for future in as_completed(futures):
            try:
                report, kraken_output = future.result()
                if report is not None:
                    all_kraken_reports.append(report)
                    if kraken_output:
                        kraken_output_files.append(kraken_output)
            except Exception as e:
                logging.error(f"Threaded task failed: {e}")

    if all_kraken_reports:
        combined_kraken = pd.concat(all_kraken_reports)
        kraken_summary_file = summary_dir / "combined_kraken2_species_abundance.csv"
        combined_kraken.to_csv(kraken_summary_file, index=False)
        print(f"✅ Saved combined Kraken2 report to {kraken_summary_file}")


    if len(kraken_output_files) > 1:
        combined_krona_output = krona_dir / "combined_krona.html"
        combined_kraken_temp = krona_dir / "combined_kraken_temp.txt"
        with open(combined_kraken_temp, 'w') as outfile:
            for fname in kraken_output_files:
                with open(fname) as infile:
                    outfile.write(infile.read())
        try:
            print("\n🌐 Generating combined Krona report...")
            run_cmd(f"ktImportTaxonomy -q 1 -t 2 -c -o {combined_krona_output} {combined_kraken_temp}")
            print(f"✅ Saved combined Krona report to {combined_krona_output}")
        except Exception as e:
            logging.error(f"Failed to generate combined Krona report: {e}")
        finally:
            combined_kraken_temp.unlink(missing_ok=True)

    run_cmd(f"multiqc {results_dir}/qc -o {results_dir}/multiqc --force --ignore violin")

	
    mapped_reads_summary = summary_dir / "mapped_reads_summary.tsv"
    assembly_stats_summary = summary_dir / "assembly_stats_summary.tsv"

    if mapped_reads_counts:
        pd.DataFrame(mapped_reads_counts).to_csv(mapped_reads_summary, sep="\t", index=False)
    if assembly_stats:
        pd.DataFrame(assembly_stats).to_csv(assembly_stats_summary, sep="\t", index=False)

    print(f"\n🎉 Pipeline finished. Results saved in {results_dir}")
if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ Pipeline failed: {str(e)}")
        sys.exit(1)
