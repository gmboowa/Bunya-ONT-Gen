# NOTE: Please activate the conda environment manually before running this script:
# $ conda activate Bunya-ONT-Gen

import os
import subprocess
import shutil
from pathlib import Path
import pandas as pd
import ssl
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import re

# Setup logging
logging.basicConfig(
    filename='pipeline.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

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
assembly_stats = []  # List to store assembly statistics

os.environ["SNPEFF_HOME"] = str(Path.home() / "snpEff")

MAPPED_READS_SUMMARY = Path("./results/mapped_reads_summary.tsv")
ASSEMBLY_STATS_SUMMARY = Path("./results/assembly_stats_summary.tsv")  # Path for assembly stats output

SUMMARY_DIR = Path("./results/summary")
KRONA_DIR = Path("./results/krona")

for dir_path in [SUMMARY_DIR, KRONA_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

def run_cmd(cmd):
    logging.info(f"Running command: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}\n{e}")
        raise

def ensure_java_version():
    try:
        result = subprocess.run(["java", "-version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        version_match = re.search(r'version "(?P<ver>\d+\.\d+)', result.stdout)
        if version_match and float(version_match.group("ver").split('.')[0]) < 21:
            print("⚠ Java version is less than 21. snpEff may fail. Please upgrade to Java 21+.")
    except FileNotFoundError:
        print(" Java is not installed or not found in PATH.")
        raise

def extract_sample_name(fastq_path):
    return Path(fastq_path).stem.replace(".fastq", "")

def should_call_variants(fq_path):
    path = Path(fq_path)
    return path.exists() and path.stat().st_size > 0

def setup_directories():
    for subdir in ["qc", "nanoplot", "summary", "krona", "multiqc"]:
        Path(f"./results/{subdir}").mkdir(parents=True, exist_ok=True)

def ensure_dependencies():
    logging.info("Checking dependencies")
    missing = [pkg for pkg, cmd in TOOL_MAPPING.items() if shutil.which(cmd) is None]
    if missing:
        run_cmd(f"conda install -y -c bioconda -c conda-forge {' '.join(missing)}")

def check_or_build_kraken2_db():
    kraken_db_path = Path("./kraken2_viral_db")
    if not (kraken_db_path / "taxo.k2d").exists():
        print("⚠ Kraken2 DB not found. Attempting to build...")
        kraken_db_path.mkdir(parents=True, exist_ok=True)
        run_cmd(f"kraken2-build --download-taxonomy --db {kraken_db_path}")
        run_cmd(f"yes | kraken2-build --download-library viral --db {kraken_db_path}")
        run_cmd(f"kraken2-build --build --db {kraken_db_path}")
        run_cmd(f"kraken2-build --clean --db {kraken_db_path}")
    return kraken_db_path

def quality_control(fastq_list):
    for fq in fastq_list:
        run_cmd(f"fastqc {fq} -o ./results/qc")
        run_cmd(f"nanoplot --fastq {fq} -o ./results/nanoplot")

def summarize_kraken2_report(report_path, summary_dir):
    sample_name = Path(report_path).parent.name
    summary_dir = Path(summary_dir)
    summary_dir.mkdir(parents=True, exist_ok=True)
    try:
        df = pd.read_csv(report_path, sep="\t", header=None, names=[
            "Percentage", "Reads_covered", "Reads_assigned",
            "Taxonomic_rank", "NCBI_ID", "Taxon_name"])
        df = df[df["Taxonomic_rank"] == "S"]
        df["Taxon_name"] = df["Taxon_name"].str.strip()
        df = df.groupby("Taxon_name")[["Percentage", "Reads_covered"]].sum().sort_values(by="Reads_covered", ascending=False)
        df.to_csv(summary_dir / f"{sample_name}_species_abundance_summary.csv")
    except Exception as e:
        logging.error(f"Kraken2 summary failed for {sample_name}: {e}")

def get_assembly_stats(assembly_fasta, sample, assembly_type):
    """Calculate assembly statistics using samtools and custom Python code"""
    try:
        # Get basic stats with samtools
        stats_cmd = f"samtools faidx {assembly_fasta}"
        run_cmd(stats_cmd)
        
        # Read the .fai file to get contig lengths
        fai_file = f"{assembly_fasta}.fai"
        if not Path(fai_file).exists():
            raise FileNotFoundError(f"Index file {fai_file} not found")
            
        contig_lengths = []
        with open(fai_file) as f:
            for line in f:
                parts = line.strip().split('\t')
                contig_lengths.append(int(parts[1]))
        
        if not contig_lengths:
            return None
            
        contig_lengths.sort(reverse=True)
        total_length = sum(contig_lengths)
        num_contigs = len(contig_lengths)
        
        # Calculate N50
        half_length = total_length / 2
        cumulative_length = 0
        n50 = 0
        for length in contig_lengths:
            cumulative_length += length
            if cumulative_length >= half_length:
                n50 = length
                break
                
        # Calculate L50
        l50 = 0
        cumulative_length = 0
        for length in contig_lengths:
            cumulative_length += length
            l50 += 1
            if cumulative_length >= half_length:
                break
                
        # Calculate average contig length
        avg_length = total_length / num_contigs
        
        # Find largest contig
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

def configure_snpeff(ref_fasta):
    genome_id = Path(ref_fasta).stem
    snpeff_home = Path(os.environ.get("SNPEFF_HOME", str(Path.home() / "snpEff")))
    jar_path = next(snpeff_home.rglob("snpEff.jar"), None)

    if not jar_path:
        logging.error("snpEff.jar not found")
        raise FileNotFoundError("snpEff.jar not found")

    data_dir = snpeff_home / "data" / genome_id
    config_file = snpeff_home / "snpEff.config"
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

def extract_sample(fq, ref_fasta, sample_dir, sample, genome_id, jar_path, config_file):
    kraken_db = check_or_build_kraken2_db()

    bam = sample_dir / f"{sample}_aligned.bam"
    sorted_bam = sample_dir / f"{sample}_aligned.sorted.bam"
    mapped_fq = sample_dir / f"{sample}_mapped_reads.fastq"
    dedup_fq = sample_dir / f"{sample}_mapped_reads_dedup.fastq"
    ann_vcf = sample_dir / f"{sample}_variants_annotated.vcf"
    vcf = sample_dir / f"{sample}_variants.vcf.gz"
    kraken_report = sample_dir / "kraken2_report.txt"
    
    # Create distinct directories for each assembly type
    denovo_dir = sample_dir / "denovo_assembly"
    reference_dir = sample_dir / "reference_assembly"
    variant_txt = sample_dir / f"{sample}_variant_summary.txt"

    try:
        # Mapping and read processing
        run_cmd(f"minimap2 -ax map-ont {ref_fasta} {fq} | samtools view -Sb - > {bam}")
        run_cmd(f"samtools sort -o {sorted_bam} {bam}")
        run_cmd(f"samtools index {sorted_bam}")
        run_cmd(f"samtools fastq -F 4 {sorted_bam} > {mapped_fq}")
        run_cmd(f"awk '{{if(NR%4==1){{$0=sprintf(\"@%s\", NR/4)}} print}}' {mapped_fq} > {dedup_fq}")

        read_count = sum(1 for _ in open(dedup_fq)) // 4
        mapped_reads_counts.append({"Sample": sample, "Extracted_Mapped_Reads": dedup_fq.name, "Read_Count": read_count})

        # Taxonomic classification
        run_cmd(f"kraken2 --db {kraken_db} --threads 4 --report {kraken_report} --output {sample_dir}/kraken2_output.txt {dedup_fq}")
        summarize_kraken2_report(kraken_report, "./results/summary")

        if should_call_variants(dedup_fq):
            # Step 1: De novo assembly in its own directory
            print(f"\n🔬 Performing DE NOVO assembly for {sample}")
            denovo_dir.mkdir(exist_ok=True)
            run_cmd(f"flye --nano-raw {dedup_fq} --out-dir {denovo_dir} --threads 4")
            
            # Get de novo assembly statistics
            denovo_assembly_fasta = denovo_dir / "assembly.fasta"
            if denovo_assembly_fasta.exists():
                stats = get_assembly_stats(denovo_assembly_fasta, sample, "de_novo")
                if stats:
                    assembly_stats.append(stats)
            
            # Step 2: Reference-based assembly in its own directory
            print(f"\n🧬 Performing REFERENCE-BASED assembly for {sample}")
            reference_dir.mkdir(exist_ok=True)
            run_cmd(f"medaka_consensus -i {dedup_fq} -d {ref_fasta} -o {reference_dir} -t 4 -m r941_min_high_g360")
            
            # Get reference-based assembly statistics
            ref_assembly_fasta = reference_dir / "consensus.fasta"
            if ref_assembly_fasta.exists():
                stats = get_assembly_stats(ref_assembly_fasta, sample, "reference_based")
                if stats:
                    assembly_stats.append(stats)
            
            # Variant calling (keeps original location)
            run_cmd(f"bcftools mpileup -f {ref_fasta} {sorted_bam} | bcftools call -mv -Oz -o {vcf}")
            run_cmd(f"bcftools index {vcf}")
            if jar_path and config_file:
                run_cmd(f"java -Xmx4g -jar {jar_path} ann -noStats -v {genome_id} -c {config_file} {vcf} > {ann_vcf}")
                run_cmd(f"bcftools stats {vcf} > {variant_txt}")
    except Exception as e:
        logging.error(f"Processing sample {sample} failed: {e}")

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-inputs", required=True)
    parser.add_argument("-reference", required=True)
    args = parser.parse_args()

    ensure_java_version()
    setup_directories()
    ensure_dependencies()

    fastq_list = [line.strip() for line in open(args.inputs) if line.strip()]
    quality_control(fastq_list)
    genome_id, jar_path, config_file = configure_snpeff(args.reference)

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = []
        for fq in fastq_list:
            sample = extract_sample_name(fq)
            sample_dir = Path("./results") / sample
            sample_dir.mkdir(exist_ok=True)
            futures.append(executor.submit(extract_sample, fq, args.reference, sample_dir, sample, genome_id, jar_path, config_file))
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logging.error(f"Threaded task failed: {e}")

    run_cmd("multiqc ./results/qc -o ./results/multiqc --force")

    if mapped_reads_counts:
        pd.DataFrame(mapped_reads_counts).to_csv(MAPPED_READS_SUMMARY, sep="\t", index=False)
    
    # Save assembly statistics
    if assembly_stats:
        pd.DataFrame(assembly_stats).to_csv(ASSEMBLY_STATS_SUMMARY, sep="\t", index=False)

    print("\n Pipeline finished. Results saved in ./results")

if __name__ == "__main__":
    main()
