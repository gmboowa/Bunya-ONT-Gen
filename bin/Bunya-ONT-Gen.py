#!/usr/bin/env python3
# NOTE: Please activate the conda environment manually before running this script:
# $ conda activate BunyaGen

import os
import sys
import subprocess
import shutil
from pathlib import Path
import pandas as pd
import ssl
from Bio import Entrez
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import re
import glob
import tarfile
import urllib.request

# ----------------------------- Logging ------------------------------------------
logging.basicConfig(
    filename='pipeline.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# ----------------------------- NCBI / SSL ---------------------------------------
ssl._create_default_https_context = ssl._create_unverified_context
Entrez.email = "gmboowa@gmail.com"

# ----------------------------- Tools list ---------------------------------------
BIOCONDA_TOOLS = [
    "fastqc", "nanoplot", "minimap2", "samtools", "bcftools", "medaka",
    "multiqc", "spades.py", "kraken2", "mafft", "fasttree", "seqtk", "flye", "krona", "snpEff", "seqkit"
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
    "snpEff": "snpEff",
    "seqkit": "seqkit",
}

# ----------------------------- Global storage -----------------------------------
mapped_reads_counts = []
assembly_stats = []  # List to store assembly statistics

# ----------------------------- Small helpers ------------------------------------
def run_cmd(cmd, cwd=None):
    """Run a shell command with optional working directory; log & raise on error."""
    logging.info(f"Running command: {cmd} (cwd={cwd})")
    try:
        subprocess.run(cmd, shell=True, check=True, cwd=cwd)
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {cmd}\n{e}")
        raise

def run_cmd_capture(cmd, cwd=None) -> str:
    """Run a command and return stdout; raise on error."""
    logging.info(f"Running (capture) command: {cmd} (cwd={cwd})")
    out = subprocess.check_output(cmd, shell=True, text=True, stderr=subprocess.STDOUT, cwd=cwd)
    return out

def ensure_java_version():
    try:
        result = subprocess.run(["java", "-version"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        m = re.search(r'version\s+"(\d+)', result.stdout)
        if m and int(m.group(1)) < 21:
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
    missing_execs = []
    for key, exe in TOOL_MAPPING.items():
        if shutil.which(exe) is None:
            missing_execs.append(key)
    if missing_execs:
        # map keys to conda package names where needed
        pkgs = []
        for key in missing_execs:
            if key == "spades":
                pkgs.append("spades")
            else:
                pkgs.append(key)
        print(f"Installing missing tools via conda: {' '.join(pkgs)}")
        run_cmd(f"conda install -y -c bioconda -c conda-forge {' '.join(pkgs)}")

# ----------------------------- snpEff helpers -----------------------------------
def _first_existing(*paths) -> str | None:
    for p in paths:
        if p and Path(p).exists():
            return str(p)
    return None

def find_snpeff_paths() -> tuple[str | None, str | None, bool]:
    """
    Returns (snpeff_invoker, snpeff_config, uses_cli).
    """
    env_jar = os.environ.get("SNPEFF_JAR")
    env_cfg = os.environ.get("SNPEFF_CONFIG")
    if env_jar and Path(env_jar).exists():
        cfg = env_cfg if env_cfg and Path(env_cfg).exists() else str(Path(env_jar).with_name("snpEff.config"))
        return env_jar, cfg if Path(cfg).exists() else None, False

    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        candidates = sorted(glob.glob(os.path.join(conda_prefix, "share", "snpeff*", "snpEff.jar")))
        if candidates:
            jar = candidates[0]
            cfg = str(Path(jar).with_name("snpEff.config"))
            return jar, cfg if Path(cfg).exists() else None, False

    home_jar = Path.home() / "snpEff" / "snpEff.jar"
    home_cfg = Path.home() / "snpEff" / "snpEff.config"
    if home_jar.exists():
        return str(home_jar), str(home_cfg) if home_cfg.exists() else None, False

    cli = shutil.which("snpEff")
    if cli:
        cfg = None
        if conda_prefix:
            maybe_cfg = sorted(glob.glob(os.path.join(conda_prefix, "share", "snpeff*", "snpEff.config")))
            if maybe_cfg:
                cfg = maybe_cfg[0]
        return cli, cfg, True

    return None, None, False

def snpeff_build(invoker: str, uses_cli: bool, genome_id: str, config_file: str | None) -> None:
    if uses_cli:
        cmd = f"snpEff build -genbank -v {genome_id}"
        if config_file:
            cmd += f" -c {config_file}"
        run_cmd(cmd)
    else:
        cmd = f"java -Xmx4g -jar {invoker} build -genbank -v {genome_id}"
        if config_file:
            cmd += f" -c {config_file}"
        run_cmd(cmd)

def snpeff_annotate(invoker: str, uses_cli: bool, genome_id: str, config_file: str | None, vcf_in: str, out_vcf: str) -> None:
    if uses_cli:
        cmd = f"snpEff ann -noStats -v {genome_id}"
        if config_file:
            cmd += f" -c {config_file}"
        cmd += f" {vcf_in} > {out_vcf}"
        run_cmd(cmd)
    else:
        cmd = f"java -Xmx4g -jar {invoker} ann -noStats -v {genome_id}"
        if config_file:
            cmd += f" -c {config_file}"
        cmd += f" {vcf_in} > {out_vcf}"
        run_cmd(cmd)

def _looks_like_accession(genome_id: str) -> bool:
    return re.match(r"^[A-Z]{1,6}_?\d+(\.\d+)?$", genome_id, re.IGNORECASE) is not None

def configure_snpeff(ref_fasta):
    """
    Build a local, writable snpEff workspace in ./snpeff_work and attempt DB build
    if a GBK is available for the reference ID.
    """
    genome_id = Path(ref_fasta).stem

    invoker, _, uses_cli = find_snpeff_paths()
    if not invoker:
        logging.error("snpEff not found (jar or CLI). Install snpEff or set $SNPEFF_JAR/$SNPEFF_CONFIG.")
        return genome_id, None, None, False

    work_dir = Path.cwd() / "snpeff_work"
    data_root = work_dir
    config_file = work_dir / "snpEff.config"
    data_dir = data_root / "data" / genome_id
    data_dir.mkdir(parents=True, exist_ok=True)

    lines = []
    if config_file.exists():
        lines = config_file.read_text().splitlines()
    if not any(l.strip().lower().startswith("data.dir") for l in lines):
        lines.append(f"data.dir = {data_root}")
    genome_line = f"{genome_id}.genome : {genome_id}"
    if not any(l.strip().startswith(f"{genome_id}.genome") for l in lines):
        lines.append(genome_line)
    config_file.write_text("\n".join(lines) + "\n")

    sequences_fa = data_dir / "sequences.fa"
    try:
        if sequences_fa.exists() or sequences_fa.is_symlink():
            sequences_fa.unlink()
    except Exception:
        pass
    try:
        sequences_fa.symlink_to(Path(ref_fasta).resolve())
    except Exception:
        shutil.copy2(ref_fasta, sequences_fa)

    gbk_file = data_dir / "genes.gbk"
    if not gbk_file.exists() and _looks_like_accession(genome_id):
        logging.info(f"Attempting to download GenBank for {genome_id}")
        try:
            handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
            txt = handle.read()
            handle.close()
            if "LOCUS" in txt and "ORIGIN" in txt:
                gbk_file.write_text(txt)
            else:
                logging.warning(f"Downloaded GBK for {genome_id} seems incomplete; skipping snpEff build.")
        except Exception as e:
            logging.warning(f"Could not download GBK for {genome_id}: {e}")

    if gbk_file.exists():
        try:
            snpeff_build(invoker, uses_cli, genome_id, str(config_file))
            if not any((data_dir).glob("*.bin")):
                logging.warning(f"snpEff database for {genome_id} not created in {data_dir}; continuing without annotation.")
        except Exception as e:
            logging.warning(f"snpEff build failed for {genome_id}: {e}")
    else:
        logging.info(f"No GBK available for {genome_id}; snpEff build skipped. Variant annotation will be skipped.")

    return genome_id, invoker, str(config_file), uses_cli

# --------- Krona taxonomy helpers (FULLY PATCHED) --------------------------------
def _taxonomy_present(d: Path) -> bool:
    """Krona mainly expects taxonomy.tab; accept that OR names.dmp+nodes.dmp."""
    if not d.exists():
        return False
    if (d / "taxonomy.tab").exists():
        return True
    return (d / "names.dmp").exists() and (d / "nodes.dmp").exists()

def _copy_tree(src: Path, dst: Path):
    dst.mkdir(parents=True, exist_ok=True)
    for item in src.rglob("*"):
        rel = item.relative_to(src)
        out = dst / rel
        if item.is_dir():
            out.mkdir(parents=True, exist_ok=True)
        else:
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(item, out)

def _promote_taxdump_files(root: Path) -> bool:
    """If names.dmp/nodes.dmp ended up in a subfolder (e.g. taxdump/), move to root."""
    moved = False
    for item in root.rglob("names.dmp"):
        nodes = item.parent / "nodes.dmp"
        if nodes.exists():
            for f in (item, nodes):
                dest = root / f.name
                try:
                    if dest.exists():
                        dest.unlink()
                except Exception:
                    pass
                shutil.move(str(f), str(dest))
            moved = True
            break
    return moved

def _download_ncbi_taxdump(dest_dir: Path):
    """Fetch NCBI taxdump and extract names.dmp/nodes.dmp into dest_dir."""
    url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    tmp_tgz = dest_dir / "taxdump.tar.gz"
    dest_dir.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as r, open(tmp_tgz, "wb") as f:
        shutil.copyfileobj(r, f)
    with tarfile.open(tmp_tgz, "r:gz") as tf:
        for member in tf.getmembers():
            base = Path(member.name).name
            if base in {"names.dmp", "nodes.dmp"}:
                tf.extract(member, path=dest_dir)
    try:
        tmp_tgz.unlink()
    except Exception:
        pass
    _promote_taxdump_files(dest_dir)

def _build_taxonomy_tab_from_dump(dest_dir: Path):
    """
    Build minimal taxonomy.tab from names.dmp + nodes.dmp:
      columns: taxid \t parent \t rank \t name
      use only 'scientific name' from names.dmp
    """
    names = dest_dir / "names.dmp"
    nodes = dest_dir / "nodes.dmp"
    out_tab = dest_dir / "taxonomy.tab"
    if not (names.exists() and nodes.exists()):
        return

    sci_name = {}
    with open(names, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                sci_name[parts[0]] = parts[1]

    with open(nodes, "r", encoding="utf-8", errors="ignore") as nf, \
         open(out_tab, "w", encoding="utf-8") as out:
        for line in nf:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 3:
                tid, parent, rank = parts[0], parts[1], parts[2]
                name = sci_name.get(tid, "")
                out.write(f"{tid}\t{parent}\t{rank}\t{name}\n")

def ensure_krona_ready():
    """
    Ensure Krona tools and taxonomy are usable.
    Priority:
      1) Use ./krona_taxonomy with taxonomy.tab if present
      2) Try Krona updater in that folder
      3) If only names/nodes exist, build taxonomy.tab locally
      4) Copy conda taxonomy if available
      5) Download NCBI taxdump and build taxonomy.tab
    """
    # Tools
    kt_import_tax = shutil.which("ktImportTaxonomy") or (Path(os.environ.get("CONDA_PREFIX",""))/"bin"/"ktImportTaxonomy")
    kt_import_text = shutil.which("ktImportText") or (Path(os.environ.get("CONDA_PREFIX",""))/"bin"/"ktImportText")
    if not (shutil.which(str(kt_import_tax)) or shutil.which(str(kt_import_text))):
        print("Installing Krona via conda...")
        run_cmd("conda install -y -c bioconda -c conda-forge krona")
        kt_import_tax = shutil.which("ktImportTaxonomy")
        kt_import_text = shutil.which("ktImportText")
        if not (kt_import_tax or kt_import_text):
            raise RuntimeError("Krona tools not found after installation.")

    local_tax = Path.cwd() / "krona_taxonomy"
    local_tax.mkdir(parents=True, exist_ok=True)
    os.environ["KRONA_TAXONOMY"] = str(local_tax)

    # If already OK
    if _taxonomy_present(local_tax):
        return {"ktImportTaxonomy": str(kt_import_tax) if kt_import_tax else None,
                "ktImportText": str(kt_import_text) if kt_import_text else None,
                "taxonomy_dir": str(local_tax)}

    # Try Krona updater (without flags for compatibility)
    updater = (shutil.which("ktUpdateTaxonomy") or
               shutil.which("ktUpdateTaxonomy.sh") or
               shutil.which("updateTaxonomy.sh"))
    if updater:
        try:
            run_cmd(updater, cwd=local_tax)
        except Exception as e:
            logging.warning(f"Krona updater failed in {local_tax}: {e}")

    # Promote and build taxonomy.tab if needed
    if not (local_tax / "taxonomy.tab").exists():
        _promote_taxdump_files(local_tax)
        if (local_tax / "names.dmp").exists() and (local_tax / "nodes.dmp").exists():
            _build_taxonomy_tab_from_dump(local_tax)

    if (local_tax / "taxonomy.tab").exists():
        print(f"✅ Krona taxonomy ready in {local_tax}")
        return {"ktImportTaxonomy": str(kt_import_tax) if kt_import_tax else None,
                "ktImportText": str(kt_import_text) if kt_import_text else None,
                "taxonomy_dir": str(local_tax)}

    # Copy from conda if present
    for d in [Path(os.environ.get("CONDA_PREFIX",""))/"opt"/"krona"/"taxonomy",
              Path(os.environ.get("CONDA_PREFIX",""))/"share"/"krona"/"taxonomy"]:
        if d and d.exists():
            try:
                _copy_tree(d, local_tax)
                if not (local_tax / "taxonomy.tab").exists():
                    _promote_taxdump_files(local_tax)
                    if (local_tax / "names.dmp").exists() and (local_tax / "nodes.dmp").exists():
                        _build_taxonomy_tab_from_dump(local_tax)
                if (local_tax / "taxonomy.tab").exists():
                    print(f"✅ Krona taxonomy ready in {local_tax} (copied from {d})")
                    return {"ktImportTaxonomy": str(kt_import_tax) if kt_import_tax else None,
                            "ktImportText": str(kt_import_text) if kt_import_text else None,
                            "taxonomy_dir": str(local_tax)}
            except Exception as e:
                logging.warning(f"Copying Krona taxonomy from {d} failed: {e}")

    # Final fallback: download and build taxonomy.tab
    try:
        _download_ncbi_taxdump(local_tax)
        if (local_tax / "names.dmp").exists() and (local_tax / "nodes.dmp").exists():
            _build_taxonomy_tab_from_dump(local_tax)
    except Exception as e:
        logging.warning(f"Direct taxdump download failed: {e}")

    if (local_tax / "taxonomy.tab").exists():
        print(f"✅ Krona taxonomy ready in {local_tax} (built from NCBI dumps)")
    else:
        print(f"⚠ Krona taxonomy not present in {local_tax}. Krona charts may warn/fail.")

    return {"ktImportTaxonomy": str(kt_import_tax) if kt_import_tax else None,
            "ktImportText": str(kt_import_text) if kt_import_text else None,
            "taxonomy_dir": str(local_tax)}

# ----------------------------- QC / Kraken2 --------------------------------------
MAPPED_READS_SUMMARY = Path("./results/mapped_reads_summary.tsv")
ASSEMBLY_STATS_SUMMARY = Path("./results/assembly_stats_summary.tsv")

SUMMARY_DIR = Path("./results/summary")
KRONA_DIR = Path("./results/krona")
for dir_path in [SUMMARY_DIR, KRONA_DIR]:
    dir_path.mkdir(parents=True, exist_ok=True)

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
        df["Taxon_name"] = df["Taxon_name"].str.strip()
        df = df.groupby("Taxon_name")[["Percentage", "Reads_covered"]].sum().sort_values(by="Reads_covered", ascending=False)
        df.to_csv(summary_dir / f"{sample_name}_species_abundance_summary.csv")
    except Exception as e:
        logging.error(f"Kraken2 summary failed for {sample_name}: {e}")

def _build_krona_text_from_kraken_report(report_path: Path, out_txt: Path):
    """
    Create a simple 2-column Krona text input:
      count<tab>Taxon name
    """
    rows = []
    with open(report_path, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            perc, cov, assigned, rank, taxid, name = parts[:6]
            name = name.strip()
            try:
                count = int(cov)
            except Exception:
                continue
            if count > 0 and name:
                rows.append((count, name))
    # write
    with open(out_txt, "w") as out:
        for c, n in rows:
            out.write(f"{c}\t{n}\n")

def make_krona_chart(kt_tools: dict, sample_dir: Path, sample: str):
    """
    Build a Krona chart from the Kraken2 report using ktImportText.
    """
    report = sample_dir / "kraken2_report.txt"
    if not report.exists():
        return
    krona_txt = sample_dir / f"{sample}.krona.txt"
    _build_krona_text_from_kraken_report(report, krona_txt)

    out_html = KRONA_DIR / f"{sample}.krona.html"
    kt_import_text = kt_tools.get("ktImportText") or shutil.which("ktImportText")
    if not kt_import_text:
        logging.warning("ktImportText not found; skipping Krona chart.")
        return
    # Use KRONA_TAXONOMY env prepared by ensure_krona_ready()
    cmd = f'"{kt_import_text}" -o "{out_html}" "{krona_txt}"'
    try:
        run_cmd(cmd)
    except Exception as e:
        logging.warning(f"Krona chart generation failed for {sample}: {e}")

# ----------------------------- Assembly stats ------------------------------------
def get_assembly_stats(assembly_fasta, sample, assembly_type):
    """Calculate assembly statistics using samtools and custom Python code"""
    try:
        stats_cmd = f"samtools faidx {assembly_fasta}"
        run_cmd(stats_cmd)
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

        half_length = total_length / 2
        cumulative_length = 0
        n50 = 0
        for length in contig_lengths:
            cumulative_length += length
            if cumulative_length >= half_length:
                n50 = length
                break

        l50 = 0
        cumulative_length = 0
        for length in contig_lengths:
            cumulative_length += length
            l50 += 1
            if cumulative_length >= half_length:
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

# ------------------ Adaptive de novo helpers (Flye / SPAdes) ---------------------
def _flye_assembled_ok(denovo_dir: Path) -> bool:
    asm = denovo_dir / "assembly.fasta"
    if not asm.exists() or asm.stat().st_size == 0:
        return False
    try:
        with open(asm, "r") as fh:
            for line in fh:
                if line.startswith(">"):
                    return True
    except Exception:
        return False
    return False

def _run_flye_adaptive(clean_fa: Path, denovo_dir: Path, threads: int = 4) -> bool:
    """
    Try Flye with allowed decreasing --min-overlap values (1000..10000).
    Returns True if an assembly was produced; False otherwise.
    """
    denovo_dir.mkdir(exist_ok=True)
    for ovlp in (4000, 3000, 2000, 1500, 1200, 1000):
        try:
            logging.info(f"Trying Flye with --min-overlap {ovlp}")
            run_cmd(f"flye --nano-raw {clean_fa} --out-dir {denovo_dir} --threads {threads} --min-overlap {ovlp}")
            if _flye_assembled_ok(denovo_dir):
                logging.info(f"Flye succeeded with --min-overlap {ovlp}")
                return True
            else:
                logging.warning(f"Flye completed but no disjointigs with --min-overlap {ovlp}")
        except Exception as e:
            logging.warning(f"Flye failed with --min-overlap {ovlp}: {e}")
        # Clean directory for next attempt
        try:
            for p in denovo_dir.iterdir():
                if p.is_file():
                    p.unlink()
                elif p.is_dir():
                    shutil.rmtree(p, ignore_errors=True)
        except Exception:
            pass
    logging.error("All Flye attempts failed or produced no contigs.")
    return False

def _seqkit_stats_fastq(fq_path: Path) -> dict:
    """
    Use seqkit stats -T to get read count and N50. Returns dict with keys:
      reads, sum_len, avg_len, n50
    If seqkit is missing or parsing fails, returns empty dict.
    """
    try:
        tsv = run_cmd_capture(f"seqkit stats -T {fq_path}")
        lines = [ln for ln in tsv.strip().splitlines() if ln.strip()]
        if len(lines) >= 2:
            header = lines[0].split('\t')
            vals = lines[1].split('\t')
            d = dict(zip(header, vals))
            out = {}
            out["reads"] = int(d.get("num_seqs", "0"))
            out["sum_len"] = int(d.get("sum_len", "0"))
            out["avg_len"] = float(d.get("avg_len", "0")) if d.get("avg_len") else 0.0
            n50_val = d.get("N50", "0").replace(",", "")
            out["n50"] = int(float(n50_val)) if n50_val else 0
            return out
    except Exception as e:
        logging.warning(f"Could not get seqkit stats for {fq_path}: {e}")
    return {}

def _run_spades_fallback(dedup_fq: Path, denovo_dir: Path) -> bool:
    """
    Run SPAdes for short/sparse reads. Returns True if contigs produced.
    """
    try:
        denovo_dir.mkdir(exist_ok=True)
        run_cmd(f"spades.py --only-assembler -s {dedup_fq} -o {denovo_dir} -k 21,33,55")
        contigs = denovo_dir / "contigs.fasta"
        if contigs.exists() and contigs.stat().st_size > 0:
            asm = denovo_dir / "assembly.fasta"
            try:
                if asm.exists():
                    asm.unlink()
            except Exception:
                pass
            try:
                asm.symlink_to(contigs.resolve())
            except Exception:
                shutil.copy2(contigs, asm)
            return True
        logging.warning("SPAdes finished but no contigs.fasta produced.")
    except Exception as e:
        logging.warning(f"SPAdes fallback failed: {e}")
    return False

# ----------------------------- Core per-sample work ------------------------------
def extract_sample(fq, ref_fasta, sample_dir, sample, genome_id, snpeff_invoker, snpeff_config, snpeff_uses_cli, kt_tools):
    kraken_db = check_or_build_kraken2_db()

    bam = sample_dir / f"{sample}_aligned.bam"
    sorted_bam = sample_dir / f"{sample}_aligned.sorted.bam"
    mapped_fq = sample_dir / f"{sample}_mapped_reads.fastq"
    dedup_fq = sample_dir / f"{sample}_mapped_reads_dedup.fastq"
    ann_vcf = sample_dir / f"{sample}_variants_annotated.vcf"
    vcf = sample_dir / f"{sample}_variants.vcf.gz"
    kraken_report = sample_dir / "kraken2_report.txt"

    denovo_dir = sample_dir / "denovo_assembly"
    reference_dir = sample_dir / "reference_assembly"
    variant_txt = sample_dir / f"{sample}_variant_summary.txt"

    try:
        # Mapping and read processing
        run_cmd(f"minimap2 -ax map-ont {ref_fasta} {fq} | samtools view -Sb - > {bam}")
        run_cmd(f"samtools sort -o {sorted_bam} {bam}")
        run_cmd(f"samtools index {sorted_bam}")
        # drop unmapped(4), secondary(256), supplementary(2048) = 2308
        run_cmd(f"samtools fastq -F 2308 {sorted_bam} > {mapped_fq}")
        # deduplicate by sequence to avoid malformed/duplicate entries
        run_cmd(f"seqkit rmdup -s {mapped_fq} > {dedup_fq}")

        read_count = sum(1 for _ in open(dedup_fq)) // 4
        mapped_reads_counts.append({"Sample": sample, "Extracted_Mapped_Reads": dedup_fq.name, "Read_Count": read_count})

        # Taxonomic classification
        run_cmd(f"kraken2 --db {kraken_db} --threads 4 --report {kraken_report} --output {sample_dir}/kraken2_output.txt {dedup_fq}")
        summarize_kraken2_report(kraken_report, "./results/summary")
        make_krona_chart(kt_tools, sample_dir, sample)

        if should_call_variants(dedup_fq):
            # Prepare cleaned FASTA for Flye (remove gaps/zero-length; convert to FASTA)
            clean_fa = sample_dir / f"{sample}_mapped_reads.clean.fasta"
            run_cmd(f"seqkit seq -m 1 -g {dedup_fq} | seqkit fq2fa > {clean_fa}")

            # Decide assembler by N50
            stats = _seqkit_stats_fastq(Path(dedup_fq))
            n50 = stats.get("n50", 0)
            logging.info(f"Read stats for {sample}: {stats}")

            # De novo assembly
            print(f"\n🔬 Performing DE NOVO assembly for {sample}")
            de_novo_ok = False
            if n50 >= 1000:
                de_novo_ok = _run_flye_adaptive(clean_fa, denovo_dir, threads=4)
                if not de_novo_ok:
                    logging.info("Flye unsuccessful; switching to SPAdes fallback.")
                    de_novo_ok = _run_spades_fallback(Path(dedup_fq), denovo_dir)
            else:
                logging.info(f"N50 ({n50}) < 1000; using SPAdes fallback directly.")
                de_novo_ok = _run_spades_fallback(Path(dedup_fq), denovo_dir)

            # De novo assembly statistics
            denovo_assembly_fasta = denovo_dir / "assembly.fasta"
            if de_novo_ok and denovo_assembly_fasta.exists():
                stats = get_assembly_stats(denovo_assembly_fasta, sample, "de_novo")
                if stats:
                    assembly_stats.append(stats)
            elif not de_novo_ok:
                logging.info(f"De novo assembly not produced for {sample}; continuing with reference-based pipeline.")

            # REFERENCE-BASED assembly
            print(f"\n🧬 Performing REFERENCE-BASED assembly for {sample}")
            reference_dir.mkdir(exist_ok=True)
            run_cmd(f"medaka_consensus -i {dedup_fq} -d {ref_fasta} -o {reference_dir} -t 4")

            # Reference-based assembly statistics
            ref_assembly_fasta = reference_dir / "consensus.fasta"
            if ref_assembly_fasta.exists():
                stats = get_assembly_stats(ref_assembly_fasta, sample, "reference_based")
                if stats:
                    assembly_stats.append(stats)

            # Variant calling
            run_cmd(f"bcftools mpileup -f {ref_fasta} {sorted_bam} | bcftools call -mv -Oz -o {vcf}")
            run_cmd(f"bcftools index {vcf}")

            # snpEff annotation only if available
            if snpeff_invoker and snpeff_config and any((Path("snpeff_work") / "data" / genome_id).glob("*.bin")):
                snpeff_annotate(snpeff_invoker, snpeff_uses_cli, genome_id, snpeff_config, str(vcf), str(ann_vcf))
                run_cmd(f"bcftools stats {vcf} > {variant_txt}")
            else:
                logging.info("Skipping snpEff annotation (DB not available).")
    except Exception as e:
        logging.error(f"Processing sample {sample} failed: {e}")

# ------------------ Kraken2 DB auto-detect + build (viral-only) ------------------
def _kraken_db_has_files(db_dir: Path) -> bool:
    return (db_dir / "taxo.k2d").exists() or (db_dir / "hash.k2d").exists()

def check_or_build_kraken2_db() -> Path:
    """
    Locate a viral Kraken2 DB, or build one if not found.
    Search order:
      1) $KRAKEN2_DB_PATH
      2) ./kraken2_viral_db
      3) $CONDA_PREFIX/share/kraken2/viral
      4) ~/kraken2/viral
    If none found, builds a viral-only DB in ./kraken2_viral_db.
    Env:
      KRAKEN_THREADS or KRAKEN2_THREADS (int, default: cpu_count)
      KRAKEN_USE_FTP=1  (add --use-ftp)
      KRAKEN_NO_CLEAN=1 (skip --clean)
    """
    env_db = os.environ.get("KRAKEN2_DB_PATH")
    if env_db:
        db = Path(env_db)
        if _kraken_db_has_files(db):
            logging.info(f"Using Kraken2 DB from $KRAKEN2_DB_PATH: {db}")
            return db
        else:
            logging.warning(f"$KRAKEN2_DB_PATH set to {db}, but DB files not found; will continue searching.")

    local_db = Path("./kraken2_viral_db")
    if _kraken_db_has_files(local_db):
        logging.info(f"Using local Kraken2 DB: {local_db}")
        return local_db

    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        conda_db = Path(conda_prefix) / "share" / "kraken2" / "viral"
        if _kraken_db_has_files(conda_db):
            logging.info(f"Using conda-provided Kraken2 viral DB: {conda_db}")
            return conda_db

    home_db = Path.home() / "kraken2" / "viral"
    if _kraken_db_has_files(home_db):
        logging.info(f"Using home Kraken2 viral DB: {home_db}")
        return home_db

    print("⚠ Kraken2 viral DB not found. Attempting to build...")
    builder = shutil.which("kraken2-build")
    if not builder:
        raise FileNotFoundError("kraken2-build not found on PATH. Please install kraken2 in this environment.")

    threads = (
        os.environ.get("KRAKEN_THREADS")
        or os.environ.get("KRAKEN2_THREADS")
        or str(os.cpu_count() or 2)
    )
    use_ftp = os.environ.get("KRAKEN_USE_FTP", "0") == "1"
    no_clean = os.environ.get("KRAKEN_NO_CLEAN", "0") == "1"

    local_db.mkdir(parents=True, exist_ok=True)
    ftp_opt = "--use-ftp" if use_ftp else ""

    run_cmd(f"{builder} --download-taxonomy --db {local_db} {ftp_opt}".strip())
    run_cmd(f"{builder} --download-library viral --db {local_db} {ftp_opt}".strip())
    run_cmd(f"{builder} --build --db {local_db} --threads {threads}")
    if not no_clean:
        run_cmd(f"{builder} --clean --db {local_db}")

    inspector = shutil.which("kraken2-inspect")
    if inspector:
        try:
            run_cmd(f"{inspector} --db {local_db} | head -n 5")
        except Exception:
            pass

    if not _kraken_db_has_files(local_db):
        raise RuntimeError(f"Kraken2 DB build appears to have failed in {local_db}")
    logging.info(f"Built Kraken2 viral DB at {local_db}")
    return local_db

# ----------------------------- Main ---------------------------------------------
def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-inputs", required=True, help="Text file with paths to FASTQ(.gz)")
    parser.add_argument("-reference", required=True, help="Reference FASTA for mapping/medaka/snpEff")
    args = parser.parse_args()

    ensure_java_version()
    setup_directories()
    ensure_dependencies()

    # Krona readiness (install + taxonomy setup)
    kt_tools = ensure_krona_ready()

    fastq_list = [line.strip() for line in open(args.inputs) if line.strip()]
    quality_control(fastq_list)
    genome_id, snpeff_invoker, snpeff_config, snpeff_uses_cli = configure_snpeff(args.reference)

    with ThreadPoolExecutor(max_workers=2) as executor:
        futures = []
        for fq in fastq_list:
            sample = extract_sample_name(fq)
            sample_dir = Path("./results") / sample
            sample_dir.mkdir(exist_ok=True)
            futures.append(executor.submit(
                extract_sample,
                fq, args.reference, sample_dir, sample,
                genome_id, snpeff_invoker, snpeff_config, snpeff_uses_cli, kt_tools
            ))
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                logging.error(f"Threaded task failed: {e}")

    run_cmd("multiqc ./results/qc -o ./results/multiqc --force")

    if mapped_reads_counts:
        pd.DataFrame(mapped_reads_counts).to_csv(MAPPED_READS_SUMMARY, sep="\t", index=False)
    if assembly_stats:
        pd.DataFrame(assembly_stats).to_csv(ASSEMBLY_STATS_SUMMARY, sep="\t", index=False)

    print("\n Pipeline finished. Results saved in ./results")

if __name__ == "__main__":
    main()
