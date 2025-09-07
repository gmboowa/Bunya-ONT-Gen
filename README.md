# Bunya-ONT-Gen: A pipeline for viral ONT metagenomic analysis

**Bunya-ONT-Gen** is a modular, end-to-end pipeline for characterizing viral genomes from Oxford Nanopore (ONT) metagenomic data. It covers host read removal, QC, assembly, variant/annotation & a **single concatenated phylogeny** built from the three Bunyamwera segments.

---

## Table of Contents
1. [Features](#features)  
2. [Environments & Installation](#environments--installation)  
   - [Genomics environment](#genomics-environment)  
   - [Phylogeny environment](#phylogeny-environment)  
   - [Optional Medaka env (macOS Apple Silicon)](#optional-medaka-env-macos-apple-silicon)  
3. [Host read removal (Hostile workflow)](#host-read-removal-hostile-workflow)
4. [Results layout](#results-layout)  
5. [Whole-genome analysis](#whole-genome-analysis)  
6. [Concatenated phylogeny from VCFs](#concatenated-phylogeny-from-vcfs)  
   - [Script: `build_bunya_concat_tree.sh`](#script-build_bunya_concat_treesh)  
   - [Inputs](#inputs)  
   - [Run examples](#run-examples)  
   - [Outputs](#outputs)  
   - [Interpreting the tree & showing bootstrap in iTOL](#interpreting-the-tree--showing-bootstrap-in-itol)  
7. [Troubleshooting](#troubleshooting)  
8. [Reproducibility notes](#reproducibility-notes)  
9. [License](#license)

---

## Features

- Optimized for **ONT** long-read data.
- **Human host read removal** before downstream analysis.
- QC summary with FastQC/NanoPlot/MultiQC.
- *De novo* assembly & polishing.
- Variant calling & functional annotation.
- **Concatenated phylogeny across the three Bunyamwera segments** (L/M/S) with partitioned models & bootstrap support.

---

## Environments & Installation

### Genomics environment

```bash
# Create and activate main analysis env
conda env create -f Bunya-ONT-Gen.yml
conda activate Bunya-ONT-Gen

# Ensure core tools are present (safe re-install)
conda install -y -c conda-forge -c bioconda \
  pandas biopython openjdk \
  prodigal diamond fastqc nanoplot minimap2 samtools bcftools multiqc \
  spades kraken2 mafft fasttree seqtk seqkit flye krona snpeff medaka
```

### Phylogeny environment

```bash
conda env create -f Bunya_phylogeny.yml
conda activate Bunya_phylogeny

# If the YAML is not available, you can install directly:
mamba install -c bioconda -c conda-forge \
  bcftools samtools mafft seqkit iqtree fasttree entrez-direct
```

### Optional Medaka env (macOS Apple Silicon)

```bash
# Some Apple Silicon systems need a separate env for Medaka
CONDA_SUBDIR=osx-64 mamba create -n medaka_v1 -c bioconda -c conda-forge medaka=1.7.3
```

---

## Host read removal (Hostile workflow)

Remove human host reads (e.g., **T2T-CHM13v2.0**) from **single-end ONT** data using `hostile + minimap2`.

```bash
conda create -n hostile -c bioconda ncbi-datasets-cli minimap2 hostile
conda activate hostile
```

**Required args**
| Arg | Description |
| --- | --- |
| `--fasta` | NCBI accession (e.g., `GCF_009914755.4`) or path to a local FASTA |
| `-i` | Text file listing paths to FASTQ files (one per line) |
| `-o` | Output directory |

**Example**
```bash
python3 hostile_clean_ont_human_minimap2.py \
  --fasta GCF_009914755.4 \
  -i Bunya_ont_sample.list.txt \
  -o Clean
```

`Bunya_ont_sample.list.txt` format:
```
~/sample1.fastq.gz
~/sample2.fastq.gz
```

---

## Whole-genome analysis

Activate the env & run the main pipeline script:

```bash
# Activate
source /opt/homebrew/Caskroom/miniforge/base/bin/activate
conda activate Bunya-ONT-Gen

# Run (example)
python Bunya-ONT-Gen.py \
  -inputs fastq_sample.txt \
  -reference Bunya.reference.fasta
```

`fastq_sample.txt` example:
```
~/Sample_01_ONT_reads.fastq.gz
~/Sample_02_ONT_reads.fastq.gz
~/Sample_03_ONT_reads.fastq.gz
```
---
## Whole-genome analysis results layout

```bash
results/
├── 00_RawDataQC/          # FastQC/NanoPlot reports
├── 01_Assemblies/         # Flye & Medaka outputs
├── 02_Variants/           # VCF files & annotations (snpEff)
├── 03_Phylogeny/
│   ├── alignment.fasta    # MAFFT MSA
│   ├── iqtree_results/    # Tree files with support
│   └── phylogeny.pdf      # ggtree visualization (optional)
├── 04_Taxonomy/           # Kraken2/Krona reports
└── reports/               # MultiQC + summary
```

---

## Concatenated phylogeny from VCFs

You will build a single phylogeny using consensus sequences per segment (**L/M/S**) concatenated per sample, with partitions preserved & bootstrap support.

### Script: `build_bunya_concat_tree.sh`

- **What it does**
  - For each VCF/sample: extract each segment → normalize variants → **consensus FASTA**.
  - Per-segment MAFFT alignment.
  - Concatenate segments by sample name (**robust to record order**).
  - IQ-TREE with **UFboot + SH-aLRT** (bootstrap), partitioned by segment.
  - **Automatic fallback** to FastTree with `-boot` if IQ-TREE fails.

- **Install the dependences**
  ```bash
  mamba install -c bioconda -c conda-forge \
    bcftools samtools mafft seqkit iqtree fasttree entrez-direct
  ```

- **Place the script** at the repo root:
  ```
  ./build_bunya_concat_tree.sh
  ```
  Make it executable:
  ```bash
  chmod +x build_bunya_concat_tree.sh
  ```

### Inputs

1. **VCF list file** (`vcfs.txt`): one **absolute path** to a VCF per line  
   _or_ `sample</>sample.vcf[.gz]`.

   Example (4 samples):
   ```
   ~/Bunyamwera/SRR34843736.vcf
   ~/Bunyamwera/SRR34843737.vcf
   ~/Bunyamwera/SRR34843738.vcf
   ~/Bunyamwera/SRR34843739.vcf
   
   ```

2. **Segment accessions** (recommended to **pin**):  
   - `MF926354.1` (L), `KP795084.1` (M), `KP795093.1` (S)

> If you don’t pass `--segments`, the script auto-detects the top 3 contigs in the first VCF; pinning avoids mismatches.

### Run examples

**With automatic NCBI reference fetch & bootstrap 1000:**
```bash
./build_bunya_concat_tree.sh \
  -l ~/vcfs.txt \
  -o bunya_tree_out_concat \
  --segments MF926354.1,KP795084.1,KP795093.1 \
  --auto-ref \
  --boot 1000
```

**With local FASTAs for L/M/S (headers will be rewritten to match VCF contigs):**
```bash
./build_bunya_concat_tree.sh \
  -l ~/vcfs.txt \
  -o bunya_tree_out_concat_localref \
  --segments MF926354.1,KP795084.1,KP795093.1 \
  --refL L.fa --refM M.fa --refS S.fa \
  --boot 1000
```

**Use IUPAC ambiguity codes in the consensus (instead of `-H 1` haplotype selection):**
```bash
./build_bunya_concat_tree.sh \
  -l ~/vcfs.txt \
  -o bunya_tree_out_concat_iupac \
  --segments MF926354.1,KP795084.1,KP795093.1 \
  --auto-ref \
  --iupac \
  --boot 1000
```

### Outputs

```
bunya_tree_out_concat/
├── aln/
│   ├── MF926354.1.aln.fa
│   ├── KP795084.1.aln.fa
│   ├── KP795093.1.aln.fa
│   ├── concat.aln.fa              # concatenated L+M+S alignment
│   └── concat.partitions          # per-segment partitions
├── cons/                          # per-sample per-segment consensus FASTAs
├── logs/
│   └── concat.iqtree.log          # shows model/UFboot/SH-aLRT runs
├── trees/
│   ├── concat.treefile            # final Newick with support values (IQ-TREE or FastTree)
│   └── concat.contree             # IQ-TREE consensus tree (if created)
└── refs/                          # fetched or provided references
```

---

## Interpreting the tree & showing bootstrap in iTOL

**Files to upload to iTOL**
- `trees/concat.treefile` (primary)  
- Optionally, `trees/concat.contree` (IQ-TREE consensus tree)

**To see support values**
- IQ-TREE: UFboot and SH-aLRT are printed on the same node label (e.g., `95/99`).
- FastTree fallback: node labels are bootstrap percentages.

**iTOL steps**
1. Create a new tree → Upload `concat.treefile`.
2. In **Style** → **Internal node labels**, choose to **show numeric labels**.
3. (Optional) Add a color strip or dataset for sample groups (by sample naming).

> If the tree looks star-like or shows no support, check `logs/concat.iqtree.log` and the `DIAG: concat.aln.fa` line printed by the script. If `unique=1` or `varSites=0`, the samples are identical across L/M/S after filtering/masking.

---


## Live reports & outputs

All static reports & pipeline outputs are published via GitHub Pages:

https://gmboowa.github.io/Bunya-ONT-Gen/


# Custom SnpEff Database for Bunyamwera Orthobunyavirus (macOS: Apple Silicon & Intel)

## Detect your Mac architecture

```bash
uname -m
# arm64  -> Apple Silicon (M1/M2/M3)
# x86_64 -> Intel
```

---

## Apple Silicon (arm64) setup

### <a name="apple-create-conda-env"></a>1) Create Conda env with SnpEff, Java & NCBI EDirect

```bash
# Initialize Conda for this shell (run once, then reopen the terminal)
conda init zsh   # or: conda init bash
```

```bash
# Create & activate env
conda create -n BunyaGen -c bioconda -c conda-forge snpeff openjdk=21 entrez-direct -y
conda activate BunyaGen

# Sanity checks
which snpEff
snpEff -version
esearch -help >/dev/null
```

## Make a writable SnpEff home

```bash
SNPEFF_SRC="$(ls -d "$CONDA_PREFIX"/share/snpeff-* | head -n1)"
cp -R "$SNPEFF_SRC" "$HOME/snpeff"
export SNPEFF_HOME="$HOME/snpeff"

# Persist for your shell (pick ONE of the following depending on your shell)
echo 'export SNPEFF_HOME="$HOME/snpeff"' >> ~/.zshrc
# or:
echo 'export SNPEFF_HOME="$HOME/snpeff"' >> ~/.bashrc
```

## Build a Bunyamwera DB from GenBank

Use the *RefSeq* accessions (L/M/S segments):

- **L**: `NC_001925.1`  
- **M**: `NC_001926.1`  
- **S**: `NC_001927.1`

```bash
DB=bunyamwera_ref

# Register a readable name in the config
echo "$DB.genome : Bunyamwera orthobunyavirus (RefSeq NC_001925/6/7)" >> "$SNPEFF_HOME/snpEff.config"

# Create data folder and fetch GenBank records
mkdir -p "$SNPEFF_HOME/data/$DB"
GBK="$SNPEFF_HOME/data/$DB/genes.gbk"
: > "$GBK"
for ACC in NC_001925.1 NC_001926.1 NC_001927.1; do
  efetch -db nucleotide -format gbwithparts -id "$ACC" >> "$GBK"
done

# Should print '3'
grep -c '^LOCUS' "$GBK"

# Build the database
snpEff build \
  -c "$SNPEFF_HOME/snpEff.config" \
  -dataDir "$SNPEFF_HOME/data" \
  -genbank -v "$DB"

# You should now have:
ls -lh "$SNPEFF_HOME/data/$DB/snpEffectPredictor.bin"
```

## Annotate your VCFs

```bash
snpEff ann \
  -c "$SNPEFF_HOME/snpEff.config" \
  -dataDir "$SNPEFF_HOME/data" \
  -i vcf \
  -s snpeff_bunyamwera.html \
  -csvStats snpeff_bunyamwera.csv \
  "$DB" \
  ~/Sample_variants.vcf \
  > Sample.snpeff.vcf
```

---

## (Optional) Force x86_64 packages on an Apple Silicon Mac

```bash
export CONDA_SUBDIR=osx-64
```

```bash
conda create -n BunyaGen -c bioconda -c conda-forge snpeff openjdk=21 entrez-direct -y
conda activate BunyaGen

which snpEff
snpEff -version
esearch -help >/dev/null
```

## Make a writable SnpEff home (x86_64 env)

```bash
SNPEFF_SRC="$(ls -d "$CONDA_PREFIX"/share/snpeff-* | head -n1)"
cp -R "$SNPEFF_SRC" "$HOME/snpeff"
export SNPEFF_HOME="$HOME/snpeff"
echo 'export SNPEFF_HOME="$HOME/snpeff"' >> ~/.bashrc
```

## Build a Bunyamwera DB from GenBank (x86_64 env)

```bash
DB=bunyamwera_ref
echo "$DB.genome : Bunyamwera orthobunyavirus (RefSeq NC_001925/6/7)" >> "$SNPEFF_HOME/snpEff.config"

mkdir -p "$SNPEFF_HOME/data/$DB"
GBK="$SNPEFF_HOME/data/$DB/genes.gbk"
: > "$GBK"
for ACC in NC_001925.1 NC_001926.1 NC_001927.1; do
  efetch -db nucleotide -format gbwithparts -id "$ACC" >> "$GBK"
done
grep -c '^LOCUS' "$GBK"

snpEff build \
  -c "$SNPEFF_HOME/snpEff.config" \
  -dataDir "$SNPEFF_HOME/data" \
  -genbank -v "$DB"

ls -lh "$SNPEFF_HOME/data/$DB/snpEffectPredictor.bin"
```

### <a name="intel-annotate"></a>4) Annotate your VCFs (either arch)

```bash
# STDIN method
snpEff ann -c "$SNPEFF_HOME/snpEff.config" -dataDir "$SNPEFF_HOME/data" \
  -s snpeff_bunyamwera.html -csvStats snpeff_bunyamwera.csv \
  "$DB" < ~/Sample.vcf > Sample.snpeff.vcf

# or positional arg with input format
snpEff ann -c "$SNPEFF_HOME/snpEff.config" -dataDir "$SNPEFF_HOME/data" \
  -i vcf -s snpeff_bunyamwera.html -csvStats snpeff_bunyamwera.csv \
  "$DB" ~/Sample.vcf > Sample.snpeff.vcf
```

---

## One-shot builder script

Save as `build_snpeff_bunyamwera.sh` and run once per machine.

```bash
#!/usr/bin/env bash
set -euo pipefail

DB="${1:-bunyamwera_ref}"

# Require snpEff + entrez-direct in current conda env
command -v snpEff >/dev/null || { echo "Install snpEff in this env"; exit 1; }
command -v efetch  >/dev/null || { echo "Install entrez-direct in this env"; exit 1; }

SNPEFF_SRC="$(ls -d "$CONDA_PREFIX"/share/snpeff-* | head -n1)"
SNPEFF_HOME="${SNPEFF_HOME:-$HOME/snpeff}"
[[ -d "$SNPEFF_HOME" ]] || cp -R "$SNPEFF_SRC" "$SNPEFF_HOME"

grep -q "^$DB\.genome" "$SNPEFF_HOME/snpEff.config" \
  || echo "$DB.genome : Bunyamwera orthobunyavirus (RefSeq NC_001925/6/7)" >> "$SNPEFF_HOME/snpEff.config"

mkdir -p "$SNPEFF_HOME/data/$DB"
GBK="$SNPEFF_HOME/data/$DB/genes.gbk"
: > "$GBK"
for ACC in NC_001925.1 NC_001926.1 NC_001927.1; do
  echo "Fetching $ACC ..."
  efetch -db nucleotide -format gbwithparts -id "$ACC" >> "$GBK"
done

echo "Records: $(grep -c '^LOCUS' "$GBK") (expect 3)"
snpEff build -c "$SNPEFF_HOME/snpEff.config" -dataDir "$SNPEFF_HOME/data" -genbank -v "$DB"

echo "Built $DB at $SNPEFF_HOME/data/$DB"
echo "Annotate with:"
echo "  snpEff ann -c $SNPEFF_HOME/snpEff.config -dataDir $SNPEFF_HOME/data $DB < input.vcf > output.vcf"
```

Make it executable & run:

```bash
chmod +x build_snpeff_bunyamwera.sh
./build_snpeff_bunyamwera.sh
```

---

## Common pitfalls & fixes

- **“Unknown parameter … .vcf”**  
  Put the genome ID before the VCF or use STDIN. Also consider `-i vcf`.

- **`Cannot read .../genes.gbk`**  
  File must be named **`genes.gbk`** (not `genes.gb`). Use `-format gbwithparts` when fetching.

- **No annotations / many `INTERGENIC` calls**  
  Your VCF `CHROM` names must match **`NC_001925.1`**, **`NC_001926.1`**, **`NC_001927.1`**. Remap with:
  ```bash
  # map.txt lines: old_name <TAB> NC_001925.1   (etc.)
  bcftools annotate --rename-chrs map.txt -o fixed.vcf -O v input.vcf
  ```

- **Permission denied editing SnpEff config**  
  Don’t edit the conda share dir. Use a writable copy in `~/snpeff` and pass `-c` & `-dataDir`.

- **Apple Silicon vs Intel**  
  Commands are the same. If mixing tools that lack native arm64 builds, you can set `CONDA_SUBDIR=osx-64` before creating a separate env (rarely needed for SnpEff/EDirect).

- **Java memory (unlikely for tiny viral DBs)**  
  You can set `JAVA_TOOL_OPTIONS="-Xmx2g"` if ever needed.

---

## Verification checklist

```bash
# DB registered?
grep '^bunyamwera_ref\.genome' "$SNPEFF_HOME/snpEff.config"

# Predictor exists?
ls -lh "$SNPEFF_HOME/data/bunyamwera_ref/snpEffectPredictor.bin"

# Minimal test variant (on L segment)
printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nNC_001925.1\t100\t.\tA\tG\t.\t.\t.\n" > test.vcf
snpEff ann -c "$SNPEFF_HOME/snpEff.config" -dataDir "$SNPEFF_HOME/data" bunyamwera_ref < test.vcf > test.ann.vcf
```

---

## Use with BunyaGen

Point your BunyaGen annotation step to `bunyamwera_ref` and the config/data under `~/snpeff`:

```bash
snpEff ann \
  -c "$SNPEFF_HOME/snpEff.config" -dataDir "$SNPEFF_HOME/data" \
  bunyamwera_ref < your_sample.vcf > your_sample.snpeff.vcf
```


---
## Troubleshooting

| Symptom | Fix |
|---|---|
| `Missing dependency` | `mamba install -c bioconda -c conda-forge <tool>` |
| IQ-TREE crashes/segfault on macOS | Script auto-retries with `-T 1`, then `GTR+G`, then **FastTree -boot** |
| iTOL shows no bootstrap | Make sure you uploaded `concat.treefile` & turned on **Internal node labels** |
| Alignment invariant | Check `DIAG: concat.aln.fa` → if `varSites=0`, there’s no phylogenetic signal |
| Wrong segment IDs | Always pass `--segments MF926354.1,KP795084.1,KP795093.1` to pin L/M/S |
| Consensus empty | Script auto-repairs from reference; see `logs/*.cons.log` |

---

## Reproducibility notes

- Use **absolute paths** in `vcfs.txt`.  
- Pin segment accessions with `--segments` to avoid auto-detect mistakes.  
- For deterministic runs, keep the same tool versions (export your `conda env export > env.lock.yml`).

---

