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
  -l vcfs.txt \
  -o bunya_tree_out_concat_localref \
  --segments MF926354.1,KP795084.1,KP795093.1 \
  --refL L.fa --refM M.fa --refS S.fa \
  --boot 1000
```

**Use IUPAC ambiguity codes in the consensus (instead of `-H 1` haplotype selection):**
```bash
./build_bunya_concat_tree.sh \
  -l vcfs.txt \
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

## Troubleshooting

| Symptom | Fix |
|---|---|
| `Missing dependency` | `mamba install -c bioconda -c conda-forge <tool>` |
| IQ-TREE crashes/segfault on macOS | Script auto-retries with `-T 1`, then `GTR+G`, then **FastTree -boot** |
| iTOL shows no bootstrap | Make sure you uploaded `concat.treefile` and turned on **Internal node labels** |
| Alignment invariant | Check `DIAG: concat.aln.fa` → if `varSites=0`, there’s no phylogenetic signal |
| Wrong segment IDs | Always pass `--segments MF926354.1,KP795084.1,KP795093.1` to pin L/M/S |
| Consensus empty | Script auto-repairs from reference; see `logs/*.cons.log` |

---

## Reproducibility notes

- Use **absolute paths** in `vcfs.txt`.  
- Pin segment accessions with `--segments` to avoid auto-detect mistakes.  
- For deterministic runs, keep the same tool versions (export your `conda env export > env.lock.yml`).

---

## License

MIT

---

## Appendix — run command used in the paper/analysis

```bash
./build_bunya_concat_tree.sh \
  -l /Users/gmboowa/Bunya-ONT-Gen/vcfs.txt \
  -o bunya_tree_out_concat_real \
  --segments MF926354.1,KP795084.1,KP795093.1 \
  --auto-ref \
  --boot 1000
```
