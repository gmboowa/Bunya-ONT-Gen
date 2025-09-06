# Bunya-ONT-Gen: A pipeline for viral ONT metagenomic analysis

**Bunya-ONT-Gen** is a modular, end-to-end analysis pipeline designed for the characterization of viral genomes from Oxford Nanopore long-read metagenomic data. It is tailored for emerging & re-emerging viruses such as Bunyaviruses & supports everything from raw read quality control to phylogenetic tree construction.

---

## Key Features

- Optimized for Oxford Nanopore (ONT) long-read data  
- Human host read removal  
- Quality control and filtering of metagenomic reads  
- De novo assembly of viral genomes  
- Taxonomic classification using multiple tools  
- Functional gene annotation  
- Phylogenetic tree for evolutionary inference  

---

## Installation

### Environment Setup for Genomic Analysis

```bash
conda env create -f Bunya-ONT-Gen.yml
conda activate Bunya-ONT-Gen
conda install -y -c conda-forge -c bioconda \
  pandas biopython openjdk \
  prodigal diamond fastqc nanoplot minimap2 samtools bcftools multiqc \
  spades kraken2 mafft fasttree seqtk seqkit flye krona snpeff medaka

CONDA_SUBDIR=osx-64 mamba create -n medaka_v1 -c bioconda -c conda-forge medaka=1.7.3

```

### Hostile: Clean Human Reads Pipeline

This pipeline removes human host reads (e.g., T2T-CHM13v2.0) from single-end ONT data using `hostile` & `minimap2`. It supports automatic reference downloading via NCBI & is built for reproducibility & ease of use.

```bash
# Install all dependencies in a clean conda environment
conda create -n hostile -c bioconda ncbi-datasets-cli minimap2 hostile
conda activate hostile
```

**Required arguments:**

| Argument      | Description |
|---------------|-------------|
| `--fasta`     | NCBI accession (e.g., GCF_009914755.4) or path to a local FASTA reference |
| `-i, --input` | Text file listing paths to FASTQ files (one per line) |
| `-o, --output`| Directory where cleaned results will be stored |

**Basic usage:**

```bash
# Example using T2T-CHM13v2.0 human reference
python3 hostile_clean_ont_human_minimap2.py --fasta GCF_009914755.4 -i sample.list.txt -o Clean
```

**Input format:**  
The input list file (`Bunya_ont_sample.list.txt`) should contain absolute paths to single-end FASTQ files, one per line:

```bash
~/sample1.fastq.gz
~/sample2.fastq.gz
```

---

## Usage

### Run Whole-Genome Analysis

```bash
source /opt/homebrew/Caskroom/miniforge/base/bin/activate

conda activate BunyaGen

```

```bash
python Bunya-ONT-Gen.py -inputs fastq_sample.txt -reference Bunya.reference.fasta 
```

**Example `sample_list.txt`:**

```bash
~/Sample_01_ONT_reads.fastq.gz
~/Sample_02_ONT_reads.fastq.gz
~/Sample_03_ONT_reads.fastq.gz
```

---

### Environment Setup for Phylogenetic Analysis

```bash
conda env create -f Bunya_phylogeny.yml
conda activate Bunya_phylogeny
```

### Run Phylogenetic Inference

```bash
Bunya_phylogeny.py -r multifasta_Bunyavirus.fasta -i fasta_sample.txt -o results -t 8 -b 100
```

**Example `fasta_sample.txt`:**

```bash
~/Sample_01_assembly.fasta
~/Sample_02_assembly.fasta
~/Sample_03_assembly.fasta
```

---

## Steps

1. **Input**
   
   - FASTQ files from ONT sequencers

3. **Human Reads Removal**
   
   - **Hostile**: Removes host sequences from short & long read (meta)genomes

5. **Quality Control**  
   - **Fastp**: Adapter trimming & quality filtering  
   - **NanoStat**: Summary statistics of read quality  
   - **NanoPlot**: Visualization of read metrics  
   - **Filtlong**: Length-based read filtering  

6. **Assembly**  
   - **Flye**, **Canu**, or **Raven** for long-read genome assembly  

7. **Classification & Annotation**  
   - **Kraken2**, **Kaiju**, **Centrifuge**, or **Bracken** for taxonomic classification  
   - **Medaka** and **Racon** for polishing ONT assemblies  

8. **Functional Annotation**  
   - **Prodigal** for gene prediction  
   - **DIAMOND** for fast similarity search against protein databases  

9. **Phylogenetic Analysis**  
   - **MAFFT** for multiple sequence alignment  
   - **IQ-TREE** for phylogenetic tree inference with model selection  

---

## Inputs & Outputs

### Inputs for Genomic Analysis

| File Type        | Description           |
|------------------|-----------------------|
| `*.fastq`        | ONT sequence files    |
| `reference.fasta`| Bunya reference genome|

### Inputs for Phylogenetic Analysis

| Parameter        | Description                       | Default   |
|------------------|-----------------------------------|-----------|
| `-inputs`        | File with FASTQ paths             | Required  |
| `-reference`     | Bunya reference genome            | Required  |
| `-threads`       | CPU threads for parallel steps    | `8`       |
| `-bootstrap`     | Phylogenetic bootstrap replicates | `100`    |
| `-min_coverage`  | Consensus calling threshold       | `20x`     |
| `-tree_model`    | IQ-TREE substitution model        | `MFP`     |
| `-aln_consensus` | trimAl conservation threshold     | `60%`     |

---

### Output Structure

```bash
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

| Issue             | Recommendation |
|-------------------|----------------|
| Memory errors     | Reduce thread count (e.g., `-threads 4`) |
| Assembly failures | Inspect quality reports in `results/qc/nanoplot/` |
| Dependency issues | Update Conda with `conda env update -f Bunya-ONT-Gen.yml` |

---
