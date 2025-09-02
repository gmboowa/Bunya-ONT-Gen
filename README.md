
# Bunya-ONT-Gen: A pipeline for viral ONT metagenomic analysis

**Bunya-ONT-Gen** is a modular, end-to-end analysis pipeline designed for the characterization of viral genomes from Oxford Nanopore long-read metagenomic data. It is tailored for emerging and re-emerging viruses such as Bunyaviruses and supports everything from raw read quality control to phylogenetic tree construction.

---

## Key features

- Optimized for Oxford Nanopore (ONT) long-read data
- Human host read removal
- Quality control and filtering of metagenomic reads
- De novo assembly of viral genomes
- Taxonomic classification using multiple tools
- Functional gene annotation
- Phylogenetic tree for evolutionary inference


---

## Workflow overview

![BunyaGen Workflow](BunyaGen_Workflow.png)

> The diagram illustrates nanopore metagenomic analysis pipeline, ending with phylogenetic tree construction using MAFFT and IQ-TREE.

### Steps:


#### 1. **Input**
- FASTQ files from ONT sequencers

#### 2. **Human reads removal**
- **Hostile**: Removes host sequences from short and long read (meta)genomes,
  
#### 3. **Quality control**
- **Fastp**: Adapter trimming and quality filtering
- **NanoStat**: Summary statistics of read quality
- **NanoPlot**: Visualization of read metrics
- **Filtlong**: Length-based read filtering

#### 4. **Assembly**
- Tools such as **Flye**, **Canu**, or **Raven** for long-read genome assembly

#### 5. **Classification & annotation**
- **Kraken2**, **Kaiju**, **Centrifuge**, or **Bracken** for taxonomic classification
- **Medaka** and **Racon** for polishing ONT assemblies

#### 6. **Functional annotation**
- **Prodigal** for gene prediction
- **DIAMOND** for fast similarity search against protein databases

#### 7. **Phylogenetic analysis**
- **MAFFT** for multiple sequence alignment
- **IQ-TREE** for phylogenetic tree inference with model selection

---

## Installation

# Hostile clean human reads pipeline

This pipeline removes human host reads (e.g., T2T-CHM13v2.0) from single-end Oxford Nanopore sequencing (ONT) data using `hostile` and `minimap2`. It supports automatic reference downloading via NCBI and is built for reproducibility and ease of use.

---


```
# Install all dependencies in a clean conda environment:
conda create -n hostile -c bioconda ncbi-datasets-cli minimap2 hostile
conda activate hostile



Required arguments

Argument        Description
--fasta	         NCBI accession (e.g., GCF_009914755.4) or path to a local FASTA reference
-i, --input      Text file listing paths to FASTQ files (one per line)
-o, --output     Directory where cleaned results will be stored


```
Basic usage

```

# General usage
python3 hostile_clean_ont_human_minimap2.py \
  --fasta  <ACCESSION_OR_PATH> \
  -i       <INPUT_LIST> \
  -o       <OUTPUT_DIR>

# Example using T2T-CHM13v2.0 human reference
python3 hostile_clean_ont_human_minimap2.py \
  --fasta  GCF_009914755.4 \
  -i       sample.list.txt \
  -o       Clean

```
```

python3 hostile_clean_ont_human_minimap2.py --fasta GCF_009914755.4 -i sample.list.txt -o Clean

```

```

Input format

The input list file (Bunya_ont_sample.list.txt) should contain absolute or relative paths to single-end FASTQ files, one per line:


/path/to/sample1.fastq.gz
/path/to/sample2.fastq.gz

```

```

### Environment setup for genomic analysis

```
conda env create -f Bunya-ONT-Gen.yml
conda activate Bunya-ONT-Gen
conda install -c bioconda -c conda-forge fastqc nanoplot minimap2 samtools bcftools medaka \
multiqc spades kraken2 mafft fasttree seqtk flye krona snpeff -y

```

### Setup Kraken2 database

```bash
bash scripts/setup_kraken_db.sh

```

---

## Usage

### Run whole-genome analysis

```
python Bunya-ONT-Gen.py -inputs fastq_sample.txt -reference Bunya.reference.fasta -threads 8

```
### Example sample_list.txt

```
/home/user/Bunya_data/Sample_01_ONT_reads.fastq.gz
/home/user/Bunya_data/Sample_02_ONT_reads.fastq.gz
/home/user/Bunya_data/Sample_03_ONT_reads.fastq.gz

```
---

### Environment setup for phylogenetic analysis

```
conda env create -f Bunya_phylogeny.yml
conda activate Bunya_phylogeny

```

### Run phylogenetic inference

```

Bunya_phylogeny.py -r Bunyavirus.gbk -i fasta_sample.txt -o results -t 8 -b 1000

```

### Example fasta_sample.txt

```

/home/user/Bunya_assemblies/Sample_01_assembly.fasta
/home/user/Bunya_assemblies/Sample_02_assembly.fasta
/home/user/Bunya_assemblies/Sample_03_assembly.fasta

```
---

## Inputs and outputs

### Inputs for genomic analysis

```

| File Type        | Description                |
|------------------|----------------------------|
| `*.fastq`        | ONT sequence files         |
| `reference.fasta`| Bunya reference genome      |

### Inputs for Phylogenetic analysis

| Parameter         | Description                          | Default       |
|-------------------|--------------------------------------|---------------|
| `-inputs`         | File with FASTQ paths                | Required      |
| `-reference`      | Bunya reference genome               | Required      |
| `-threads`        | CPU threads for parallel steps       | `8`           |
| `-bootstrap`      | Phylogenetic bootstrap replicates    | `1000`        |
| `-min_coverage`   | Consensus calling threshold          | `20x`         |
| `-tree_model`     | IQ-TREE substitution model           | `MFP (auto)`  |
| `-aln_consensus`  | trimAl conservation threshold        | `60%`         |

---

### Output structure


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
| Assembly failures   | Inspect quality reports in `results/qc/nanoplot/`               |
| Dependency problems | Update Conda with `conda env update -f Bunya-ONT-Gen.yml`            |


---


