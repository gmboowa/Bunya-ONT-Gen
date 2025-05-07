
# BunyaGen: A Pipeline for Viral ONT Metagenomic Analysis

**BunyaGen** is a modular, end-to-end analysis pipeline designed for the characterization of viral genomes from Oxford Nanopore long-read metagenomic data. It is tailored for emerging and re-emerging viruses such as Bunyaviruses and supports everything from raw read quality control to phylogenetic tree construction.

---

## Key Features

- Optimized for Oxford Nanopore (ONT) long-read data
- Quality control and filtering of metagenomic reads
- De novo assembly of viral genomes
- Taxonomic classification using multiple tools
- Functional gene annotation
- Phylogenetic analysis pipeline for evolutionary inference
- Flexible tool selection for each stage
- Easy customization for additional pathogens or workflows

---

## Workflow Overview

![BunyaGen Workflow](BunyaGen_Workflow.png)

> The diagram illustrates a comprehensive ONT metagenomic analysis pipeline, ending with phylogenetic tree construction using MAFFT and IQ-TREE.

### Detailed Steps:

#### 1. **Input**
- FASTQ files from ONT sequencers

#### 2. **Quality Control**
- **Fastp**: Adapter trimming and quality filtering
- **NanoStat**: Summary statistics of read quality
- **NanoPlot**: Visualization of read metrics
- **Filtlong**: Length-based read filtering

#### 3. **Assembly**
- Tools such as **Flye**, **Canu**, or **Raven** for long-read genome assembly

#### 4. **Classification & Annotation**
- **Kraken2**, **Kaiju**, **Centrifuge**, or **Bracken** for taxonomic classification
- **Medaka** and **Racon** for polishing ONT assemblies

#### 5. **Functional Annotation**
- **Prodigal** for gene prediction
- **DIAMOND** for fast similarity search against protein databases

#### 6. **Phylogenetic Analysis**
- **MAFFT** for multiple sequence alignment
- **IQ-TREE** for phylogenetic tree inference with model selection

---

## Installation

You can install all dependencies via conda using the provided `installer.yml`:

```bash
conda env create -f installer.yml
conda activate bunyagen
```

---

## Usage

Run the pipeline with the following structure:

```bash
bash run_bunyagen.sh \
  -i path/to/reads.fastq \
  -o output_directory \
  --min_length 1000 \
  --threads 8
```

Refer to `docs/config.yaml` to customize tool parameters and file paths.

---

## Output

- QC reports (`.html`, `.json`, `.pdf`)
- Assembled genome(s) (`.fasta`)
- Taxonomic classification table
- Predicted genes and protein functions
- Multiple sequence alignment and `.nwk` tree file
- Publication-ready phylogenetic tree visualizations

---

## Citing BunyaGen

If you use this pipeline, please cite this repository and the relevant tools individually.

---

## License

MIT License. See [LICENSE](LICENSE) for more information.

---

## Maintainer

**Gerald Mboowa**  
Bioinformatics Implementation Expert  
Broad Institute | Africa CDC  
Email: [your_email@example.com](mailto:your_email@example.com)
