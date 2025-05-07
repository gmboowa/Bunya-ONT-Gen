#!/bin/bash

# BunyaGen pipeline runner
# Example: bash run_bunyagen.sh -i input.fastq -o output_dir --min_length 1000 --threads 8

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) input="$2"; shift ;;
        -o|--output) output="$2"; shift ;;
        --min_length) min_length="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

mkdir -p "$output/qc"

echo "Running Fastp..."
fastp -i "$input" -o "$output/qc/filtered.fastq" --length_required "$min_length" -w "$threads"

echo "Running NanoStat..."
NanoStat --fastq "$output/qc/filtered.fastq" -o "$output/qc/nanostat.txt"

echo "Running NanoPlot..."
NanoPlot --fastq "$output/qc/filtered.fastq" -o "$output/qc/nanoplot" --threads "$threads"

echo "Running Filtlong..."
filtlong --min_length "$min_length" "$output/qc/filtered.fastq" > "$output/qc/filtlong.fastq"

# Add additional steps for assembly, annotation, and phylogenetics as needed
echo "Pipeline step placeholders complete. Extend as needed."
