#!/bin/bash
#Author: Hamid D. Ismail, Ph.D.
#Book title: Bioinformatics of Autoimmune Diseases
# Check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 sra_ids.txt"
    exit 1
fi
input_file="$1"
# Check if fasterq-dump is available
if ! command -v fasterq-dump &> /dev/null; then
    echo "Error: fasterq-dump is not installed or not in PATH."
    exit 1
fi
# Set number of threads
threads=8
# Download paired-end FASTQ files using fasterq-dump
while IFS= read -r sra_id; do
    if [ -n "$sra_id" ]; then
        echo "Downloading $sra_id with $threads threads..."
        fasterq-dump --split-files --threads "$threads" --verbose "$sra_id"
        gzip "${sra_id}_1.fastq" "${sra_id}_2.fastq"
    fi
done < "$input_file"

echo "Download completed."
