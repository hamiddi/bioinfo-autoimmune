#!/bin/bash
# Check if input file is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 sra_ids.txt"
    exit 1
fi
input_file="$1"
# Check if fastq-dump is available
if ! command -v fastq-dump &> /dev/null; then
    echo "Error: fastq-dump is not installed or not in PATH."
    exit 1
fi
# Download paired-end FASTQ files
while IFS= read -r sra_id; do
    if [ -n "$sra_id" ]; then
        echo "Downloading $sra_id..."
        fastq-dump --split-files --gzip "$sra_id"
    fi
done < "$input_file"
echo "Download completed."
