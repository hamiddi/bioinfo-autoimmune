"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script renames FASTQ files that follow the Illumina CASAVA 1.8+ naming convention 
(e.g., `sample1_S1_L001_R1_001.fastq.gz`) into a simpler, cleaner format 
(e.g., `sample1_1.fastq.gz` and `sample1_2.fastq.gz`). This renaming helps standardize filenames 
for easier downstream processing in bioinformatics pipelines.

Required Python Packages:
- os
- re

Input:
- FASTQ files located in a target directory (default: `data/raw`)
- Expected filename pattern: `{sampleID}_S##_L###_R[12]_001.fastq.gz`

Output:
- Renamed FASTQ files using the format `{sampleID}_1.fastq.gz` and `{sampleID}_2.fastq.gz`
- Original files are renamed in place

Usage:
Run the script as a standalone utility to rename all matching FASTQ files:
```bash
python rename_fastq_files.py

Note:
Update the target_directory variable if your FASTQ files are in a different location.
This script assumes all files are single-lane, paired-end reads with CASAVA formatting.
"""



import os
import re

def rename_fastq_files(directory):
    casava_pattern = re.compile(r"(.+)_S\d+_L\d{3}_R([12])_001\.fastq\.gz")

    for filename in os.listdir(directory):
        if filename.endswith(".fastq.gz"):
            match = casava_pattern.match(filename)
            if match:
                sample_id = match.group(1)
                read_direction = match.group(2)
                new_name = f"{sample_id}_{read_direction}.fastq.gz"
                old_path = os.path.join(directory, filename)
                new_path = os.path.join(directory, new_name)
                os.rename(old_path, new_path)
                print(f"Renamed: {filename} ➜ {new_name}")

if __name__ == "__main__":
    target_directory = "data/raw"  # change this path if needed
    rename_fastq_files(target_directory)

