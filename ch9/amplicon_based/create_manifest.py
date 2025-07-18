"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This Python script generates a QIIME 2-compatible manifest file for paired-end FASTQ sequencing data. 
It scans the `data/raw/` directory for files ending in `_1.fastq.gz` and `_2.fastq.gz`, associates 
each pair by sample ID, and writes a CSV manifest that includes the absolute paths to the forward 
and reverse read files.

This manifest file is commonly used to import data into QIIME 2 or other downstream bioinformatics pipelines 
that require explicit pairing of sequencing reads.

Required Python Packages:
- os
- glob

Input:
- Paired-end FASTQ files in `data/raw/`, named as `{sample_id}_1.fastq.gz` and `{sample_id}_2.fastq.gz`

Output:
- Manifest file: `meta/manifest.csv` with columns:
  - sample-id
  - forward-absolute-filepath
  - reverse-absolute-filepath

Example Output Row:
sample123,/absolute/path/to/sample123_1.fastq.gz,/absolute/path/to/sample123_2.fastq.gz

Notes:
- Ensure that all paired files follow the expected naming convention.
- This script assumes all FASTQ files are compressed with `.gz` and are located in the `data/raw/` directory.
"""



import os
import glob

raw_dir = os.path.abspath("data/raw")
manifest_file = "meta/manifest.csv"

samples = {}

for f in glob.glob(os.path.join(raw_dir, "*.fastq.gz")):
    basename = os.path.basename(f)
    if "_1.fastq.gz" in basename:
        sample_id = basename.replace("_1.fastq.gz", "")
        samples.setdefault(sample_id, {})["forward"] = f
    elif "_2.fastq.gz" in basename:
        sample_id = basename.replace("_2.fastq.gz", "")
        samples.setdefault(sample_id, {})["reverse"] = f

with open(manifest_file, "w") as out:
    out.write("sample-id,forward-absolute-filepath,reverse-absolute-filepath\n")
    for sample, paths in samples.items():
        out.write(f"{sample},{paths['forward']},{paths['reverse']}\n")

print(f"Manifest saved to: {manifest_file}")

