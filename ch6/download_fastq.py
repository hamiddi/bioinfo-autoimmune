"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autoimmune Diseases
############################################################

This script downloads FASTQ files from NCBI's SRA using `fasterq-dump`.

Usage:
1. Prepare a plain text file named `ids.txt` in the same directory as this script.
   - Each line of the file should contain one SRA Run ID (e.g., SRR12345678).
2. Run the script: `python3 download_fastqs.py`
   - The script will create a directory called `fastqs` (if it does not already exist).
   - It will then download the corresponding FASTQ files into this directory.
   - If the `fastqs` directory already exists, the script will exit with a warning.

Requirements:
- `fasterq-dump` must be installed and accessible in your system's PATH.
- Ensure you have sufficient disk space and internet bandwidth.
- The script uses 18 threads for downloading.

"""
import os
import subprocess
import sys

# Define input file and output directory
run_ids_file = "ids.txt"
output_dir = "data/raw"

# Check if fastqs directory exists
if os.path.exists(output_dir):
    print(f"Warning: Directory '{output_dir}' already exists. Exiting.")
    sys.exit(1)
else:
    os.makedirs(output_dir)
    print(f"Created directory '{output_dir}'.")

# Read run IDs from the file
try:
    with open(run_ids_file, 'r') as f:
        run_ids = [line.strip() for line in f if line.strip()]
except FileNotFoundError:
    print(f"Error: File '{run_ids_file}' not found.")
    sys.exit(1)
# Download FASTQ files using fasterq-dump
for run_id in run_ids:
    print(f"Downloading {run_id}...")
    try:
        subprocess.run([
            "fasterq-dump", run_id,
            "-O", output_dir,
            "--skip-technical",
            "--split-files",
            "--progress",
            "--threads", "18"
        ], check=True)
        print(f"Finished downloading {run_id}")
    except subprocess.CalledProcessError as e:
        print(f"Error downloading {run_id}: {e}")

print("All downloads completed.")
