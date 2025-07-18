"""
Book title: Bioinformatics of Autoimmune Diseases
Author: Hamid D. Ismail, Ph.D.
############################################################
This script reads genomic region information from a CSV file and fetches corresponding nucleotide sequences
from the NCBI Nucleotide database using the Entrez API. The sequences are then saved as FASTA files in a specified
output directory.

Usage Instructions:
1. Ensure you have internet access and the Biopython library installed (`pip install biopython`).
2. Replace 'your_email@example.com' with your actual email address to comply with NCBI Entrez API usage policy.
3. Prepare a CSV file named 'extracted_data.csv' in the same directory as this script.
   The CSV must include the following columns: chraccver, chrstart, chrstop, and optionally: uid, name, mim,
   exoncount, organism, taxid.
4. Run the script using Python. It will:
   - Create an output directory named 'fasta_sequences' (if it doesn't already exist),
   - Fetch sequences for valid entries from the CSV,
   - Write each sequence to an individual FASTA file with a descriptive header.
5. FASTA files will be saved under 'fasta_sequences/' using 'uid' as the filename.

Note: Any missing or invalid entries in the CSV will be skipped with an error message.
"""

import os
import csv
from Bio import Entrez, SeqIO

# Set your email for NCBI Entrez
Entrez.email = "your_email@example.com"  # Change this to your actual email

# Define input CSV file path
input_file = "extracted_data.csv"

# Define output directory for FASTA files
output_dir = "fasta_sequences"
os.makedirs(output_dir, exist_ok=True)

def fetch_sequence(accver, start, stop):
    """Fetch sequence from NCBI using accession number and coordinates."""
    try:
        handle = Entrez.efetch(db="nucleotide", id=accver, rettype="fasta", strand=1, seq_start=start, seq_stop=stop)
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return record
    except Exception as e:
        print(f"Error fetching {accver}:{start}-{stop}: {e}")
        return None

with open(input_file, newline='') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',')  # Adjust delimiter if necessary

    for row in reader:
        # Check if required fields are present
        chraccver = row.get("chraccver", "").strip()
        chrstart = row.get("chrstart", "").strip()
        chrstop = row.get("chrstop", "").strip()

        if chraccver and chrstart and chrstop and chraccver != "NA" and chrstart != "NA" and chrstop != "NA":
            try:
                start = int(chrstart)
                stop = int(chrstop)
            except ValueError:
                print(f"Skipping invalid coordinates: {chraccver}, {chrstart}, {chrstop}")
                continue

            # Fetch sequence
            sequence_record = fetch_sequence(chraccver, start, stop)
            if sequence_record:
                # Build FASTA defline
                defline_parts = [
                    f"uid={row.get('uid', 'NA')}",
                    f"name={row.get('name', 'NA')}",
                    f"mim={row.get('mim', 'NA')}",
                    f"chraccver={chraccver}",
                    f"chrstart={chrstart}",
                    f"chrstop={chrstop}",
                    f"exoncount={row.get('exoncount', 'NA')}",
                    f"organism={row.get('organism', 'NA')}",
                    f"taxid={row.get('taxid', 'NA')}"
                ]
                fasta_defline = ">" + "|".join(defline_parts)

                # Define output file path
                fasta_filename = f"{row.get('uid', 'unknown')}.fasta"
                fasta_filepath = os.path.join(output_dir, fasta_filename)

                # Write to FASTA file
                with open(fasta_filepath, "w") as fasta_file:
                    fasta_file.write(fasta_defline + "\n")
                    fasta_file.write(str(sequence_record.seq) + "\n")

                print(f"Saved: {fasta_filepath}")

print("Processing complete.")

