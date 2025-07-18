# Author: Hamid D. Ismail, Ph.D.
# Book: Bioinformatics of Autoimmune Diseases
# This script parses a FASTA file using Biopython, extracts each sequence and its identifier, 
# and prints them to the console.

from Bio import SeqIO
fasta_file = "example.fasta"
# Parse the FASTA file and print sequences
for record in SeqIO.parse(fasta_file, "fasta"):
    id  = record.id
    seq = record.seq
    print(f"Sequence ID: {id}")
    print(f"Sequence: {seq}")

