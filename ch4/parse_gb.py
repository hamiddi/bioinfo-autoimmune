# Author: Hamid D. Ismail, Ph.D.
# Book: Bioinformatics of Autoimmune Diseases
# This script parses a GenBank file using Biopython, extracts the sequence ID and sequence,
# and prints the sequence along with its length.
from Bio import SeqIO
# Parse the GenBank file
genbank_file = "example.gb"
with open(genbank_file, "r") as handle:
    record = SeqIO.read(handle, "genbank")
# Print sequence
id =  record.id
seq = record.seq
print("Sequence ID:", id)
print("Sequence Length:", len(seq))
print("Sequence:\n", seq)
