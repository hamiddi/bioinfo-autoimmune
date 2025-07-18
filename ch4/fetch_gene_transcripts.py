"""
Author: Hamid D. Ismail, Ph.D.
Book title: Bioinformatics of Autimmune Diseases
##########################################################
This script retrieves RefSeq transcript and protein sequences for a given NCBI Gene ID 
using NCBI's E-Utilities (via Biopython's Entrez module), and saves them as FASTA files.

How to use:
1. Replace the 'gene_id' variable in the main() function with the NCBI Gene ID of interest.
2. Set your email address in the 'Entrez.email' field (required by NCBI for API usage).
3. Run the script. It will:
   - Retrieve RefSeq transcript IDs and protein IDs associated with the gene.
   - Download the corresponding nucleotide and protein FASTA sequences.
   - Save them in files named 'gene_<gene_id>_transcripts.fasta' and 'gene_<gene_id>_proteins.fasta'.

Dependencies:
- Biopython (install via: pip install biopython)

Example:
If gene_id = "3105" (ESR1 gene), it will generate:
- gene_3105_transcripts.fasta
- gene_3105_proteins.fasta
"""
import os
from Bio import Entrez, SeqIO
Entrez.email = "your_email@example.com"
def get_refseq_ids(gene_id):
    transcript_ids, protein_ids = [], []
    # Query NCBI Gene Database
    handle = Entrez.elink(dbfrom="gene", 
                          db="nucleotide",
                          id=gene_id, 
                          linkname="gene_nuccore_refseqrna")
    records = Entrez.read(handle)
    handle.close()
    # Extract transcript RefSeq IDs
    if records and records[0]['LinkSetDb']:
        linkset_db_tr = records[0]['LinkSetDb'][0]['Link']
        transcript_ids = [
                          link["Id"]
                          for link in linkset_db_tr
                         ] 
    # Query for protein IDs
    handle = Entrez.elink(dbfrom="gene", db="protein", id=gene_id,
                          linkname="gene_protein_refseq")
    records = Entrez.read(handle)
    handle.close()
    # Extract protein RefSeq IDs
    if records and records[0]['LinkSetDb']:
        linkset_db_pr = records[0]['LinkSetDb'][0]['Link']
        protein_ids = [
                       link["Id"]
                       for link in linkset_db_pr
                      ]
    return transcript_ids, protein_ids
def fetch_fasta(seq_ids, db, output_file):
    if not seq_ids:
        print(f"No sequences found in {db} database.")
        return
    handle = Entrez.efetch(db=db, 
                          id=",".join(seq_ids), 
                          rettype="fasta", 
                          retmode="text")
    fasta_data = handle.read()
    handle.close()
    with open(output_file, "w") as f:
        f.write(fasta_data)
    print(f"Saved {db} sequences to {output_file}")
def main(gene_id):
    transcript_ids, protein_ids = get_refseq_ids(gene_id)
    # Save transcripts
    fetch_fasta(transcript_ids, "nucleotide", f"gene_{gene_id}_transcripts.fasta")
    # Save proteins
    fetch_fasta(protein_ids, "protein", f"gene_{gene_id}_proteins.fasta")
if __name__ == "__main__":
    gene_id = "3105"
    main(gene_id)

