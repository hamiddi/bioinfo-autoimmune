"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
Program: UniProt to CDS FASTA Converter

Description:
This script retrieves the coding DNA sequences (CDS) associated with a given UniProt protein ID.
It performs the following steps:
1. Retrieves the gene name associated with the UniProt ID using the UniProt REST API.
2. Searches NCBI for the corresponding RefSeq nucleotide accession (NM_ format).
3. Extracts the CDS features from the RefSeq record.
4. Saves the CDS sequences to a FASTA file named after the gene.

Requirements:
- Biopython
- requests

Setup:
- Replace `Entrez.email` with your email address to comply with NCBI usage policies.

Usage:
- Set the `uniprot_id` variable near the end of the script to the desired UniProt ID.
- Run the script. If successful, a FASTA file (e.g., HLA-A.fasta) will be created in the working directory.

Example:
    uniprot_id = "P04439"

Output:
    A FASTA file containing one or more CDS sequences with a defline format:
    >UniProtID|GeneName|RefSeqID|CDS_index
"""

import requests
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Set your email for NCBI access
Entrez.email = "your_email@example.com"

def uniprot_to_gene(uniprot_id):
    """
    Retrieves the gene name from UniProt API.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Error retrieving data for {uniprot_id}")
        return None
    
    data = response.json()
    
    # Extract gene name
    gene_name = None
    if "genes" in data and data["genes"]:
        gene_name = data["genes"][0].get("geneName", {}).get("value")
    
    return gene_name

def gene_to_refseq_nucleotide(gene_name, organism="Homo sapiens"):
    """
    Retrieves the RefSeq Nucleotide sequence (NM_ accession) for a given gene name from NCBI.
    """
    try:
        search_query = f"{gene_name}[Gene Name] AND {organism}[Organism] AND refseq[filter]"
        handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            print(f"No RefSeq Nucleotide sequence found for {gene_name} in {organism}")
            return None
        
        # Fetch the first NM_ accession found
        for nucleotide_id in record["IdList"]:
            handle = Entrez.efetch(db="nucleotide", id=nucleotide_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            # Check if the accession starts with NM_
            if record.id.startswith("NM_"):
                return record.id
        
    except Exception as e:
        print(f"Error retrieving RefSeq Nucleotide ID for {gene_name}: {e}")
    
    return None

def get_cds_from_refseq(refseq_id):
    """
    Retrieves the CDS (Coding DNA Sequence) for a given RefSeq Nucleotide ID from NCBI.
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        cds_list = []
        for feature in record.features:
            if feature.type == "CDS":
                cds_seq = feature.extract(record.seq)
                cds_list.append(str(cds_seq))
        
        return cds_list if cds_list else None

    except Exception as e:
        print(f"Error retrieving CDS for {refseq_id}: {e}")
        return None

def write_fasta(uniprot_id, gene_name, refseq_id, cds_sequences):
    """
    Writes the CDS sequences to a FASTA file with a proper defline.
    The file name is set to `gene_name.fasta`.
    """
    output_file = f"{gene_name}.fasta"  # Set the file name based on gene name

    records = []
    for i, cds_seq in enumerate(cds_sequences, 1):
        header = f">{uniprot_id}|{gene_name}|{refseq_id}|CDS_{i}"
        record = SeqRecord(Seq(cds_seq), id=header, description="")
        records.append(record)
    
    with open(output_file, "w") as fasta_file:
        SeqIO.write(records, fasta_file, "fasta")

    print(f"FASTA file '{output_file}' created successfully.")

# Example usage
uniprot_id = "P04439"  # Example UniProt ID for HLA-A
gene_name = uniprot_to_gene(uniprot_id)

if gene_name:
    print(f"Gene Name: {gene_name}")
    
    # Step 1: Get RefSeq Nucleotide ID (NM_ accession)
    refseq_nucleotide_id = gene_to_refseq_nucleotide(gene_name)

    if refseq_nucleotide_id:
        print(f"Mapped RefSeq Nucleotide ID: {refseq_nucleotide_id}")

        # Step 2: Retrieve CDS
        cds_sequences = get_cds_from_refseq(refseq_nucleotide_id)
        if cds_sequences:
            # Step 3: Write to FASTA with gene name as file name
            write_fasta(uniprot_id, gene_name, refseq_nucleotide_id, cds_sequences)
    else:
        print("No RefSeq Nucleotide ID found.")
else:
    print("No gene name found.")

