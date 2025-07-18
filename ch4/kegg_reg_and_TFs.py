"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This program searches the KEGG database for pathways associated with a given autoimmune disease,
extracts the genes involved in those pathways, and identifies potential transcription factors 
among them. The results are saved to two text files: one listing all genes and another listing 
only transcription factors.

Usage:
- By default, the script is set to analyze "rheumatoid arthritis".
- To analyze a different disease, uncomment the input() line in the main() function and comment 
  out the hardcoded `disease_name` assignment.
- Output files will be saved in the current directory with names based on the disease name:
  e.g., "rheumatoid_arthritis_genes.txt" and "rheumatoid_arthritis_tfs.txt".

Requirements:
- Install the `bioservices` Python package (e.g., using `pip install bioservices`).
"""

from bioservices import KEGG
import re

def find_disease_pathway(disease_name):
    """
    Searches KEGG for pathways related to a specific autoimmune disease.
    """
    kegg = KEGG()
    disease_data = kegg.find("pathway", disease_name)

    if not disease_data:
        return {}

    pathways = {
        line.split("\t")[0].replace("path:", "").replace("map", "hsa"): line.split("\t")[1]
        for line in disease_data.split("\n") if line
    }

    return pathways

def get_genes_from_pathway(pathway_id):
    """
    Retrieves all genes associated with a given KEGG pathway.
    """
    kegg = KEGG()
    linked_genes = kegg.link("hsa", pathway_id)

    if not linked_genes:
        return []

    genes = [
        line.split("\t")[1].replace("hsa:", "") for line in linked_genes.split("\n") if line
    ]
    return genes

def get_transcription_factors(genes):
    """
    Identifies potential transcription factors from a gene list using KEGG.
    """
    kegg = KEGG()
    transcription_factors = set()

    for gene in genes:
        gene_info = kegg.get(f"hsa:{gene}")
        if gene_info:
            if re.search(r"transcription factor|TF|STAT|FOXP|NF-kB", gene_info, re.IGNORECASE):
                transcription_factors.add(gene)

    return transcription_factors

def save_to_file(filename, data):
    """
    Saves extracted genes or transcription factors to a text file.
    """
    with open(filename, "w") as f:
        for item in data:
            f.write(item + "\n")

def main():
    #disease_name = input("Enter an autoimmune disease (e.g., rheumatoid arthritis, lupus): ").strip()
    disease_name = "rheumatoid arthritis"
    safe_disease_name = disease_name.replace(" ", "_").lower()  # Create a safe filename

    print(f"\nSearching for KEGG pathways related to {disease_name}...\n")
    pathways = find_disease_pathway(disease_name)

    if not pathways:
        print(f"No KEGG pathways found for {disease_name}.")
        return

    print(f"KEGG Pathways for {disease_name}:")
    for pathway_id, name in pathways.items():
        print(f"- {name} ({pathway_id})")

    all_genes = set()
    all_tfs = set()

    for pathway_id in pathways.keys():
        genes = get_genes_from_pathway(pathway_id)
        all_genes.update(genes)

    # Identify transcription factors from the extracted genes
    all_tfs = get_transcription_factors(all_genes)

    # Save results to files
    genes_filename = f"{safe_disease_name}_genes.txt"
    tfs_filename = f"{safe_disease_name}_tfs.txt"

    save_to_file(genes_filename, all_genes)
    save_to_file(tfs_filename, all_tfs)

    print(f"\nGenes involved in regulation for {disease_name} saved to: {genes_filename}")
    print(f"Identified Transcription Factors saved to: {tfs_filename}")

if __name__ == "__main__":
    main()

