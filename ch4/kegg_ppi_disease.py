"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
Protein-Protein Interaction (PPI) Retrieval for Autoimmune Disease Genes using KEGG and STRING APIs

Description:
This script retrieves protein-protein interaction (PPI) data for genes associated with KEGG pathways 
linked to a specified autoimmune disease. It uses the KEGG REST API to find disease-associated pathways 
and genes, converts gene IDs to STRING-compatible identifiers, and then queries the STRING database 
for known PPIs. Results are saved in both CSV and JSON formats.

Usage:
1. Set the autoimmune disease name in the `disease_name` variable at the bottom of the script.
2. Run the script using Python 3:
   $ python script_name.py
3. The script will:
   - Retrieve KEGG pathways related to the disease.
   - Extract associated human genes.
   - Convert KEGG gene IDs to STRING-compatible IDs.
   - Query STRING for protein-protein interactions.
   - Save the resulting PPI data to `ppi_results.csv` and `ppi_results.json`.

Dependencies:
- requests
- biopython
- xml (ElementTree)
- json
- csv

Note:
Ensure you have an internet connection as this script relies on external web APIs (KEGG and STRING).
"""
import requests
from Bio import Entrez
import xml.etree.ElementTree as ET
import json
import csv

def get_kegg_pathways(disease_name):
    """Retrieve KEGG pathways associated with a given autoimmune disease."""
    base_url = "https://rest.kegg.jp/find/disease/"
    response = requests.get(base_url + disease_name)
    
    if response.status_code != 200:
        print("Error retrieving data from KEGG.")
        return []
    
    pathways = []
    for line in response.text.strip().split('\n'):
        fields = line.split('\t')
        if len(fields) > 1:
            disease_id = fields[0].split(':')[1]
            pathways.append(disease_id)
    return pathways

def get_kegg_pathway_genes(pathway_id):
    """Retrieve genes associated with a KEGG pathway."""
    base_url = f"https://rest.kegg.jp/link/hsa/{pathway_id}"
    response = requests.get(base_url)
    
    if response.status_code != 200:
        print(f"Error retrieving genes for pathway {pathway_id}.")
        return []
    
    genes = set()
    for line in response.text.strip().split('\n'):
        fields = line.split('\t')
        if len(fields) > 1:
            gene_id = fields[1].split(':')[1]
            genes.add(gene_id)
    return list(genes)

def convert_kegg_to_string_ids(genes):
    """Convert KEGG gene IDs to STRING-compatible IDs using KEGG and STRING API."""
    ncbi_gene_ids = []
    
    # Convert KEGG gene IDs to NCBI Gene IDs
    for gene in genes:
        response = requests.get(f"https://rest.kegg.jp/conv/ncbi-geneid/hsa:{gene}")
        if response.status_code == 200 and response.text.strip():
            for line in response.text.splitlines():
                fields = line.split("\t")
                if len(fields) == 2:
                    ncbi_id = fields[1].split(":")[1]  # Extract NCBI Gene ID
                    ncbi_gene_ids.append(ncbi_id)
    
    if not ncbi_gene_ids:
        print("Failed to retrieve NCBI Gene IDs.")
        return []
    
    print(f"Converted KEGG IDs to NCBI Gene IDs: {ncbi_gene_ids}")  # Debugging
    
    # Convert NCBI Gene IDs to STRING IDs using STRING API
    string_ids = []
    string_api_url = "https://string-db.org/api/json/get_string_ids?identifiers=" + "%0d".join(ncbi_gene_ids) + "&species=9606"
    response = requests.get(string_api_url)
    
    if response.status_code == 200 and response.json():
        for entry in response.json():
            if "stringId" in entry:
                string_ids.append(entry["stringId"])
    
    if not string_ids:
        print("Failed to retrieve STRING-compatible IDs.")
        return []
    
    print(f"Converted NCBI Gene IDs to STRING-compatible IDs: {string_ids}")  # Debugging
    return string_ids

def get_ppi_from_stringdb(genes):
    """Retrieve protein-protein interactions from STRING database with debugging."""
    if not genes:
        print("No genes available for PPI search.")
        return []
    
    genes = [gene.strip() for gene in genes if gene.strip()]
    
    if not genes:
        print("No valid STRING-compatible IDs found after cleaning.")
        return []
    
    string_api_url = "https://string-db.org/api/xml/network?species=9606&identifiers="
    query = "%0d".join(genes)  # STRING API expects %0d as separator for multiple IDs
    
    print(f"Querying STRING API with STRING-compatible IDs: {query}")  # Debugging
    response = requests.get(string_api_url + query)
    
    if response.status_code != 200:
        print("Error retrieving PPI data from STRING DB.")
        print(f"Response Code: {response.status_code}")
        print(f"Response Text: {response.text}")  # Show error message from STRING
        return []
    
    ppis = []
    root = ET.fromstring(response.text)
    for interaction in root.findall(".//record"):
        protein_a = interaction.find(".//preferredName_A").text
        protein_b = interaction.find(".//preferredName_B").text
        ppis.append((protein_a, protein_b))
    
    return ppis

def save_ppi_to_csv(ppi_data, filename="ppi_results.csv"):
    """Save PPI data to a CSV file."""
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Protein A", "Protein B"])
        writer.writerows(ppi_data)
    print(f"PPI data saved to {filename}")

def save_ppi_to_json(ppi_data, filename="ppi_results.json"):
    """Save PPI data to a JSON file."""
    with open(filename, "w") as file:
        json.dump(ppi_data, file, indent=4)
    print(f"PPI data saved to {filename}")

if __name__ == "__main__":
    disease_name = "rheumatoid arthritis"  # Example autoimmune disease
    pathways = get_kegg_pathways(disease_name)
    
    if not pathways:
        print("No pathways found.")
    else:
        print(f"Found {len(pathways)} pathways.")
        all_genes = set()
        for pathway in pathways:
            genes = get_kegg_pathway_genes(pathway)
            all_genes.update(genes)
        
        print(f"Total unique genes: {len(all_genes)}")
        if all_genes:
            converted_genes = convert_kegg_to_string_ids(list(all_genes))
            ppi_data = get_ppi_from_stringdb(converted_genes)
            print(f"Retrieved {len(ppi_data)} protein-protein interactions.")
            for ppi in ppi_data[:10]:  # Display first 10 interactions
                print(ppi)
            
            # Save results
            save_ppi_to_csv(ppi_data)
            save_ppi_to_json(ppi_data)

