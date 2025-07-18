"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script retrieves gene and variant information associated with a given autoimmune disease.

Usage:
- By default, the script is set to use "rheumatoid arthritis" as the disease name.
- To analyze a different disease, uncomment the input line in the `main()` function and provide the desired disease name.
- The script performs the following:
    1. Retrieves the KEGG disease entry.
    2. Extracts associated human genes.
    3. Converts KEGG gene IDs to gene symbols.
    4. Retrieves variant data from ClinVar, GWAS Catalog, and Ensembl for each gene.
    5. Saves the results to a CSV file named "genetic_variants.csv".

Requirements:
- Python 3.x
- Internet connection (uses REST APIs: KEGG, MyGene.info, NCBI ClinVar, GWAS Catalog, Ensembl)

Output:
- A CSV file containing KEGG gene IDs, gene symbols, and associated variant identifiers.
"""

import requests
import re
import csv

def get_disease_entry(disease_name):
    """Retrieve KEGG disease entry ID for a given disease name."""
    url = "http://rest.kegg.jp/find/disease/" + disease_name.replace(" ", "%20")
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.strip().split("\n")
        if lines:
            return lines[0].split("\t")[0]  # Return the first disease ID found
    return None

def get_genes_associated_with_disease(disease_id):
    """Retrieve genes associated with a given KEGG disease ID."""
    url = f"http://rest.kegg.jp/get/{disease_id}"
    response = requests.get(url)
    if response.status_code == 200:
        genes = re.findall(r"HSA:\d+", response.text)
        return list(set(genes))  # Return unique gene entries
    return []

def convert_hsa_to_gene_symbol(hsa_id):
    """Convert KEGG HSA gene ID to a gene symbol using MyGene.info API."""
    url = f"https://mygene.info/v3/gene/{hsa_id}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data.get("symbol", hsa_id)  # Return gene symbol if found
    return hsa_id

def get_gene_variants_clinvar(gene_symbol):
    """Retrieve variant information from ClinVar for a given gene symbol."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={gene_symbol}[GENE]&retmode=json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        variant_ids = data.get("esearchresult", {}).get("idlist", [])
        return variant_ids
    return []

def get_disease_variants_gwas(disease_name):
    """Retrieve disease-associated SNPs from GWAS Catalog."""
    url = f"https://www.ebi.ac.uk/gwas/rest/api/associations?query={disease_name}"
    response = requests.get(url, headers={"Accept": "application/json"})
    if response.status_code == 200:
        data = response.json()
        variants = [entry["rsId"] for entry in data.get("_embedded", {}).get("associations", []) if "rsId" in entry]
        return variants
    return []

def get_gene_variants_ensembl(gene_symbol):
    """Retrieve genetic variants from Ensembl for a given gene symbol."""
    url = f"https://rest.ensembl.org/overlap/id/{gene_symbol}?feature=variation;content-type=application/json"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        variants = [variant["id"] for variant in data if "id" in variant]
        return variants
    return []

def main():
    #disease_name = input("Enter the name of the autoimmune disease: ")
    disease_name = "rheumatoid arthritis"
    disease_id = get_disease_entry(disease_name)
    
    if not disease_id:
        print("No disease entry found for the given name.")
        return
    
    print(f"KEGG Disease ID: {disease_id}")
    genes = get_genes_associated_with_disease(disease_id)
    
    if not genes:
        print("No associated genes found.")
        return
    
    print("Associated Genes and Variants:")
    disease_variants = get_disease_variants_gwas(disease_name)
    if disease_variants:
        print(f"Disease-associated SNPs from GWAS Catalog: {', '.join(disease_variants)}")
    
    data_to_save = []
    for gene in genes:
        hsa_id = gene.split(":")[1]  # Extract gene ID
        gene_symbol = convert_hsa_to_gene_symbol(hsa_id)
        print(f"{gene} ({gene_symbol})")
        
        clinvar_variants = get_gene_variants_clinvar(gene_symbol)
        gwas_variants = get_disease_variants_gwas(disease_name)
        ensembl_variants = get_gene_variants_ensembl(gene_symbol)
        
        all_variants = set(clinvar_variants + gwas_variants + ensembl_variants)
        variant_str = " | ".join(all_variants) if all_variants else "No known variants"
        print(f"  Variants: {variant_str}")
        
        data_to_save.append([gene, gene_symbol, variant_str])
    
    with open("genetic_variants.csv", "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["KEGG Gene ID", "Gene Symbol", "Variants"])
        writer.writerows(data_to_save)
    
    print("Results saved to genetic_variants.csv")

if __name__ == "__main__":
    main()

