"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script retrieves biological pathway, gene, and enzyme information related to a specified disease 
using the KEGG REST API. It is particularly useful for exploring the molecular mechanisms associated 
with diseases such as autoimmune disorders.

How to use:
- Set the `disease_name` variable in the `main()` function to the disease of interest 
  (e.g., "rheumatoid arthritis").
- Alternatively, uncomment the `input()` line to allow user input at runtime.
- Run the script using Python 3: `python script_name.py`
- The program will:
    1. Find the KEGG disease ID for the given disease name.
    2. Retrieve associated pathways, human genes, and enzymes from KEGG.
    3. Save the results in both `results.json` and `results.tsv` formats.

Dependencies:
- Python 3
- requests
- json
- csv

Note:
- Ensure you have an internet connection, as the script fetches data from KEGG's public API.
"""

import requests
import json
import csv

def get_disease_info(disease_name):
    search_url = f"https://rest.kegg.jp/find/disease/{disease_name}"
    response = requests.get(search_url)
    
    if response.status_code != 200 or not response.text:
        print("No disease information found.")
        return None
    
    disease_data = response.text.split('\n')[0]
    if '\t' not in disease_data:
        print("Unexpected response format.")
        return None

    disease_id = disease_data.split('\t')[0]
    print(f"Found KEGG disease ID: {disease_id}")
    return disease_id

def get_pathways_for_disease(disease_id):
    pathways_url = f"https://rest.kegg.jp/link/pathway/{disease_id}"
    response = requests.get(pathways_url)
    
    if response.status_code != 200 or not response.text:
        print("No pathways found for the disease.")
        return []
    
    pathways = [line.split('\t')[1] for line in response.text.strip().split('\n') if '\t' in line]
    print(f"Identified pathways: {pathways}")
    return pathways

def get_genes_for_pathways(pathways):
    genes = set()
    for pathway in pathways:
        pathway_url = f"https://rest.kegg.jp/link/hsa/{pathway}"
        response = requests.get(pathway_url)
        
        if response.status_code == 200 and response.text:
            for line in response.text.strip().split('\n'):
                if '\t' in line:
                    genes.add(line.split('\t')[1])

    if not genes:
        print("No genes found for the pathways.")
    else:
        print(f"Identified genes: {genes}")
    return genes

def filter_enzymes(genes):
    enzymes = {}
    for gene in genes:
        enzyme_url = f"https://rest.kegg.jp/link/enzyme/{gene}"
        response = requests.get(enzyme_url)

        if response.status_code == 200 and response.text:
            for line in response.text.strip().split('\n'):
                if '\t' in line:
                    enzyme_id = line.split('\t')[1]
                    enzyme_name = get_enzyme_name(enzyme_id)
                    enzymes[enzyme_id] = enzyme_name

    if not enzymes:
        print("No enzymes found for the genes.")
    else:
        print(f"Enzymes linked to genes: {enzymes}")
    return enzymes

def get_enzyme_name(enzyme_id):
    enzyme_url = f"https://rest.kegg.jp/get/{enzyme_id}"
    response = requests.get(enzyme_url)
    
    if response.status_code == 200 and response.text:
        for line in response.text.strip().split('\n'):
            if line.startswith("NAME"):
                # Extract enzyme name correctly (may have multiple names)
                enzyme_name = line.split("        ")[1].strip().split(";")[0]
                return enzyme_name
    return "Unknown"

def save_results(disease_name, pathways, genes, enzymes):
    data = {
        "disease": disease_name,
        "pathways": pathways,
        "genes": list(genes),
        "enzymes": enzymes
    }
    
    # Save to JSON
    with open("results.json", "w") as json_file:
        json.dump(data, json_file, indent=4)
    
    # Save to TSV
    with open("results.tsv", "w", newline="") as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerow(["Category", "Data"])
        writer.writerow(["Disease", disease_name])
        writer.writerow(["Pathways", ", ".join(pathways)])
        writer.writerow(["Genes", ", ".join(genes)])
        writer.writerow(["Enzymes", ", ".join([f"{eid} ({ename})" for eid, ename in enzymes.items()])])
    
    print("Results saved to results.json and results.tsv")

def main():
    #disease_name = input("Enter an autoimmune disease (e.g., rheumatoid arthritis): ").strip()
    disease_name = "rheumatoid arthritis"
    disease_id = get_disease_info(disease_name)
    
    if disease_id:
        pathways = get_pathways_for_disease(disease_id)
        if pathways:
            genes = get_genes_for_pathways(pathways)
            enzymes = filter_enzymes(genes)
            
            print("\nFinal Results:")
            print(f"Pathways involved: {pathways}")
            print(f"Genes involved: {genes}")
            print(f"Enzymes involved: {enzymes}")
            
            save_results(disease_name, pathways, genes, enzymes)
        else:
            print("No pathways retrieved.")
    else:
        print("No disease information found.")

if __name__ == "__main__":
    main()

