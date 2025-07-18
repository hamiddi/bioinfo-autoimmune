"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
###########################################################
KEGG Autoimmune Disease Extractor
----------------------------------
This script retrieves structured KEGG data for a specified autoimmune disease, including:
- KEGG disease description and ID
- Associated KEGG pathways
- Human gene IDs and names involved in those pathways

Usage:
    1. Set the `disease` variable at the bottom of the script to the desired disease name
       (e.g., "rheumatoid arthritis", "systemic lupus erythematosus").
    2. Run the script using Python 3 (internet connection required).
    3. The script will:
        - Query the KEGG API to collect disease-related data
        - Save the data to a JSON file (e.g., rheumatoid_arthritis.json)
        - Save a formatted plain-text summary (e.g., rheumatoid_arthritis.txt)

Note:
    - The script includes sleep delays to comply with KEGG’s rate-limiting.
    - Ensure you have the `requests` module installed (`pip install requests`).
"""
import requests
import re
import time
import json

def fetch_kegg_disease_info(disease_name):
    """
    Fetches KEGG disease description, pathway IDs, and associated gene names/IDs for a given autoimmune disease.
    
    Parameters:
        disease_name (str): Name of the autoimmune disease (e.g., "rheumatoid arthritis")
    
    Saves:
        - JSON file with disease name.
        - Text file with structured format.
    
    Returns:
        dict: A structured dictionary containing KEGG disease information, pathways, and associated genes.
    """
    results = {"KEGG_Disease": None, "KEGG_Pathways": [], "KEGG_Genes": []}
    disease_filename = disease_name.replace(" ", "_")  # Ensure filename safety

    #Step 1: Find KEGG Disease ID
    kegg_disease_url = f"https://rest.kegg.jp/find/disease/{disease_name}"
    
    try:
        response = requests.get(kegg_disease_url)
        response.raise_for_status()
        diseases = response.text.strip().split("\n")

        disease_id = None
        if diseases:
            first_disease = diseases[0].split("\t")
            if len(first_disease) > 1:
                disease_id = first_disease[0].replace("ds:", "")
                disease_desc = first_disease[1]
                results["KEGG_Disease"] = {"disease_id": disease_id, "description": disease_desc}
        
    except requests.exceptions.RequestException as e:
        print(f"KEGG API error: {e}")
        return results

    #Step 2: Find Pathways Associated with the Disease
    if disease_id:
        kegg_pathway_url = f"https://rest.kegg.jp/link/pathway/{disease_id}"
        
        try:
            response = requests.get(kegg_pathway_url)
            response.raise_for_status()
            pathways = response.text.strip().split("\n")

            pathway_ids = []
            for line in pathways:
                parts = line.split("\t")
                if len(parts) > 1:
                    pathway_id = parts[1].replace("path:", "")
                    pathway_ids.append(pathway_id)
                    results["KEGG_Pathways"].append({"pathway_id": pathway_id})
        
        except requests.exceptions.RequestException as e:
            print(f"❌ KEGG Pathway API error: {e}")

    #Step 3: Extract Genes from Each Pathway
    gene_ids = set()
    for pathway_id in pathway_ids:
        kegg_pathway_gene_url = f"https://rest.kegg.jp/link/hsa/{pathway_id}"
        
        try:
            response = requests.get(kegg_pathway_gene_url)
            response.raise_for_status()
            gene_data = response.text.strip().split("\n")

            for line in gene_data:
                parts = line.split("\t")
                if len(parts) > 1 and "hsa:" in parts[1]:  
                    gene_id = parts[1].replace("hsa:", "")
                    gene_ids.add(gene_id)  

        except requests.exceptions.RequestException as e:
            print(f"❌ KEGG Pathway Gene API error for {pathway_id}: {e}")

    #Step 4: Retrieve Gene Names with Improved Extraction
    gene_data_list = []
    for gene_id in gene_ids:
        kegg_gene_url = f"https://rest.kegg.jp/get/hsa:{gene_id}"
        
        try:
            response = requests.get(kegg_gene_url)
            response.raise_for_status()
            gene_info = response.text

            # Extract gene name
            match = re.search(r"NAME\s+([^\n;]+)", gene_info)
            gene_name = match.group(1).strip() if match else "Unknown"

            # If still "Unknown", try alternative extraction
            if gene_name == "Unknown":
                alt_match = re.search(r"DESCRIPTION\s+([^\n;]+)", gene_info)
                if alt_match:
                    gene_name = alt_match.group(1).strip()

            gene_data_list.append(f"{gene_id}: {gene_name}")
            results["KEGG_Genes"].append({"gene_id": gene_id, "gene_name": gene_name})

            # Delay to avoid KEGG rate-limiting
            time.sleep(1)

        except requests.exceptions.RequestException as e:
            print(f"KEGG Gene API error for {gene_id}: {e}")

    #Save JSON file
    json_filename = f"{disease_filename}.json"
    with open(json_filename, "w", encoding="utf-8") as json_file:
        json.dump(results, json_file, indent=4)
    print(f"JSON data saved to {json_filename}")

    #Save formatted text file
    text_filename = f"{disease_filename}.txt"
    with open(text_filename, "w", encoding="utf-8") as text_file:
        text_file.write(f"disease_id: {results['KEGG_Disease']['disease_id']}\n")
        text_file.write(f"description: {results['KEGG_Disease']['description']}\n")
        text_file.write(f"pathway_id: {' | '.join([p['pathway_id'] for p in results['KEGG_Pathways']])}\n")
        text_file.write(f"KEGG_Genes: {' | '.join(gene_data_list)}\n")
    print(f"Formatted text saved to {text_filename}")

    return results

#Example usage
disease = "rheumatoid arthritis"
fetch_kegg_disease_info(disease)

