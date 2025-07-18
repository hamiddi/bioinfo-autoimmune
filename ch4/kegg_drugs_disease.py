"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script fetches drug-target interaction data for a specified autoimmune disease
from the KEGG DRUG database using KEGG's REST API. It retrieves drug entries linked
to the disease, extracts their target proteins, and saves the data in both JSON and TSV formats.

How to Use:
1. Ensure you have an active internet connection.
2. Install the required Python libraries (requests).
3. Call the `fetch_kegg_drug_target_interactions()` function with the name of an autoimmune disease
   (e.g., "rheumatoid arthritis").
4. The function will save the results to:
    - A JSON file named "drug_targets.json"
    - A TSV file named "drug_targets.tsv"
5. The function also returns the parsed data as a Python list of dictionaries.

Example:
    data = fetch_kegg_drug_target_interactions("rheumatoid arthritis")

Dependencies:
- requests
- json
- csv
- xml.etree.ElementTree (standard library)
"""
import requests
import json
import csv
import xml.etree.ElementTree as ET

def fetch_kegg_drug_target_interactions(disease_name, json_filename="drug_targets.json", tsv_filename="drug_targets.tsv"):
    """
    Fetch drug-target interactions for a given autoimmune disease from KEGG DRUG database and save the data.
    
    Parameters:
        disease_name (str): Name of the autoimmune disease.
        json_filename (str): Name of the JSON file to save results.
        tsv_filename (str): Name of the TSV file to save results.
    """
    base_url = "http://rest.kegg.jp"
    
    # Search for disease entry in KEGG
    disease_search_url = f"{base_url}/find/disease/{disease_name}"
    response = requests.get(disease_search_url)
    
    if response.status_code != 200:
        raise ValueError("Failed to retrieve data from KEGG.")
    
    disease_entries = response.text.strip().split("\n")
    if not disease_entries:
        return []
    
    disease_ids = [entry.split("\t")[0] for entry in disease_entries if "H" in entry]
    
    drug_target_data = []
    
    for disease_id in disease_ids:
        # Get drugs associated with the disease
        drugs_url = f"{base_url}/link/drug/{disease_id}"
        response = requests.get(drugs_url)
        
        if response.status_code != 200:
            continue
        
        drug_entries = response.text.strip().split("\n")
        drug_ids = [entry.split("\t")[1] for entry in drug_entries]
        
        for drug_id in drug_ids:
            # Get target information for the drug
            target_url = f"{base_url}/get/{drug_id}"
            response = requests.get(target_url)
            
            if response.status_code != 200:
                continue
            
            drug_info = response.text
            
            # Extract relevant details
            lines = drug_info.split("\n")
            drug_name = "Unknown"
            targets = []
            
            current_target = ""
            
            for line in lines:
                if line.startswith("NAME"):
                    drug_name = line.split("   ")[-1].strip()
                elif line.startswith("TARGET"):
                    current_target = line.split("   ")[-1].strip()
                    targets.append({"target": current_target, "details": []})
                elif current_target and line.startswith("  "):
                    targets[-1]["details"].append(line.strip())
            
            drug_target_data.append({
                "disease": disease_name,
                "disease_id": disease_id,
                "drug_id": drug_id,
                "drug_name": drug_name,
                "targets": targets
            })
    
    # Save data as JSON
    with open(json_filename, "w", encoding="utf-8") as json_file:
        json.dump(drug_target_data, json_file, indent=4)
    
    # Save data as TSV
    with open(tsv_filename, "w", encoding="utf-8", newline="") as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter="\t")
        tsv_writer.writerow(["disease", "disease_id", "drug_id", "drug_name", "targets"])
        
        for entry in drug_target_data:
            targets_str = "; ".join([t["target"] + " (" + ", ".join(t["details"]) + ")" if t["details"] else t["target"] for t in entry["targets"]])
            tsv_writer.writerow([entry["disease"], entry["disease_id"], entry["drug_id"], entry["drug_name"], targets_str])
    
    return drug_target_data

# Example usage
disease = "rheumatoid arthritis"
data = fetch_kegg_drug_target_interactions(disease)
for entry in data:
    print(entry)

