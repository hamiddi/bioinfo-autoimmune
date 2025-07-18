"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script fetches a list of autoimmune diseases from the KEGG database using the KEGG REST API.
It filters the results against a predefined list of known autoimmune disorders and saves the
matched entries in two formats:

1. `autoimmune_diseases.json`: A structured JSON file containing KEGG disease IDs and names.
2. `autoimmune_diseases.txt`: A plain-text file with readable disease information.

Usage:
    - Simply run the script with Python 3.
    - Make sure you have an active internet connection.
    - Required library: `requests` (install via `pip install requests` if not already installed).

Output files will be saved in the current working directory.
"""
import requests
import json

# Predefined list of autoimmune diseases to match from KEGG
autoimmune_diseases_list = [
    "Rheumatoid arthritis",
    "Systemic lupus erythematosus",
    "Multiple sclerosis",
    "Type 1 diabetes mellitus",
    "Psoriasis",
    "Ankylosing spondylitis",
    "Sjogren syndrome",
    "Inflammatory bowel disease",
    "Ulcerative colitis",
    "Crohn's disease",
    "Autoimmune thyroid disease",
    "Hashimoto thyroiditis",
    "Graves disease",
    "Celiac disease",
    "Autoimmune hemolytic anemia",
    "Goodpasture syndrome",
    "Myasthenia gravis",
    "Guillain-Barre syndrome"
]

def fetch_immune_disease_pathways():
    """
    Fetches KEGG disease pathways for known autoimmune diseases.

    Saves:
        - JSON file: Contains KEGG disease IDs and names for autoimmune diseases.
        - Text file: Readable list of diseases.

    Returns:
        dict: A structured dictionary containing autoimmune disease pathways.
    """
    results = {"Autoimmune_Diseases": []}
    kegg_disease_url = "https://rest.kegg.jp/list/disease"

    try:
        response = requests.get(kegg_disease_url)
        response.raise_for_status()
        lines = response.text.split("\n")

        for line in lines:
            if line.strip():
                parts = line.split("\t")
                if len(parts) == 2:
                    disease_id = parts[0].replace("ds:", "").strip()
                    disease_name = parts[1].strip()

                    # Match against known autoimmune diseases (exact match)
                    if disease_name in autoimmune_diseases_list:
                        results["Autoimmune_Diseases"].append(
                            {"disease_id": disease_id, "disease_name": disease_name}
                        )

    except requests.exceptions.RequestException as e:
        print(f"KEGG API error: {e}")
        return results

    if results["Autoimmune_Diseases"]:
        json_filename = "autoimmune_diseases.json"
        with open(json_filename, "w", encoding="utf-8") as json_file:
            json.dump(results, json_file, indent=4)

        text_filename = "autoimmune_diseases.txt"
        with open(text_filename, "w", encoding="utf-8") as text_file:
            for disease in results["Autoimmune_Diseases"]:
                text_file.write(f"disease_id: {disease['disease_id']}\n")
                text_file.write(f"disease_name: {disease['disease_name']}\n\n")

        print(f"Data successfully saved to {json_filename} and {text_filename}")

    else:
        print("No autoimmune diseases found in KEGG.")

    return results

fetch_immune_disease_pathways()

