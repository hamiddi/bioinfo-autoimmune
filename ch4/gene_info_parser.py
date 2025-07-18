"""
Book title: Bioinformatics of Autoimmune Diseases
Author: Hamid D. Ismail, Ph.D.
############################################################
This script fetches gene summaries from the NCBI Gene database using the Entrez E-Utilities API.
It reads gene IDs from a text file (one ID per line), retrieves their summaries in JSON format,
parses and formats the content into readable text, and saves each gene summary to a separate text
file in a specified output directory.

Instructions:
1. Ensure you have the Biopython library installed (`pip install biopython`).
2. Create a text file (e.g., "SLE_gene_ids.txt") containing one NCBI Gene ID per line.
3. Set the `input_file` and `output_dir` variables at the bottom of this script to match your files.
4. Run the script. It will:
   - Fetch data from NCBI for each gene ID.
   - Format the summary.
   - Save the results as .txt files in the output directory.

Note:
- Replace `"your_email@example.com"` with your actual email address in the `fetch_gene_summary` function.
  NCBI requires this to monitor API usage.
"""

import os
from Bio import Entrez
import json

# Function to fetch gene summary from NCBI in JSON format
def fetch_gene_summary(gene_id, email="your_email@example.com"):
    Entrez.email = email  # Provide a valid email
    try:
        handle = Entrez.esummary(db="gene", id=gene_id, retmode="json")
        gene_summary = json.load(handle)
        handle.close()
        return gene_summary
    except Exception as e:
        print(f"Error fetching gene {gene_id}: {e}")
        return None

# Function to format JSON data into key-value pairs on the same line
def format_json_to_text(json_obj):
    text_output = ""
    
    def parse_json(data, prefix=""):
        nonlocal text_output
        if isinstance(data, dict):
            for key, value in data.items():
                new_prefix = f"{prefix}{key} "
                parse_json(value, new_prefix)
        elif isinstance(data, list):
            values = "; ".join(str(item) for item in data) if data else "NA"
            text_output += f"{prefix.strip()}: {values}\n"
        else:
            text_output += f"{prefix.strip()}: {data if data else 'NA'}\n"
    
    parse_json(json_obj)
    return text_output

# Function to save gene summary to a file
def save_gene_summary(gene_id, output_dir):
    gene_data = fetch_gene_summary(gene_id)
    
    if gene_data:
        try:
            gene_info = gene_data["result"].get(str(gene_id), {})
            formatted_text = format_json_to_text(gene_info)

            # Create a file named after the gene ID
            file_path = os.path.join(output_dir, f"{gene_id}.txt")
            with open(file_path, "w", encoding="utf-8") as file:
                file.write(formatted_text)
            
            print(f"Saved: {file_path}")
        except KeyError:
            print(f"Error: Gene data not properly structured for ID {gene_id}")
    else:
        print(f"Failed to fetch gene summary for ID {gene_id}")

# Function to process multiple gene IDs from a file
def batch_fetch_genes(input_file, output_dir):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read gene IDs from file
    with open(input_file, "r", encoding="utf-8") as file:
        gene_ids = [line.strip() for line in file if line.strip()]

    print(f"Fetching data for {len(gene_ids)} genes...")
    
    for gene_id in gene_ids:
        save_gene_summary(gene_id, output_dir)

# Example usage
input_file = "SLE_gene_ids.txt"  # File containing gene IDs, one per line
output_dir = "SLE_gene_ids"  # Directory for output files

batch_fetch_genes(input_file, output_dir)

