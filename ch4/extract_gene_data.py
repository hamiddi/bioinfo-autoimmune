"""
Book title: Bioinformatics of Autoimmune Diseases
Author: Hamid D. Ismail, Ph.D.
############################################################
This Python script extracts specific gene and genomic information from a collection of NCBI-style .txt files
and compiles it into a structured CSV file. It is designed to work with .txt files formatted with key-value 
pairs, typically containing metadata about genes.

USAGE:
1. Place all the .txt files to be processed in a directory (default is "SLE_gene_ids").
2. Modify the 'directory' variable at the bottom of the script if your files are in a different folder.
3. Run the script using a Python interpreter (e.g., `python extract_gene_info.py`).
4. The script will parse each .txt file, extract relevant fields (e.g., gene name, organism, chromosome position),
   and save the results to a CSV file named "extracted_data.csv" in the current directory.

Note:
- The script handles nested "genomicinfo" fields and maps renamed keys like "organism scientificname" to "organism".
- Fields with missing data are filled with "NA".
"""

import os
import ast
import csv

# Define the keys to extract (mapping "organism scientificname" → "organism", "organism taxid" → "taxid")
KEYS = [
    "uid", "name", "description", "chromosome", "maplocation", 
    "otherdesignations", "nomenclaturesymbol", "nomenclaturename", "mim", 
    "chraccver", "chrstart", "chrstop", "exoncount", 
    "organism", "taxid"
]

# Mapping for renamed fields
FIELD_MAPPING = {
    "organism scientificname": "organism",
    "organism taxid": "taxid"
}

def parse_txt_file(file_path):
    """Extract relevant information from a .txt file."""
    data = {key: "NA" for key in KEYS}  # Default values as NA
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            
            key, value = line.split(": ", 1) if ": " in line else (line, "")
            
            # Handle direct key matches
            if key in KEYS:
                data[key] = value

            # Handle renamed fields
            elif key in FIELD_MAPPING:
                data[FIELD_MAPPING[key]] = value
            
            # Handle genomicinfo parsing
            elif key == "genomicinfo":
                try:
                    genomic_info = ast.literal_eval(value)
                    if isinstance(genomic_info, dict):
                        data["chraccver"] = str(genomic_info.get("chraccver", "NA"))
                        chrstart = genomic_info.get("chrstart", "NA")
                        chrstop = genomic_info.get("chrstop", "NA")
                        exoncount = genomic_info.get("exoncount", "NA")

                        # Ensure chrstart is actually smaller than chrstop
                        if chrstart != "NA" and chrstop != "NA":
                            chrstart, chrstop = int(chrstart), int(chrstop)
                            if chrstart > chrstop:  # Fix the reversed values
                                chrstart, chrstop = chrstop, chrstart

                        data["chrstart"] = str(chrstart)
                        data["chrstop"] = str(chrstop)
                        data["exoncount"] = str(exoncount)
                except (SyntaxError, ValueError):
                    data["chraccver"] = data["chrstart"] = data["chrstop"] = data["exoncount"] = "NA"
    
    return data

def extract_data_from_directory(directory):
    """Process all .txt files in the given directory."""
    output_data = []
    
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            file_path = os.path.join(directory, filename)
            extracted_data = parse_txt_file(file_path)
            output_data.append(extracted_data)
    
    return output_data

def save_to_csv(data, output_file):
    """Save extracted data to a CSV file."""
    with open(output_file, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=KEYS)
        writer.writeheader()
        writer.writerows(data)



# Directory containing .txt files
directory = "SLE_gene_ids"  # Change this to your target directory
output_csv = "extracted_data.csv"

# Extract data and save to CSV
extracted_data = extract_data_from_directory(directory)
save_to_csv(extracted_data, output_csv)

print(f"Data extraction complete. Saved to {output_csv}")




