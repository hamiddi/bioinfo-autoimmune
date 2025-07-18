"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This script reads a JSON file containing biological records, corrects formatting issues if the JSON
objects are not properly enclosed in a list, parses the content, and writes selected fields to a CSV file.
Each entry includes accession, protein, gene, organism, and keywords, which are concatenated with a pipe (|) symbol.
"""
import json
import csv
# Define file paths
json_file_path = "selected_fields.json"  # Update the path if necessary
csv_file_path = "output.csv"
# Read the JSON file and fix formatting issues
with open(json_file_path, "r", encoding="utf-8") as file:
    content = file.read()
# Ensure JSON objects are correctly structured as a list
content = "[" + content.replace("}\n{", "},\n{") + "]"
try:
    data = json.loads(content)  # Convert corrected content to JSON
except json.JSONDecodeError as e:
    print(f"Error parsing JSON: {e}")
    exit(1)
# Open CSV file for writing
with open(csv_file_path, "w", newline="", encoding="utf-8") as csv_file:
    fieldnames =["Accession","Protein","Gene","Organism","Keywords"]
    writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
    writer.writeheader()  # Write CSV header
    # Process each JSON object and write to CSV
    for entry in data:
        keywords = "|".join([kw["name"] for kw in entry.get("Keywords", [])])
        writer.writerow({
            "Accession": entry["Accession"],
            "Protein": entry["Protein"],
            "Gene": entry["Gene"],
            "Organism": entry["Organism"],
            "Keywords": keywords
        })
print(f"CSV file has been created successfully: {csv_file_path}")

