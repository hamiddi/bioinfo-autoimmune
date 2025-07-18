# Author: Hamid D. Ismail, Ph.D.
# Book: Bioinformatics of Autoimmune Diseases
# This script reads a JSON file containing genomic information, extracts key biological attributes 
# such as organism, chromosome, and sequence data, and prints them along with related publications.
import json
# Define the file path for the JSON file
file_path = "example.json"
# Load and parse the JSON file
with open(file_path, "r", encoding="utf-8") as file:
    data = json.load(file)
# Extract relevant details
organism = data.get("Organism", "Not Found")
common_name = data.get("CommonName", "Not Found")
chromosome = data.get("Chromosome", "Not Found")
location = data.get("Location", "Not Found")
sequence = data.get("Sequence", "Not Found")
publications = data.get("Publications", [])
# Display the extracted information
print(f"Organism: {organism}")
print(f"Common Name: {common_name}")
print(f"Chromosome: {chromosome}")
print(f"Location: {location}")
print(f"Sequence: {sequence}\n")
print("Publications:")
for pub in publications:
    print(f"  PMID: {pub['PMID']}")
    print(f"  Title: {pub['Title']}\n")

