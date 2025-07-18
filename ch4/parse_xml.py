"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases

This script parses an XML file containing genomic information, extracts key fields such as organism,
chromosome, and sequence, and retrieves publication data including PMID and titles. The extracted 
information is then printed in a structured format.
"""

import xml.etree.ElementTree as ET
# Load and parse the XML file
file_path = "example.xml"
tree = ET.parse(file_path)
root = tree.getroot()

# Extract essential data
data = {
    "Organism": root.find("Organism").text,
    "Common Name": root.find("CommonName").text,
    "Chromosome": root.find("Chromosome").text,
    "Location": root.find("Location").text,
    "Sequence": root.find("Sequence").text
}
# Extract publication data
pubs = {}
for pub in root.findall("Publications/Publication"):
    pmid = pub.find("PMID").text
    title = pub.find("Title").text
    pubs[pmid] = title
for key in data:
    print(f"{key}: {data[key]}")
for key in pubs:
    print(f"{key} : {pubs[key]}")

