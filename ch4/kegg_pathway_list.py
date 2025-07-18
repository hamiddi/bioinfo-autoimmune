"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script retrieves the complete list of KEGG pathways using the bioservices package
and saves the results to a text file named 'kegg_pathways.txt'.

Usage:
- Ensure that the 'bioservices' package is installed. You can install it via pip:
    pip install bioservices
- Run the script in a Python environment:
    python script_name.py
- The list of KEGG pathways will be printed to the console and saved to 'kegg_pathways.txt'
  in the current working directory.

Output:
- A text file 'kegg_pathways.txt' containing the KEGG pathway IDs and descriptions.
"""
from bioservices import KEGG
def kegg_pathway_list():
  kegg = KEGG()
  pathways = kegg.list("pathway")
  pathways_list = [line for line in pathways.split("\n")]
  with open("kegg_pathways.txt", "w") as f:
     f.write("\n".join(pathways_list))
  print("\n".join(pathways_list))
  print("Pathway list saved to kegg_pathways.txt")
kegg_pathway_list()
