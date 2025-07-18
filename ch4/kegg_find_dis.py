"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script allows you to search for disease entries in the KEGG database using a keyword.
Usage:
- Run the script and provide a keyword related to a disease (e.g., "autoimmune").
- The program will search the KEGG disease database and print all matching entries.

Requirements:
- Install the 'bioservices' library before running the script:
    pip install bioservices

Example:
    The function call `find_kegg_dis("autoimmune")` will print KEGG disease entries related to autoimmune disorders.
"""
from bioservices import KEGG
def find_kegg_dis(keyword):
   kegg = KEGG()
   diseases = kegg.find("disease", keyword)
   print(diseases)
find_kegg_dis("autoimmune")
