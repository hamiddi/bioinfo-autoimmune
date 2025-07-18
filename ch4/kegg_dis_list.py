"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
This script fetches the complete list of diseases from the KEGG (Kyoto Encyclopedia of Genes and Genomes) database
using the BioServices Python package and saves the list to a local text file named 'kegg_diseases.txt'.

How to use:
1. Ensure that the required packages are installed:
   - bioservices
   - requests (though not used in this script, it's imported)

   You can install them using pip:
   pip install bioservices requests

2. Run the script using Python:
   python fetch_kegg_diseases.py

3. After execution, the disease list will be saved in the current directory as 'kegg_diseases.txt'.
"""

from bioservices import KEGG
def fetch_Kegg_dis_list():
  kegg = KEGG()
  disease_list = kegg.list("disease")
  with open("kegg_diseases.txt", "w") as f:
     f.write(disease_list)
  print("Disease list saved to kegg_diseases.txt")
fetch_Kegg_dis_list()
