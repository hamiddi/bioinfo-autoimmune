"""
Author: Hamid D. Ismail, Ph.D.
Book: Bioinformatics of Autoimmune Diseases
############################################################
Protein-Protein Interaction Network Visualizer

This script reads a TSV or CSV file containing protein-protein interaction (PPI) data
and generates a graphical representation of the interaction network.

Usage Instructions:
1. Ensure your input file is named 'ppi_results.csv' and is located in the same directory as this script.
2. The file must contain at least two columns labeled 'Protein A' and 'Protein B',
   each representing an interaction between two proteins.
3. Run the script using a Python interpreter with required packages installed:
       python ppi_network_draw.py
4. The script will generate and display a network plot, and save it as 'ppi_network.png' in the current directory.

Required Libraries:
- pandas
- networkx
- matplotlib

Note:
If your file uses tabs as delimiters instead of commas, change the line:
    df = pd.read_csv("ppi_results.csv", sep=",")
to:
    df = pd.read_csv("ppi_results.csv", sep="\t")
"""

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Load TSV file
df = pd.read_csv("ppi_results.csv", sep=",")

# Create an undirected graph
G = nx.Graph()

# Add edges from the dataframe
for _, row in df.iterrows():
    G.add_edge(row['Protein A'], row['Protein B'])

# Draw the network
plt.figure(figsize=(14, 12))
pos = nx.spring_layout(G, seed=42)

nx.draw_networkx_nodes(G, pos, node_color='skyblue', node_size=1000)
nx.draw_networkx_edges(G, pos, width=2)
nx.draw_networkx_labels(G, pos, font_size=14, font_weight='bold')  # Increased font size

plt.title("Protein-Protein Interaction Network", fontsize=20)  # Larger title
plt.axis('off')
plt.tight_layout()
plt.savefig("ppi_network.png")
plt.show()

