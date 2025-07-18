#!/bin/bash
# Book title: Bioinformatics of Autoimmune Diseases
# Author: Hamid D. Ismail, Ph.D.
URL="https://rest.uniprot.org/uniprotkb/search"
PARAM="?query=autoimmune+disease&format=json"
curl -X GET \
     "${URL}${PARAM}" \
     -H "Accept: application/json" \
     > uniprot_autoimmune.json
jq '
  .results[] |
  {
    Accession: .primaryAccession,
    Protein: .proteinDescription.recommendedName.fullName.value,
    Gene: .genes[0].geneName.value,
    Organism: .organism.commonName,
    Keywords: .keywords
  }
' uniprot_autoimmune.json > selected_fields.json

