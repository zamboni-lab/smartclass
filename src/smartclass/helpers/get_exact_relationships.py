"""TODO.

TODO.
"""

# MISSING DEF
# __all__ = [
#     "get_exact_relationships",
# ]
from __future__ import annotations

import logging
import re

# Define the path to your OBO file
obo_file_path = "data/ChemOnt_2_1.obo"

# Create a dictionary to store exact synonyms
chebi_synonyms: dict[str, list[tuple[str, str]]] = {}

# Open and read the OBO file
with open(obo_file_path) as obo_file:
    current_id = None
    for line in obo_file:
        line = line.strip()
        if line.startswith("id: "):
            current_id = line.split(": ")[1]
        elif line.startswith("synonym: ") and "EXACT" in line:
            synonym_match = re.match(r"""synonym: "(.*?)" EXACT ChEBI_TERM \[CHEBI:(\d+)\]""", line)
            if synonym_match and current_id:
                synonym = synonym_match.group(1)
                chebi_id = synonym_match.group(2)
                chebi_synonyms[current_id] = [
                    *chebi_synonyms.get(current_id, []),
                    (synonym, chebi_id),
                ]
logging.debug(chebi_synonyms.items())
# Print the ChEBI IDs and exact synonyms
for term_id, synonyms in chebi_synonyms.items():
    logging.debug("CHEMONT ID: %s", term_id)
    for _synonym, chebi_id in synonyms:
        logging.debug("ChEBI ID: CHEBI: %s", chebi_id)
    logging.debug("")
