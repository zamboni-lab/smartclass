"""TODO.

TODO.
"""

from __future__ import annotations

import logging
import re
import zipfile

# __all__ = [
#     "get_exact_relationships",
# ]


def get_exact_relationships(
    obo_file_path: str = "scratch/ChemOnt_2_1.obo.zip",
) -> dict[str, list[tuple[str, str]]]:
    """
    Get exact relationships.

    :param obo_file_path: Obo file path.
    :type obo_file_path: str

    :returns: Description.
    :rtype: dict[str, list[tuple[str, str]]]
    """
    chebi_synonyms: dict[str, list[tuple[str, str]]] = {}
    with zipfile.ZipFile(obo_file_path, "r") as z:
        with z.open(z.namelist()[0]) as obo_file:
            # Read and decode the entire file content
            content = obo_file.read().decode("utf-8")
            current_id = None
            for line in content.splitlines():
                line = line.strip()
                if line.startswith("id: "):
                    current_id = line.split(": ")[1]
                elif line.startswith("synonym: ") and "EXACT" in line:
                    synonym_match = re.match(
                        r"""synonym: "(.*?)" EXACT ChEBI_TERM \[CHEBI:(\d+)\]""", line
                    )
                    if synonym_match and current_id:
                        synonym = synonym_match.group(1)
                        chebi_id = synonym_match.group(2)
                        chebi_synonyms[current_id] = [
                            *chebi_synonyms.get(current_id, []),
                            (synonym, chebi_id),
                        ]
    return chebi_synonyms


if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(level=logging.DEBUG, format="%(levelname)s: %(message)s")

    # Create a dictionary to store exact synonyms
    chebi_synonyms = get_exact_relationships()

    # Print the ChEBI IDs and exact synonyms
    for term_id, synonyms in chebi_synonyms.items():
        logging.debug("CHEMONT ID: %s", term_id)
        for _synonym, chebi_id in synonyms:
            logging.debug("ChEBI ID: CHEBI: %s", chebi_id)
