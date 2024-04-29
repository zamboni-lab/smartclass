"""Extract CHEBI."""

from __future__ import annotations

import re

from smartclass.io import export_results


def extract_chebi(file_path: str = "scratch/chebi.obo", output: str = "scratch/test.tsv"):
    """
    Extract CHEBI.

    :param file_path: Input.
    :type file_path: str

    :param output: Output.
    :type output: str
    """
    # Define regular expressions to match relevant lines in the OBO file
    term_start_pattern = re.compile(r"^\[Term\]")
    id_pattern = re.compile(r"^id: (.+)")
    inchikey_pattern = re.compile(r"^property_value: .+?/chebi/inchikey \"(.+)\" xsd:string")
    is_a_pattern = re.compile(r"^is_a: (.+)")

    # Initialize a list to store dictionaries of extracted data
    term_data_list = []
    current_id = None
    current_inchikey = None
    current_is_a = None

    # Read the OBO file
    with open(file_path) as file:
        for line in file:
            line = line.strip()
            if term_start_pattern.match(line):
                if current_id:
                    term_data_list.append(
                        {
                            "id": current_id.replace("CHEBI:", "") if current_id else None,
                            "inchikey": current_inchikey if current_inchikey else None,
                            "parent": current_is_a.replace("CHEBI:", "") if current_is_a else None,
                        }
                    )
                current_id = None
                current_inchikey = None
                current_is_a = None
            elif id_match := id_pattern.match(line):
                current_id = id_match.group(1)
            elif inchikey_match := inchikey_pattern.match(line):
                current_inchikey = inchikey_match.group(1)
            elif is_a_match := is_a_pattern.match(line):
                current_is_a = is_a_match.group(1)

    # Add the last term data
    if current_id:
        term_data_list.append(
            {
                "id": current_id.replace("CHEBI:", "") if current_id else None,
                "inchikey": current_inchikey if current_inchikey else None,
                "parent": current_is_a.replace("CHEBI:", "") if current_is_a else None,
            }
        )

    # Export the list of dictionaries
    export_results(output=output, results=term_data_list)
