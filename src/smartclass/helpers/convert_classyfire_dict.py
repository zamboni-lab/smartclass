"""Converts the classyfire json into a CHEMONTID dictionary."""

from __future__ import annotations

import json
import logging

__all__ = [
    "convert_classyfire_dict",
]


def convert_classyfire_dict(
    classyfire_json: str = "scratch/classyfire.json",
    output: str = "scratch/chemontids.txt",
) -> None:
    """
    Converts the classyfire json into a CHEMONTID dictionary.

    :param classyfire_json: Path of the classyfire json.
    :type classyfire_json: str

    :param output: Output.
    :type output: str
    """
    with open(classyfire_json) as json_file:
        data = json.load(json_file)

    with open(output, "w") as txt_file:
        for item in data:
            chemont_id = item.get("chemont_id")
            if chemont_id:
                txt_file.write(chemont_id + "\n")

    logging.debug(f"CHEMONTIDs have been written to {output}")
