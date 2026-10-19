"""Converts the classyfire json into a CHEMONTID dictionary."""

from __future__ import annotations

import json

from smartclass.logging import get_logger

__all__ = [
    "convert_classyfire_dict",
]

logger = get_logger(__name__)


def convert_classyfire_dict(
    classyfire_json: str = "scratch/classyfire.json",
    output: str = "scratch/chemontids.txt",
) -> None:
    """Converts the classyfire json into a CHEMONTID dictionary.

    Parameters
    ----------
    classyfire_json : str
        Default is 'scratch/classyfire.json'.
    output : str
        Default is 'scratch/chemontids.txt'.
    """
    with open(classyfire_json) as json_file:
        data = json.load(json_file)

    with open(output, "w") as txt_file:
        for item in data:
            chemont_id = item.get("chemont_id")
            if chemont_id:
                txt_file.write(chemont_id + "\n")

    logger.debug(f"CHEMONTIDs have been written to {output}")
