"""TODO.

TODO.
"""

# MISSING DEF
# __all__ = [
#     "classyfire_to_dict",
# ]
from __future__ import annotations

import json
import logging

# Load the JSON data from the file
with open("data/classyfire.json") as json_file:
    data = json.load(json_file)

# Extract CHEMONTIDs and write them to a text file
with open("chemontids.txt", "w") as txt_file:
    for item in data:
        chemont_id = item.get("chemont_id")
        if chemont_id:
            txt_file.write(chemont_id + "\n")

logging.debug("CHEMONTIDs have been written to chemontids.txt")
