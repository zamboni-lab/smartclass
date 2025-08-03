"""TODO.

TODO.
"""

# MISSING DEF
# __all__ = [
#     "extract_smiles_lotus",
# ]
from __future__ import annotations

import gzip
import os

import polars

path_lotus = "~/Git/lotus-processor/data/processed/230106_frozen_metadata.csv.gz"
path_lotus_smiles_out = "~/Git/smartclass/data/lotus_smiles.csv.gz"

# Expand the tilde symbol in the file paths
path_lotus = os.path.expanduser(path_lotus)
path_lotus_smiles_out = os.path.expanduser(path_lotus_smiles_out)

# Read the data from path_lotus and select the "structure_smiles" column
lotus = polars.read_csv(path_lotus, columns=["structure_smiles"])

# Remove duplicate values in the "structure_smiles" column and rename it to "smiles"
lotus = lotus.unique("structure_smiles").rename({"structure_smiles": "smiles"})

# Write the resulting DataFrame to path_lotus_smiles_out as a CSV file
with gzip.open(path_lotus_smiles_out, "wt") as f:
    lotus.write_csv(f)
