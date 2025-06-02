from __future__ import annotations

import polars as pl
from tqdm import tqdm

from smartclass.chem.helpers.canonicalize_smarts import canonicalize_smarts

smarts_path = "scratch/wikidata_classes_smarts.tsv"
output_path = "scratch/wikidata_classes_smarts_canonicalized.tsv"

# Read the input
smarts_df = pl.read_csv(smarts_path, separator="\t")

# Ensure only distinct SMARTS from 'structure'
distinct_smarts_df = smarts_df.select("structure").unique()

# Convert to list for tqdm processing
structures = distinct_smarts_df["structure"].to_list()

# Run canonicalization with tqdm progress bar
smarts_canonical = []
smarts_canonical_with_mapping = []

for s in tqdm(structures, desc="Canonicalizing SMARTS"):
    smarts_canonical.append(canonicalize_smarts(s, return_mapping=False))
    smarts_canonical_with_mapping.append(canonicalize_smarts(s, return_mapping=True))

# Reconstruct DataFrame with results
result_df = pl.DataFrame(
    {
        "structure": structures,
        "smarts_canonical": smarts_canonical,
        "smarts_canonical_with_mapping": smarts_canonical_with_mapping,
    },
    strict=False,
)

# Write to file
result_df.write_csv(output_path, separator="\t")
