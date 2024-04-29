"""TODO.

TODO.
"""

# MISSING DEF
# __all__ = [
#     "match_chebi",
# ]
from __future__ import annotations

import polars
from rdkit import Chem

ids_file = "../../Downloads/ChEBI_126_classyfire_21_annotations.csv"
ids_df = polars.read_csv(ids_file)

chebi_sdf_file = "../../Downloads/ChEBI_lite_3star.sdf"
chebi_molecules = list(Chem.MultithreadedSDMolSupplier(chebi_sdf_file))

matched_molecules = []

for mol in chebi_molecules:
    if mol is not None:
        chebi_id = mol.GetProp("ChEBI ID")
        chebi_id = chebi_id.replace("CHEBI:", "")
        smiles = Chem.MolToSmiles(mol)
        matched_molecules.append((chebi_id, smiles))

matched_df = polars.DataFrame(
    {
        "CompoundID": [row[0] for row in matched_molecules],
        "SMILES": [row[1] for row in matched_molecules],
    }
)
matched_df = matched_df.with_columns(matched_df["CompoundID"].cast(ids_df["CompoundID"].dtype))
merged_df = ids_df.join(matched_df, on="CompoundID", how="inner")

output_csv_file = "matched_molecules.csv"
merged_df.write_csv(output_csv_file)
