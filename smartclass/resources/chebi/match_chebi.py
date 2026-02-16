"""Match ChEBI."""

from __future__ import annotations

import csv
import gzip
from zipfile import ZipFile

from polars import DataFrame
from rdkit.Chem import ForwardSDMolSupplier

from smartclass.chem.conversion.convert_mol_to_smiles import convert_mol_to_smiles
from smartclass.io.load_csv_from_path import load_csv_from_path  # noqa:F401


__all__ = [
    "match_chebi",
]


def match_chebi(
    ids_file_zip: str = "scratch/chebi_classyfire.csv.zip",
    ids_file: str = "ChEBI_126_classyfire_21_annotations.csv",
    chebi_sdf_file: str = "scratch/ChEBI_lite_3star.sdf.gz",
    output: str = "scratch/chebi_matched_molecules.tsv",
) -> None:
    """
    Match ChEBI.

    :param ids_file_zip: chebi_classyfire zip.
    :type ids_file_zip: str

    :param ids_file: chebi_classyfire.
    :type ids_file: str

    :param chebi_sdf_file: chebi sdf.
    :type chebi_sdf_file: str

    :param output: output.
    :type output: str
    """
    # Polars cannot read it properly directly
    with ZipFile(ids_file_zip, "r") as zip_ref:
        with zip_ref.open(ids_file) as csv_file:
            csv_content = csv_file.read().decode("utf-8")
    csv_reader = csv.DictReader(csv_content.splitlines())
    ids_df = (
        DataFrame(list(csv_reader))
        .rename({
            "Smiles": "ChemOntID",
            "ChemOntID": "ParentName",
            "ParentName": "drop",
        })
        .drop("drop")
    )

    with gzip.open(chebi_sdf_file, "rb") as sdf:
        supplier = ForwardSDMolSupplier(sdf)
        chebi_molecules = [
            (mol.GetProp("ChEBI ID").replace("CHEBI:", ""), convert_mol_to_smiles(mol))
            for mol in supplier
            if mol is not None
        ]

    matched_df = DataFrame(chebi_molecules, schema=["CompoundID", "smiles"])
    matched_df = matched_df.with_columns(
        matched_df["CompoundID"].cast(ids_df["CompoundID"].dtype),
    )
    merged_df = ids_df.join(matched_df, on="CompoundID", how="inner")

    merged_df.write_csv(output, separator="\t")


# Example usage
if __name__ == "__main__":
    match_chebi()
