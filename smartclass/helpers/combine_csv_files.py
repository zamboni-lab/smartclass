"""Combine multiple CSV files into a single CSV file."""

from __future__ import annotations

import polars

from smartclass.logging import get_logger

__all__ = [
    "combine_csv_files",
]

logger = get_logger(__name__)


def combine_csv_files(
    input_files: list[str] | str,
    output_file: str,
    separator: str = "\t",
) -> None:
    """Combine multiple CSV files into a single CSV file.

    Parameters
    ----------
    input_files : list[str] | str
        One or more paths to delimited files to combine.
    output_file : str
        Path where the combined data will be saved.
    separator : str
        Field separator used for both reading and writing (default is a tab).
    """
    if isinstance(input_files, str):
        input_files = [input_files]  # Ensure input_files is a list

    if not input_files:
        raise ValueError("input_files must contain at least one file path")

    # Create an empty list to store DataFrames
    dfs = []

    # Read each CSV file and store it in the list
    for file_path in input_files:
        df = polars.read_csv(file_path, separator=separator)
        dfs.append(df)

    # Create a merged DataFrame
    concatenated_df = dfs[0]

    # Iterate over the remaining DataFrames and concatenate them row-wise
    for df in dfs[1:]:
        concatenated_df = concatenated_df.vstack(df)

    # Export the merged DataFrame as a CSV file
    concatenated_df.write_csv(output_file, separator=separator)


if __name__ == "__main__":
    input_files = [
        "scratch/wikidata_classes_cxsmiles.tsv",
        "scratch/wikidata_classes_smarts.tsv",
        "scratch/wikidata_classes_smiles_canonical.tsv",
        "scratch/wikidata_classes_smiles_isomeric.tsv",
    ]
    output_file = "scratch/wikidata_classes_full.tsv"

    combine_csv_files(input_files, output_file)
    logger.info(output_file)
