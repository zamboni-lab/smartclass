"""Load SMILES."""

from __future__ import annotations

import polars

__all__ = [
    "load_smiles",
]


def load_smiles(input: str, column: str = "smiles") -> list[str]:
    """
    Load SMILES.

    :param input: The path to the input CSV file.
    :type input: str

    :param column: The name of the column containing the SMILES.
    :type column: str

    :returns: List of SMILES.
    :rtype: list[str]
    """
    # Read the CSV file and extract the "smiles" column
    if input.endswith(".tsv"):
        df = polars.read_csv(input, separator="\t", columns=[column])
    else:
        # Read the CSV file with default separator
        df = polars.read_csv(input, columns=[column])

    # Filter empty values
    df = df.filter(df[column].is_not_null())

    # Uncomment the following lines for debug
    # df = df.head(1048576)
    # df = df.head(524288)
    # df = df.head(262144)
    # df = df.head(131072)
    # df = df.head(65536)
    # df = df.head(32768)
    # df = df.head(16384)
    # df = df.head(8192)
    # df = df.head(4096)
    # df = df.head(2048)
    # df = df.head(1024)
    # df = df.head(512)
    # df = df.head(256)
    # df = df.head(128)
    # df = df.head(64)
    # df = df.head(32)
    # df = df.head(16)
    # df = df.head(8)
    # df = df.head(4)
    # df = df.head(2)
    # df = df.tail(1)

    # Return the unique structure SMILES as a list
    return df[column].unique().to_list()
