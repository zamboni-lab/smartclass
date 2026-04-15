"""Load a Polars DataFrame from an external tsv file with chemical classes."""

from __future__ import annotations

import polars
from polars import DataFrame  # Because of type


__all__ = ["load_external_classes_file"]


def load_external_classes_file(
    file: str,
    id_name: str = "class",
    smarts_name: str = "structure",
) -> DataFrame:
    """Load a Polars DataFrame from an external tsv file with chemical classes.

    Parameters
    ----------
    file : str
        The name of the file to load.
    id_name : str
        Default is 'class'.
    smarts_name : str
        Default is 'structure'.

    Returns
    -------
    DataFrame
        DataFrame containing the loaded data.
    """
    df = polars.read_csv(file, truncate_ragged_lines=True, separator="\t")
    return df.select([
        polars.col(id_name).alias("class_id"),
        polars.col(smarts_name).alias("class_smarts"),
    ])


if __name__ == "__main__":
    df = load_external_classes_file(file="scratch/wikidata_classes_smarts.tsv")
