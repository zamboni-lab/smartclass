"""Load a Polars DataFrame from an external tsv file with chemical classes."""

from __future__ import annotations

import polars
from polars import DataFrame  # Because of type

__all__ = ["load_external_classes_file"]


def load_external_classes_file(
    file: str, id_name: None | str = "class", smarts_name: None | str = "smarts"
) -> DataFrame:
    """
    Load a Polars DataFrame from an external tsv file with chemical classes.

    :param file: The name of the file to load.
    :type file: str

    :param id_name: The name of the column containing ids.
    :type id_name: str

    :param smarts_name: The name of the column containing SMARTS.
    :type smarts_name: str

    :returns: A Polars DataFrame containing the loaded data.
    :rtype: DataFrame
    """
    return polars.read_csv(file, truncate_ragged_lines=True, separator="\t").select(
        {
            id_name: "class_id",
            smarts_name: "class_smarts",
        }
    )


if __name__ == "__main__":
    df = load_external_classes_file(file="scratch/wikidata_classes_smarts.tsv")
