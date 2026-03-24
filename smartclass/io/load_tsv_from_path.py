"""Load tsv from path."""

from __future__ import annotations

from polars import DataFrame, read_csv


__all__ = [
    "load_tsv_from_path",
]


def load_tsv_from_path(path: str) -> DataFrame:
    """
    Load tsv from path.

    :param path: Path of the file.
    :type path: str

    :returns: A Polars DataFrame.
    :rtype: DataFrame
    """
    return read_csv(path, separator="\t")
