"""Load tsv from path."""

from __future__ import annotations

from polars import DataFrame, read_csv


__all__ = [
    "load_tsv_from_path",
]


def load_tsv_from_path(path: str) -> DataFrame:
    """Load tsv from path.

Parameters
----------
path : str
    Path of the file.

Returns
-------
DataFrame
    DataFrame.
    """
    return read_csv(path, separator="\t")
