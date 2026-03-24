"""Load csv from path."""

from __future__ import annotations

from polars import DataFrame, read_csv


__all__ = [
    "load_csv_from_path",
]


def load_csv_from_path(path: str) -> DataFrame:
    """
    Load csv from path.

    :param path: Path of the file.
    :type path: str

    :returns: A Polars DataFrame.
    :rtype: DataFrame
    """
    return read_csv(path)
