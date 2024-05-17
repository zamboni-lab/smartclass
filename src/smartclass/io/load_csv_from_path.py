"""Load csv from path."""

from __future__ import annotations

from polars import read_csv

__all__ = [
    "load_csv_from_path",
]


def load_csv_from_path(path: str) -> dict:
    """
    Load csv from path.

    :param path: Path of the file.
    :type path: str

    :returns: A dictionary.
    :rtype: dict
    """
    return read_csv(path)
