"""Load bitter SMILES data from the package file into a Polars DataFrame."""

from __future__ import annotations

from typing import TYPE_CHECKING

from smartclass.io.load_pkg_file import load_pkg_file

if TYPE_CHECKING:
    from polars import DataFrame

__all__ = ["load_pkg_classes"]


def load_pkg_bitter_smiles() -> DataFrame:
    """
    Load bitter SMILES data from the package file into a Polars DataFrame.

    :returns: DataFrame containing bitter SMILES.
    :rtype: DataFrame
    """
    # Was obtained using
    # SELECT ?smiles WHERE {
    #   ?structure wdt:P2017 ?smiles;
    #     wdt:P1552 wd:Q1517187.
    # }
    # LIMIT 100
    return load_pkg_file(file="bitter_smiles.tsv")


if __name__ == "__main__":
    df = load_pkg_classes()
