"""Load a Polars DataFrame from the package file."""

from __future__ import annotations

from typing import TYPE_CHECKING

from smartclass.io.load_pkg_file import load_pkg_file


if TYPE_CHECKING:
    from polars import DataFrame

__all__ = ["load_pkg_mia"]


def load_pkg_mia() -> DataFrame:
    """
    Load Mono Indole Alkaloids (MIA) data from the package file into a Polars DataFrame.

    :returns: DataFrame containing MIA data.
    :rtype: DataFrame
    """
    return load_pkg_file(file="mia_smarts.tsv")


if __name__ == "__main__":
    df = load_pkg_mia()
