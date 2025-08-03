"""Load a Polars DataFrame from a package file."""

from __future__ import annotations

import importlib_resources
from polars import DataFrame  # Because of type

from smartclass.io.load_csv_from_path import load_csv_from_path
from smartclass.io.load_json_from_path import load_json_from_path
from smartclass.io.load_tsv_from_path import load_tsv_from_path

__all__ = ["load_pkg_file"]


def load_pkg_file(file: str, directory: str = "smartclass.data") -> DataFrame:
    """
    Load data from a package file into a Polars DataFrame.

    :param file: The name of the file to load.
    :type file: str

    :param directory: The directory of the file to load. Default to `smartclass.data`
    :type directory: str

    :raises ValueError: If the loading of the file failed.

    :returns: A Polars DataFrame containing the loaded data.
    :rtype: DataFrame
    """
    try:
        ref = importlib_resources.files(directory) / file
        with importlib_resources.as_file(ref) as path:
            if ".tsv" in file:
                return load_tsv_from_path(path)
            elif ".csv" in file:
                return load_csv_from_path(path)
            elif ".json" in file:
                return DataFrame(load_json_from_path(path))

    except Exception as e:
        raise ValueError(f"Failed to load '{file}': {e!s}")


if __name__ == "__main__":
    df1 = load_pkg_file(file="chemontids__wdids.tsv")
    df2 = load_pkg_file(file="classes_smarts.tsv")
    df3 = load_pkg_file(file="mia_smarts.tsv")
