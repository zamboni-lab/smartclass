"""Load the latest ChEMBL data and return a SubstructLibrary."""

from __future__ import annotations

import logging
import os
import pickle
from typing import TYPE_CHECKING

from smartclass.resources.chembl.get_latest_chembl import get_latest_chembl
from smartclass.resources.chembl.latest_chembl_paths import latest_chembl_paths


if TYPE_CHECKING:
    from rdkit.Chem import rdSubstructLibrary

__all__ = [
    "load_latest_chembl",
]


def load_latest_chembl() -> rdSubstructLibrary.SubstructLibrary:
    """
    Load the latest ChEMBL data and return a SubstructLibrary.

    :raises Exception: If the download of the file failed.

    :returns: SubstructLibrary containing ChEMBL data.
    :rtype: rdSubstructLibrary.SubstructLibrary
    """
    # Get the latest ChEMBL version and paths
    version, data_path, lib_path = latest_chembl_paths()

    try:
        # Attempt to open and load the SubstructLibrary from the library path
        with open(lib_path, "rb") as inf:
            sslib = pickle.load(inf)
        logging.info(f"SubstructLibrary loaded with {len(sslib)} molecules")
    except (FileNotFoundError, pickle.PickleError) as e:
        # If loading fails or library doesn't exist, download latest ChEMBL
        if os.path.exists(lib_path):
            os.remove(lib_path)  # Delete the corrupted library file

        logging.exception(f"Error loading library: {e!s}")
        logging.info("Downloading the latest ChEMBL data...")

        try:
            get_latest_chembl()  # Try to download the latest ChEMBL data
            with open(lib_path, "rb") as inf:
                sslib = pickle.load(inf)
            logging.info(f"SubstructLibrary loaded with {len(sslib)} molecules")
        except Exception as download_error:
            logging.exception(
                f"Failed to download or load ChEMBL data: {download_error!s}",
            )
            raise download_error  # Raise the error to notify the caller

    return sslib


if __name__ == "__main__":
    sslib = load_latest_chembl()
