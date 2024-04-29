"""Download and process the latest ChEMBL data."""

from __future__ import annotations

import logging
import pickle
import time

import chembl_downloader
from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary

from smartclass.resources.chembl.latest_chembl_paths import latest_chembl_paths

__all__ = [
    "get_latest_chembl",
]


def get_latest_chembl(
    fp_len: int = 2048,
    max_atoms: int = 50,
    report_interval: int = 50000,
    tautomer_fingerprints: bool = True,
) -> None:
    """
    Download and process the latest ChEMBL data.

    :param fp_len: Fingerprint length.
    :type fp_len: int

    :param max_atoms: Maximal number of atoms.
    :type max_atoms: int

    :param report_interval: Number of molecules to be processed before reporting.
    :type report_interval: int

    :param tautomer_fingerprints: Flag indicating whether tautomer fingerprints are used.
    :type tautomer_fingerprints: bool
    """
    try:
        # Get paths and constants
        version, data_path, lib_path = latest_chembl_paths()

        # Create data list
        data = []

        # Download ChEMBL data and process molecules
        with chembl_downloader.supplier(version=version) as suppl:
            t1 = time.time()
            for i, mol in enumerate(suppl):
                if not ((i + 1) % report_interval):
                    elapsed_time = time.time() - t1
                    logging.debug(f"Processed {i + 1} molecules in {elapsed_time:.1f} seconds")

                if mol is None or mol.GetNumAtoms() > max_atoms:
                    continue

                fp = Chem.PatternFingerprint(
                    mol, fpSize=fp_len, tautomerFingerprints=tautomer_fingerprints
                )
                smi = Chem.MolToSmiles(mol)
                data.append((smi, fp))

        # Save data to file
        with open(data_path, "wb") as file:
            pickle.dump(data, file)

        # Create and save the substructure library
        mols = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
        fps = rdSubstructLibrary.TautomerPatternHolder(fp_len)

        for smi, fp in data:
            mols.AddSmiles(smi)
            fps.AddFingerprint(fp)

        library = rdSubstructLibrary.SubstructLibrary(mols, fps)

        elapsed_time = time.time() - t1
        logging.debug(
            f"Processed ChEMBL data in {elapsed_time:.2f} seconds. The library has {len(library)} molecules."
        )

        with open(lib_path, "wb") as file:
            pickle.dump(library, file)
    except Exception as e:
        logging.error(f"Failed to process ChEMBL data: {e!s}")


if __name__ == "__main__":
    get_latest_chembl()
