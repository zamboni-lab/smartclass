"""Get CHEBI."""

from __future__ import annotations

from smartclass.helpers.download_file_if_not_exists import download_file_if_not_exists

__all__ = ["get_chebi"]


def get_chebi_1(
    url: str = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo",
    output: str = "scratch/chebi.obo",
):
    """
    Get CHEBI 1.

    :param url: URL.
    :type url: str

    :param output: Output.
    :type output: str
    """
    download_file_if_not_exists(url=url, output=output)


def get_chebi_2(
    url: str = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_lite_3star.sdf.gz",
    output: str = "scratch/ChEBI_lite_3star.sdf.gz",
):
    """
    Get CHEBI 2.

    :param url: URL.
    :type url: str

    :param output: Output.
    :type output: str
    """
    download_file_if_not_exists(url=url, output=output)


def get_chebi():
    """Get CHEBI."""
    get_chebi_1()
    get_chebi_2()


if __name__ == "__main__":
    get_chebi()
