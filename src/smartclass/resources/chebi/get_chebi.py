"""Get CHEBI."""

from __future__ import annotations

from smartclass.helpers.download_file_if_not_exists import download_file_if_not_exists


def get_chebi(
    url: str = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo",
    output: str = "scratch/chebi.obo",
):
    """
    Get CHEBI.

    :param url: URL.
    :type url: str

    :param output: Output.
    :type output: str
    """
    download_file_if_not_exists(url=url, output=output)


if __name__ == "__main__":
    get_chebi()
