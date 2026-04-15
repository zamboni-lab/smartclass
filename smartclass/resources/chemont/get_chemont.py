"""Get ChemOnt."""

from __future__ import annotations

from smartclass.helpers.download_file_if_not_exists import download_file_if_not_exists


__all__ = ["get_chemont"]


def get_chemont(
    url: str = "http://classyfire.wishartlab.com/system/downloads/1_0/chemont/ChemOnt_2_1.obo.zip",
    output: str = "scratch/ChemOnt_2_1.obo.zip",
):
    """Get ChemOnt.

    Parameters
    ----------
    url : str
        Default is 'http://classyfire.wishartlab.com/system/downloads/1_0/chemont/ChemOnt_2_1.obo.zip'.
    output : str
        Default is 'scratch/ChemOnt_2_1.obo.zip'.
    """
    download_file_if_not_exists(url=url, output=output)


if __name__ == "__main__":
    get_chemont()
