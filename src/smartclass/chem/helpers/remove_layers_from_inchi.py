"""Remove layers from InChI."""

from __future__ import annotations

import re

__all__ = [
    "remove_layers_from_inchi",
]


def remove_layers_from_inchi(inchi: str, layers: list) -> str:
    """
    Remove layers from InChI.

    :param inchi: InChI.
    :type inchi: str

    :param layers: Layer(s) to remove should be in [c, h, q, p, b, t, m, s, i].
    :type layers: list

    :returns: InChI.
    :rtype: str
    """
    for layer in layers:
        inchi = re.sub(rf"/{layer}.*?/", "/", inchi)
    return inchi


if __name__ == "__main__":
    inchis_to_test = [
        "InChI=1S/C16H22O9/c1-2-7-8-3-4-22-14(21)9(8)6-23-15(7)25-16-13(20)12(19)11(18)10(5-17)24-16/h2,6-8,10-13,15-20H,1,3-5H2/t7-,8+,10-,11-,12+,13-,15+,16+/m1/s1",
    ]
    layers_to_remove = ["b", "t", "m", "s"]

    for inchi in inchis_to_test:
        remove_layers_from_inchi(inchi, layers_to_remove)
