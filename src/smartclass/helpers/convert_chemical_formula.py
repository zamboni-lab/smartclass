"""Convert chemical formula."""

from __future__ import annotations

import re

__all__ = [
    "convert_chemical_formula",
]


def convert_chemical_formula(text: str) -> str:
    """
    Convert chemical formula.

    :param text: Text.
    :type text: str

    :returns: Modified text.
    :rtype: str
    """

    def replace(match: re.Match) -> str:
        """
        Matches subscripts.

        :param match: Match object.
        :type match: re.Match

        :returns: Modified text.
        :rtype: str
        """
        number, symbol = match.groups()
        subscript_map = {
            "0": "₀",
            "1": "₁",
            "2": "₂",
            "3": "₃",
            "4": "₄",
            "5": "₅",
            "6": "₆",
            "7": "₇",
            "8": "₈",
            "9": "₉",
        }
        superscript_map = {
            "0": "⁰",
            "1": "¹",
            "2": "²",
            "3": "³",
            "4": "⁴",
            "5": "⁵",
            "6": "⁶",
            "7": "⁷",
            "8": "⁸",
            "9": "⁹",
            "+": "⁺",
            "-": "⁻",
        }

        if number:
            # Check if the number is preceded by `+` or `-`
            prev_char = text[match.start(1) - 1 : match.start(1)]
            if prev_char in "+-":
                return "".join(superscript_map[char] for char in number)
            else:
                return "".join(subscript_map[char] for char in number)
        elif symbol:
            return "".join(superscript_map[char] for char in symbol)
        return ""

    pattern = re.compile(r"(\d+)|([+-])")
    result = re.sub(pattern, replace, text)
    # move `⁺` and `⁻` to the after the number (Wikidata formatting)
    result = re.sub(r"(\⁺|\⁻)+([⁰¹²³⁴⁵⁶⁷⁸⁹]+)", r"\2\1", result)

    return result


if __name__ == "__main__":
    formulas_to_test = ["C19H36Cl2CrN3O-3", "C29H32Cl2N6O17P4-12"]

    for formula in formulas_to_test:
        convert_chemical_formula(formula)
