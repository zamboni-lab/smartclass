"""Export dict to json."""

from __future__ import annotations

import json
import logging


__all__ = [
    "export_dict_to_json",
]


def export_dict_to_json(dic: dict, output: str) -> None:
    """
    Export dict to json.

    :param dic: A dictionary.
    :type dic: dict

    :param output: Output path.
    :type output: str
    """
    with open(output, "w", encoding="utf-8") as f:
        json.dump(dic, f, ensure_ascii=False, indent=4)
        logging.debug(f"Saved dict to {output}")
