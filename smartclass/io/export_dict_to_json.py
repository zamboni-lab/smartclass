"""Export dict to json."""

from __future__ import annotations

import json

from smartclass.logging import get_logger

__all__ = [
    "export_dict_to_json",
]

logger = get_logger(__name__)


def export_dict_to_json(dic: dict, output: str) -> None:
    """Export dict to json.

    Parameters
    ----------
    dic : dict
        A dictionary.
    output : str
        Output path.
    """
    with open(output, "w", encoding="utf-8") as f:
        json.dump(dic, f, ensure_ascii=False, indent=4)
        logger.debug(f"Saved dict to {output}")
