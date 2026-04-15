"""Load json from path."""

from __future__ import annotations

import json


__all__ = [
    "load_json_from_path",
]


def load_json_from_path(path: str) -> dict:
    """Load json from path.

    Parameters
    ----------
    path : str
        Path of the file.

    Returns
    -------
    dict
        A dictionary.
    """
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return data
