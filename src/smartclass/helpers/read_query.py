"""Reads a sparql query."""

from __future__ import annotations

import requests

__all__ = [
    "read_query",
]


def read_query(query):
    if query.startswith("http://") or query.startswith("https://"):
        response = requests.get(query)
        response.raise_for_status()
        return response.text
    else:
        with open(query, "r", encoding="utf-8") as file:
            return file.read()
