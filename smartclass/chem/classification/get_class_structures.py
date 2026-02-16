"""Build class_id to class_structure mapping from a classes list."""

from __future__ import annotations


def get_class_structures(classes: list[dict[str, list[str]]]) -> dict[str, list[str]]:
    """Build a mapping from class IDs to their SMARTS structures."""
    class_structures: dict[str, list[str]] = {}
    for c in classes:
        for class_id, class_structure in c.items():
            class_structures[class_id] = class_structure
    return class_structures
