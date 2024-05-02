"""Build a dictionary of class_id to class_structure from the provided classes list."""

from __future__ import annotations


def get_class_structures(classes: list[dict[str, list[str]]]) -> dict[str, list[str]]:
    """Build a dictionary of class_id to class_structure from the provided classes list."""
    class_structures: dict[str, list[str]] = {}
    for c in classes:
        for class_id, class_structure in c.items():
            class_structures[class_id] = class_structure
    return class_structures
