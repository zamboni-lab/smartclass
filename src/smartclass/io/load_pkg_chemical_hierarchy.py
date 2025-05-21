"""Load package chemical classes hierarchy."""

from __future__ import annotations

__all__ = [
    "load_pkg_chemical_hierarchy",
]


def load_pkg_chemical_hierarchy(
    file_path: str = "scratch/wikidata_classes_taxonomy.tsv",
) -> dict:
    """
    Load package chemical classes hierarchy.

    :param file_path: The path to the file containing the hierarchy data.
    :type file_path: str

    :returns: A dictionary representing the parent-child relationships between classes.
    :rtype: dict
    """
    # First pass: Collect class URIs and parent URIs
    class_uris = set()
    parent_uris = set()

    with open(file_path) as tsv_file:
        # Skip the header line
        next(tsv_file)
        for line in tsv_file:
            parts = line.strip().split("\t")
            class_uri = parts[1]
            parent_uri = parts[2]

            # Collect class URIs and parent URIs
            class_uris.add(class_uri)
            parent_uris.add(parent_uri)

    # Create a dictionary with class URIs as keys
    class_hierarchy: dict = {class_uri: [] for class_uri in class_uris}

    # Second pass: Build the hierarchy
    with open(file_path) as tsv_file:
        # Skip the header line
        next(tsv_file)
        for line in tsv_file:
            parts = line.strip().split("\t")
            class_uri = parts[1]
            parent_uri = parts[2]

            # If the parent URI is not empty and exists in class_hierarchy, add it as a child
            if parent_uri and parent_uri in class_hierarchy:
                class_hierarchy[parent_uri].append(class_uri)

    return class_hierarchy

    return class_hierarchy


if __name__ == "__main__":
    class_hierarchy = load_pkg_chemical_hierarchy()
