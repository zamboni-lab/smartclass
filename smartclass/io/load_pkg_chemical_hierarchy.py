"""Load package chemical classes hierarchy."""

from __future__ import annotations

from pathlib import Path

from smartclass.exceptions import DataLoadingError
from smartclass.logging import get_logger

__all__ = [
    "load_pkg_chemical_hierarchy",
]

logger = get_logger(__name__)


def load_pkg_chemical_hierarchy(
    file_path: str | Path = "scratch/wikidata_classes_taxonomy.tsv",
) -> dict[str, list[str]]:
    """
    Load chemical class hierarchy from a TSV file.

    The file should have at least 3 columns, where:
    - Column 2 (index 1): Class URI
    - Column 3 (index 2): Parent URI

    :param file_path: Path to the TSV file containing hierarchy data.
    :returns: Dictionary mapping parent class URIs to lists of child class URIs.
    :raises DataLoadingError: If the file cannot be read or parsed.
    """
    file_path = Path(file_path)

    if not file_path.exists():
        raise DataLoadingError(str(file_path), reason="File not found")

    try:
        # First pass: Collect all class URIs
        class_uris: set[str] = set()

        with open(file_path, encoding="utf-8") as tsv_file:
            next(tsv_file)  # Skip header
            for line in tsv_file:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    class_uris.add(parts[1])

        # Initialize hierarchy with empty lists
        class_hierarchy: dict[str, list[str]] = {uri: [] for uri in class_uris}

        # Second pass: Build parent-child relationships
        with open(file_path, encoding="utf-8") as tsv_file:
            next(tsv_file)  # Skip header
            for line in tsv_file:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    class_uri = parts[1]
                    parent_uri = parts[2]

                    # Add class as child of parent if parent exists
                    if parent_uri and parent_uri in class_hierarchy:
                        class_hierarchy[parent_uri].append(class_uri)

        logger.debug(f"Loaded {len(class_hierarchy)} classes from hierarchy file")
        return class_hierarchy

    except OSError as e:
        raise DataLoadingError(str(file_path), reason=str(e)) from e


if __name__ == "__main__":
    hierarchy = load_pkg_chemical_hierarchy()
