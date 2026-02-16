"""Export a list of dictionaries to a CSV or TSV file."""

from __future__ import annotations

import csv
from pathlib import Path

from smartclass.exceptions import DataExportError
from smartclass.logging import get_logger


__all__ = [
    "export_results",
]

logger = get_logger(__name__)


def export_results(output: str | Path, results: list[dict]) -> None:
    """
    Export a list of dictionaries to a CSV or TSV file.

    The output format is determined by the file extension:
    - .tsv: Tab-separated values
    - .csv: Comma-separated values

    :param output: Path to the output file.
    :param results: List of dictionaries to export as rows.
    :raises DataExportError: If the export fails.
    """
    if not results:
        logger.debug("No data to export.")
        return

    output_path = Path(output)

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Determine delimiter from file extension
    delimiter = "\t" if output_path.suffix == ".tsv" else ","

    fieldnames = list(results[0].keys())

    try:
        with open(output_path, "w", newline="", encoding="utf-8") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(results)
        logger.debug(f"Data exported to {output_path}")
    except OSError as e:
        raise DataExportError(str(output_path), reason=str(e)) from e


if __name__ == "__main__":
    # Example data to export
    example_results = [
        {"id": 1, "name": "Alice", "age": 30},
        {"id": 2, "name": "Bob", "age": 25},
        {"id": 3, "name": "Charlie", "age": 35},
    ]

    # Specify the path to the output CSV file
    output_file = "output.csv"

    # Call the export_results function with the example data and file path
    export_results(output_file, example_results)
