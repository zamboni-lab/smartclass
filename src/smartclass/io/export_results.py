"""Export a list of dictionaries to a CSV file."""

from __future__ import annotations

import csv
import logging

__all__ = [
    "export_results",
]


def export_results(output: str, results: list[dict]) -> None:
    """
    Export a list of dictionaries to a CSV file.

    :param output: The path to the output CSV file.
    :type output: str

    :param results: A list of dictionaries to be exported as CSV rows.
    :type results: list[dict]
    """
    if not results:
        logging.debug("No data to export.")
        return

    fieldnames = results[0].keys()

    try:
        with open(output, "w", newline="") as file:
            if output.endswith(".tsv"):
                delimiter = "\t"
            elif output.endswith(".csv"):
                delimiter = ","
            writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter=delimiter)
            writer.writeheader()
            writer.writerows(results)
        logging.debug(f"Data exported to {output}")
    except Exception as e:
        logging.error(f"Error exporting data to {output}: {e!s}")


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
