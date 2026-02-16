"""Split a CSV file into multiple smaller CSV files."""

from __future__ import annotations

import csv
import os


__all__ = [
    "split_csv",
]


def split_csv(input_file: str, output_dir: str, lines_per_file: int = 5000) -> None:
    """
    Split a CSV file into multiple smaller CSV files.

    :param input_file: Path to the input CSV file.
    :type input_file: str

    :param output_dir: Directory to save the split CSV files.
    :type output_dir: str

    :param lines_per_file: Number of lines per output file.
    :type lines_per_file: int
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    with open(input_file, newline="") as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)  # Read and save the header

        file_count = 1
        line_count = 0
        output_file = None

        for row in reader:
            if line_count % lines_per_file == 0:
                # Close the current file if it's open
                if output_file:
                    output_file.close()
                output_file_path = os.path.join(output_dir, f"split_{file_count}.csv")
                file_count += 1
                output_file = open(output_file_path, "w", newline="")
                writer = csv.writer(output_file)
                writer.writerow(header)
            writer.writerow(row)
            line_count += 1

        # Close the last output file
        if output_file:
            output_file.close()


# Example usage:
if __name__ == "__main__":
    # Replace with your input CSV file
    input_file = "scratch/maintenance_improve_subclasses_inchikeys_2.csv"
    output_directory = "scratch/output"  # Replace with your output directory
    lines_per_file = 5000  # Replace with the number of lines per output file

    split_csv(input_file, output_directory, lines_per_file)
