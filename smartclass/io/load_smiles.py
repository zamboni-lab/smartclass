"""Load SMILES strings from CSV/TSV files."""

from __future__ import annotations

from pathlib import Path

import polars

from smartclass.exceptions import DataLoadingError
from smartclass.logging import get_logger


__all__ = [
    "load_smiles",
]

logger = get_logger(__name__)


def load_smiles(
    input: str | Path,
    column: str = "smiles",
    limit: int | None = None,
) -> list[str]:
    """Load unique SMILES strings from a CSV or TSV file.

Parameters
----------
input : str | Path
    Path to the input file (CSV or TSV).
column : str
    Default is 'smiles'.
limit : int | None
    SMILES to return (for testing). Default is None.

Returns
-------
list[str]
    SMILES strings.

Raises
------
    DataLoadingError
        If the file cannot be read or column not found.
    """
    input_path = Path(input)

    if not input_path.exists():
        raise DataLoadingError(str(input_path), reason="File not found")

    try:
        # Determine separator from file extension
        separator = "\t" if input_path.suffix == ".tsv" else ","

        df = polars.read_csv(
            input_path,
            separator=separator,
            columns=[column],
        )

        # Filter out null/empty values
        df = df.filter(df[column].is_not_null())

        # Apply limit if specified (useful for testing)
        if limit is not None and limit > 0:
            df = df.head(limit)

        unique_smiles = df[column].unique().to_list()
        logger.debug(f"Loaded {len(unique_smiles)} unique SMILES from {input_path}")
        return unique_smiles

    except polars.exceptions.ColumnNotFoundError as e:
        raise DataLoadingError(
            str(input_path),
            reason=f"Column '{column}' not found in file",
        ) from e
    except Exception as e:
        raise DataLoadingError(str(input_path), reason=str(e)) from e
