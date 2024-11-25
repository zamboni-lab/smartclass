"""Get existing classified datasets."""

from __future__ import annotations

import os

from smartclass.helpers.download_file_if_not_exists import download_file_if_not_exists
from smartclass.io.load_pkg_file import load_pkg_file

__all__ = [
    "get_existing_classified_datasets",
]


def get_existing_classified_datasets(
    existing_classsified_datasets: dict | None = None,
) -> None:
    """
    Get existing classified datasets.

    :param existing_classsified_datasets: Fixed dictionary of existing classified datasets.
    :type existing_classsified_datasets: Union[dict, None]
    """
    if existing_classsified_datasets is None:
        existing_classsified_datasets = load_pkg_file("existing_classified_datasets.json").to_dict(
            as_series=False
        )

    for name, url in existing_classsified_datasets.items():
        download_file_if_not_exists(url[0], os.path.join("scratch", name))


# Example usage:
if __name__ == "__main__":
    get_existing_classified_datasets()
