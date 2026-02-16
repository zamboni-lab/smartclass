"""Get existing classes."""

from __future__ import annotations

from smartclass.io.load_json_from_url_or_path import load_json_from_url_or_path
from smartclass.io.load_pkg_file import load_pkg_file


__all__ = [
    "get_existing_classes",
]


def get_existing_classes(existing_classes: dict | None = None) -> None:
    """
    Get existing classes.

    :param existing_classes: Fixed dictionary of existing classes.
    :type existing_classes: Union[dict, None]
    """
    if existing_classes is None:
        existing_classes = load_pkg_file("existing_classes.json").to_dict(
            as_series=False,
        )

    for name, url in existing_classes.items():
        load_json_from_url_or_path(url[0], name)


# Example usage:
if __name__ == "__main__":
    get_existing_classes()
