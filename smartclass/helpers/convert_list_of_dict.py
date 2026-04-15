"""Convert a list of dictionaries to a dictionary with possible inversion."""

from __future__ import annotations


__all__ = [
    "convert_list_of_dict",
]


def convert_list_of_dict(
    list_of_dict: list,
    key: str,
    value: str,
    invert: bool = False,
) -> dict:
    """Convert a list of dictionaries to a dictionary with possible inversion.

    Parameters
    ----------
    list_of_dict : list
        The list of dictionaries to convert.
    key : str
        The key.
    value : str
        The value.
    invert : bool
        False. Default is False.

    Returns
    -------
    dict
        A dictionary with given keys and values.
    """
    result_dict: dict = {}
    for item in list_of_dict:
        item_key = item[key]
        item_value = item[value]
        if invert:
            result_dict.setdefault(item_value, []).append(item_key)
        else:
            result_dict.setdefault(item_key, []).append(item_value)
    return result_dict


# Example usage:
if __name__ == "__main__":
    data = [
        {"name": "Alice", "age": 30},
        {"name": "Bob", "age": 25},
        {"name": "Charlie", "age": 35},
    ]

    result = convert_list_of_dict(data, key="name", value="age")
