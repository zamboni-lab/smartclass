"""TODO.

TODO.
"""

from __future__ import annotations

import json
import logging
import os

import requests


def download_file(url: str, file_name: str) -> None:
    """
    Description.

    :param url: Description.
    :type url: str

    :param file_name: Description.
    :type file_name: str
    """
    file_path = os.path.join("data", file_name)  # Construct the full file path
    if os.path.exists(file_path):
        logging.debug(f"File {file_path} already exists. Skipping download.")

    response = requests.get(url, timeout=60)
    if response.status_code == 200:
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, "wb") as f:
            f.write(response.content)
        logging.debug(f"Downloaded {url} and saved as {file_path}")
    else:
        logging.debug(f"Failed to download {url}")


def load_json_from_url(url: str, name: str):
    """
    Description.

    :param url: Description.
    :type url: str

    :param name:Description.
    :type name: str

    :returns: Description.
    :rtype: Type
    """
    file_name = f"{name}.json"
    file_path = os.path.join("data", file_name)  # Construct the full file path
    if os.path.exists(file_path):
        logging.debug(f"JSON file {file_path} already exists. Loading from the file.")
        with open(file_path, encoding="utf-8") as f:
            data = json.load(f)
        return data

    response = requests.get(url, timeout=60)
    if response.status_code == 200:
        data = json.loads(response.text)
        logging.debug(f"Downloaded and loaded {name}_json successfully.")
        # Save the JSON data to a file
        with open(file_path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=4)
        logging.debug(f"Saved {name}_json to {file_path}")
        return data
    else:
        logging.debug(f"Failed to download {name}_json.")
        return None


# MISSING DEF
# __all__ = [
#     "get_existing_classes",
# ]

# Define the URLs for downloading JSON files
json_urls = {
    "classyfire": "http://classyfire.wishartlab.com/tax_nodes.json",
    "npclassifier": "https://raw.githubusercontent.com/mwang87/NP-Classifier/master/Classifier/dict/index_v1.json",
}

# Load JSON data from URLs or from local files if they exist
for name, url in json_urls.items():
    data = load_json_from_url(url, name)
    if data is not None:
        # Process the loaded JSON data as needed
        pass
