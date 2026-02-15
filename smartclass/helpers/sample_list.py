"""Utility for random sampling from collections."""

from __future__ import annotations

from random import sample
from typing import TypeVar

__all__ = ["sample_list"]

T = TypeVar("T")


def sample_list(items: list[T], max_samples: int = 1000) -> list[T]:
    """
    Randomly sample items from a list.

    Returns at most `max_samples` items. If the list has fewer items
    than `max_samples`, returns all items in random order.

    :param items: List of items to sample from.
    :param max_samples: Maximum number of samples to return.
    :returns: List of randomly sampled items.
    """
    sample_size = min(len(items), max_samples)
    return sample(items, sample_size)
