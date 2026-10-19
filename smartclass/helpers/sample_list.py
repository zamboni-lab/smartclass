"""Utility for random sampling from collections."""

from __future__ import annotations

from random import sample

__all__ = ["sample_list"]


def sample_list[T](items: list[T], max_samples: int = 1000) -> list[T]:
    """Randomly sample items from a list.

    Returns at most ``max_samples`` items. If the list has fewer items
    than ``max_samples``, returns all items in random order.

    Parameters
    ----------
    items : list[T]
        List of items to sample from.
    max_samples : int
        Maximum number of items to return. Default is 1000.

    Returns
    -------
    list[T]
        List of randomly sampled items.
    """
    sample_size = min(len(items), max_samples)
    return sample(items, sample_size)
