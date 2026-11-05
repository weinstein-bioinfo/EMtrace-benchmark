"""Adapters that load pipeline outputs into the canonical AbundanceDF schema.

Each adapter exposes a top-level `load(path)` function that returns a
pandas.DataFrame with at least the columns in `schema.REQUIRED_COLS`.
"""

from __future__ import annotations

from typing import Callable
import pandas as pd

from . import czid, epi2me, kraken_bracken, megan_lr, metamix, emtrace


_ADAPTERS: dict[str, Callable[[str], pd.DataFrame]] = {
    "emtrace": emtrace.load,
    "czid": czid.load,
    "epi2me": epi2me.load,
    "kraken_bracken": kraken_bracken.load,
    "megan_lr": megan_lr.load,
    "metamix": metamix.load,
}


def get_adapter(name: str) -> Callable[[str], pd.DataFrame]:
    """Return the load() function for the named adapter."""
    if name in _ADAPTERS:
        return _ADAPTERS[name]
    raise ValueError(f"Unknown adapter '{name}'. Available: {sorted(_ADAPTERS)}")


def available() -> list[str]:
    """Names of adapters available in this build."""
    return sorted(_ADAPTERS)
