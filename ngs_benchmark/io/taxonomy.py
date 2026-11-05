"""Taxonomy enrichment helper for adapters that emit taxid-only output.

Wraps `emtrace.models.taxonomy.NCBITaxonomy` (the dump parser already used
by the parent pipeline). The wrapper is intentionally thin:

- `lineage(taxid)` returns the canonical schema's rank columns
  (species, genus, family, order, class, division) for one taxid.
- `enrich_with_taxonomy(df, taxid_col)` fills those columns in-place on a
  pandas DataFrame, replacing any pre-existing values.

The taxonomy DB is located via:

1. `path` kwarg (highest priority)
2. `NGS_BMARK_NCBI_TAXONOMY_DB` env var
3. `$EMTRACE_DB_DIR/taxonomy/ncbi/taxonomy.db` (default `~/databases/emtrace`),
   the EMtrace database directory

If `emtrace` is not on the Python path, we attempt to add `$EMTRACE_SRC`
(default: a sibling `../EMtrace/src` checkout). If that also fails, a clear ImportError tells the
user how to fix it.
"""

from __future__ import annotations

import functools
import logging
import os
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

from ngs_benchmark.schema import TAXRANKS, UNKNOWN

logger = logging.getLogger(__name__)


DEFAULT_DB_PATH = (
    Path(os.environ.get("EMTRACE_DB_DIR", "~/databases/emtrace")).expanduser()
    / "taxonomy" / "ncbi" / "taxonomy.db"
)
DEFAULT_EMTRACE_SRC = Path(os.environ.get("EMTRACE_SRC", "../EMtrace/src"))


def _import_ncbi_taxonomy():
    """Import emtrace.models.taxonomy.NCBITaxonomy with a sibling-checkout fallback."""
    try:
        from emtrace.models.taxonomy import NCBITaxonomy  # type: ignore[import-not-found]
        return NCBITaxonomy
    except ImportError:
        pass

    if DEFAULT_EMTRACE_SRC.is_dir() and str(DEFAULT_EMTRACE_SRC) not in sys.path:
        sys.path.insert(0, str(DEFAULT_EMTRACE_SRC))
        try:
            from emtrace.models.taxonomy import NCBITaxonomy  # type: ignore[import-not-found]
            return NCBITaxonomy
        except ImportError:
            pass

    raise ImportError(
        "ngs-benchmark needs the `emtrace` package for NCBI taxonomy lookups "
        "in comparator adapters (megan_lr, metamix, epi2me, kraken_bracken). "
        "Install it via `pip install -e /path/to/EMtrace/` or set "
        "PYTHONPATH=/path/to/EMtrace/src."
    )


def _resolve_db_path(path: Optional[Path | str]) -> Path:
    if path is not None:
        return Path(path)
    env = os.environ.get("NGS_BMARK_NCBI_TAXONOMY_DB")
    if env:
        return Path(env)
    return DEFAULT_DB_PATH


@functools.lru_cache(maxsize=4)
def get_taxonomy(path: Optional[Path | str] = None):
    """Return a cached NCBITaxonomy instance pointing at the resolved DB path.

    Cache keyed on the resolved path string so different DBs can coexist within
    one Python process (e.g. virus-only vs. all-taxa builds).
    """
    db_path = _resolve_db_path(path)
    if not db_path.exists():
        raise FileNotFoundError(
            f"NCBI taxonomy DB not found at {db_path}. "
            f"Set NGS_BMARK_NCBI_TAXONOMY_DB or pass --tax-db."
        )
    NCBITaxonomy = _import_ncbi_taxonomy()
    return NCBITaxonomy(db_path)


def lineage(taxid: int, path: Optional[Path | str] = None) -> dict[str, str]:
    """Return the canonical schema's rank values for one taxid.

    Missing ranks are reported as the literal "unknown" so downstream code
    treats them consistently with §3.1's NA-handling rule.
    """
    tax = get_taxonomy(path)
    lin = tax.get_taxonomy(int(taxid)) or {}
    division = tax.get_division(int(taxid)) or UNKNOWN
    return {
        "species": (lin.get("species") or UNKNOWN),
        "genus": (lin.get("genus") or UNKNOWN),
        "family": (lin.get("family") or UNKNOWN),
        "order": (lin.get("order") or UNKNOWN),
        "class": (lin.get("class") or UNKNOWN),
        "division": division if division else UNKNOWN,
    }


def enrich_with_taxonomy(
    df: pd.DataFrame,
    taxid_col: str = "taxid",
    path: Optional[Path | str] = None,
    overwrite: bool = True,
) -> pd.DataFrame:
    """Fill `species, genus, family, order, class, division` columns from NCBI.

    Operates in-place and returns the same frame for chaining. Existing values
    in those columns are overwritten when `overwrite=True` (the default — most
    adapters call this exactly to replace the placeholder values they wrote).
    """
    if taxid_col not in df.columns:
        raise KeyError(f"Column '{taxid_col}' not found in DataFrame")

    tax = get_taxonomy(path)
    # Batch unique taxids once — adapters frequently have repeats.
    uniq_taxids = (
        df[taxid_col].dropna().astype("Int64").drop_duplicates().tolist()
    )

    lineage_map: dict[int, dict[str, str]] = {}
    for t in uniq_taxids:
        if t is None:
            continue
        t = int(t)
        lin = tax.get_taxonomy(t) or {}
        division = tax.get_division(t) or UNKNOWN
        lineage_map[t] = {
            "species": (lin.get("species") or UNKNOWN),
            "genus": (lin.get("genus") or UNKNOWN),
            "family": (lin.get("family") or UNKNOWN),
            "order": (lin.get("order") or UNKNOWN),
            "class": (lin.get("class") or UNKNOWN),
            "division": division if division else UNKNOWN,
        }

    for rank in TAXRANKS:
        if rank in df.columns and not overwrite:
            continue

        def _lookup(t, _rank=rank):
            if pd.isna(t):
                return UNKNOWN
            return lineage_map.get(int(t), {}).get(_rank, UNKNOWN)

        df[rank] = df[taxid_col].map(_lookup).astype("string")

    return df
