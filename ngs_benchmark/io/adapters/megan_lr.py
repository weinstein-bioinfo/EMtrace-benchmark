"""Adapter for MEGAN6 / MEGAN-LR taxonomy-to-count exports.

MEGAN6 exports a "tax-to-count" CSV/TSV via `File → Export → Tax2Count`. The
file is a simple two-column table — the first column is either a taxonomic
name or an integer taxid (controlled by the "Use Names instead of IDs" option
in MEGAN). The second column is the read count assigned to that taxon by the
interval-union LCA algorithm.

We handle both flavours:

- **Taxid-first**: the integer column is treated as `taxid`; ranks are filled
  from NCBI.
- **Name-first**: we treat the string as `species` and try to resolve a taxid
  from the name (best-effort; if none, ranks remain "unknown" and the row
  scores at species rank only).
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd

from ngs_benchmark.io.taxonomy import (
    _resolve_db_path,
    enrich_with_taxonomy,
)
from ngs_benchmark.schema import REQUIRED_COLS, TAXRANKS, UNKNOWN


def _lookup_taxid_by_name(names: list[str], db_path: Path) -> dict[str, int | None]:
    """Reverse-lookup taxids by scientific name (uses the `organism` column)."""
    if not db_path.exists():
        return {n: None for n in names}
    out: dict[str, int | None] = {}
    with sqlite3.connect(db_path) as conn:
        cur = conn.cursor()
        placeholders = ",".join("?" for _ in names)
        try:
            cur.execute(
                f"SELECT organism, taxid FROM taxonomy WHERE organism IN ({placeholders})",
                names,
            )
            for org, tid in cur.fetchall():
                out[org] = int(tid)
        except sqlite3.OperationalError:
            # Legacy DB without `organism` column.
            return {n: None for n in names}
    for n in names:
        out.setdefault(n, None)
    return out


def load(path: str | Path) -> pd.DataFrame:
    # Try TSV first, then CSV.
    try:
        df = pd.read_csv(path, sep="\t", header=0, keep_default_na=True)
        if df.shape[1] < 2:
            raise ValueError("not enough columns")
    except (ValueError, pd.errors.ParserError):
        df = pd.read_csv(path, sep=",", header=0, keep_default_na=True)

    if df.shape[1] < 2:
        raise ValueError(f"{path}: MEGAN export needs at least 2 columns")

    col0, col1 = df.columns[0], df.columns[1]
    is_taxid = pd.to_numeric(df[col0], errors="coerce").notna().all()

    out = pd.DataFrame()
    out["readn"] = pd.to_numeric(df[col1], errors="coerce").fillna(0.0)
    if is_taxid:
        out["taxid"] = pd.to_numeric(df[col0], errors="coerce").astype("Int64")
    else:
        # Try to back-resolve taxid from species name.
        names = df[col0].astype("string").fillna("").tolist()
        name_map = _lookup_taxid_by_name(names, _resolve_db_path(None))
        out["taxid"] = pd.Series([name_map.get(n) for n in names], dtype="Int64")
        out["species"] = df[col0].astype("string")

    has_taxids = out["taxid"].notna().any()
    if has_taxids:
        enrich_with_taxonomy(out, taxid_col="taxid", overwrite=is_taxid)
    else:
        for rank in TAXRANKS:
            if rank not in out.columns:
                out[rank] = UNKNOWN

    missing = [c for c in REQUIRED_COLS if c not in out.columns]
    if missing:
        raise ValueError(f"{path}: megan_lr adapter missing column(s) {missing}")
    for rank in TAXRANKS:
        out[rank] = out[rank].astype("string")
    return out.reset_index(drop=True)
