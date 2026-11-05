"""Adapter for Kraken2 + optional Bracken output.

Two input modes are supported, controlled by which file the user points at:

A. Kraken2 kreport (`*kreport.txt` / `*report.txt`) — the standard 6-column
   indented hierarchy:

       %    reads_clade    reads_taxon    rank    taxid    name

   We extract rows where `rank == "S"` (species level). `reads_taxon` is the
   read count uniquely assigned to that species. The kreport's `name` field
   is leading-space-indented; we strip it for the canonical schema.

B. Kraken + Bracken merged — when the file has 7 columns and a header line
   matching Bracken (`new_est_reads` column present), we delegate to the
   EPI2ME / Bracken parser (`new_est_reads` wins over `reads_taxon`).

The canonical rank columns are filled via the NCBI taxonomy DB.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ngs_benchmark.io.taxonomy import enrich_with_taxonomy
from ngs_benchmark.schema import REQUIRED_COLS, TAXRANKS


def _detect_format(path: Path) -> str:
    """Peek at first header line; return 'bracken' or 'kreport'."""
    with open(path, encoding="utf-8") as fh:
        first = fh.readline()
    if "new_est_reads" in first or "kraken_assigned_reads" in first:
        return "bracken"
    return "kreport"


def _load_kreport(path: Path) -> pd.DataFrame:
    """Parse a Kraken2 kreport into a DataFrame at species rank."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["pct", "reads_clade", "reads_taxon", "rank", "taxid", "name"],
        dtype={"taxid": "Int64"},
        keep_default_na=False,
    )
    # Species-rank rows only. Kraken's "S" includes "S1", "S2" sub-species ranks;
    # we keep only the canonical "S".
    species = df[df["rank"].astype("string").str.upper() == "S"].copy()
    species["readn"] = pd.to_numeric(species["reads_taxon"], errors="coerce").fillna(0.0)
    species["name"] = species["name"].astype("string").str.lstrip()
    return species[["readn", "taxid", "name"]].reset_index(drop=True)


def load(path: str | Path) -> pd.DataFrame:
    p = Path(path)
    fmt = _detect_format(p)
    if fmt == "bracken":
        # Delegate to the EPI2ME adapter — it handles the Bracken schema verbatim.
        from . import epi2me
        return epi2me.load(p)

    df = _load_kreport(p)
    enrich_with_taxonomy(df, taxid_col="taxid")

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: kraken_bracken adapter missing column(s) {missing}")
    for rank in TAXRANKS:
        df[rank] = df[rank].astype("string")
    return df.reset_index(drop=True)
