"""Adapter for ONT EPI2ME wf-metagenomics output.

The per-sample Bracken output `*.bracken.tsv` (per `epi2me-labs/wf-metagenomics`)
has the standard Bracken columns:

    name                      taxonomy name
    taxonomy_id               NCBI taxid
    taxonomy_lvl              one-letter rank code (S, G, F, …)
    kraken_assigned_reads     reads directly assigned by Kraken2
    added_reads               reads re-distributed by Bracken
    new_est_reads             Bracken's re-estimated read count
    fraction_total_reads      relative abundance

We map to canonical schema with:

    readn  = new_est_reads (fallback: kraken_assigned_reads + added_reads)
    taxid  = taxonomy_id

Rank columns (species/genus/.../division) are filled via NCBITaxonomy.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ngs_benchmark.io.taxonomy import enrich_with_taxonomy
from ngs_benchmark.schema import REQUIRED_COLS, TAXRANKS


def load(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        header=0,
        dtype={"taxonomy_id": "Int64"},
        keep_default_na=True,
        na_values=["", "NA", "NaN", "None"],
    )

    if "taxonomy_id" in df.columns:
        df = df.rename(columns={"taxonomy_id": "taxid"})
    elif "taxid" not in df.columns:
        raise ValueError(
            f"{path}: no `taxonomy_id` or `taxid` column. Found: {list(df.columns)}"
        )

    if "new_est_reads" in df.columns:
        readn = pd.to_numeric(df["new_est_reads"], errors="coerce")
    elif "kraken_assigned_reads" in df.columns:
        added = pd.to_numeric(df.get("added_reads", 0), errors="coerce").fillna(0)
        kraken = pd.to_numeric(df["kraken_assigned_reads"], errors="coerce")
        readn = (kraken.fillna(0) + added)
    elif "count" in df.columns:
        readn = pd.to_numeric(df["count"], errors="coerce")
    else:
        raise ValueError(
            f"{path}: no `new_est_reads` / `kraken_assigned_reads` / `count` column"
        )
    df["readn"] = readn.fillna(0.0)

    # Restrict to species-level rows when the Bracken output mixes ranks.
    if "taxonomy_lvl" in df.columns:
        df = df[df["taxonomy_lvl"].astype("string").str.upper() == "S"].copy()

    enrich_with_taxonomy(df, taxid_col="taxid")

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path}: epi2me adapter missing required column(s) {missing}"
        )
    for rank in TAXRANKS:
        df[rank] = df[rank].astype("string")
    return df.reset_index(drop=True)
