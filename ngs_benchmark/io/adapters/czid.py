"""Adapter for CZ ID output: per-sample Taxon Report CSV, or `taxpath.tsv`.

CZ ID emits both nucleotide (nt) and protein (nr) read counts per taxon. They
are resolved as

    readn = max(nt_count, nr_count)        (NA counted as 0)

The legacy rule (`BMARK.simmet.part1.R:48-53`) took nt_count whenever it was
> 0 and nr_count only otherwise. That drops a well-supported protein-level call
whenever a single read also matched at nucleotide level: on the Buddle DNA mock
at 60,000 gc/mL, cytomegalovirus (a spiked, true virus) had nt 1 / nr 75 and was
scored as 1 read, below the 10-read threshold (Supplementary Table 2).

and drops any row with `is_phage == True`.

Taxon Report CSV (from the CZ ID web UI; handling ported from the
framework-A-code branch, commit 6c5d119): comma-separated, `tax_id` plus a
numeric `genus_tax_id`, but no resolved lineage. `tax_id` is renamed to `taxid`
and the missing rank columns are filled from the shared NCBI taxonomy DB
(`ngs_benchmark.io.taxonomy.enrich_with_taxonomy`), as for `megan_lr`/`metamix`.
CZ ID also reports a genus-level aggregate row (`tax_level=2`) beside each
species-level row (`tax_level=1`) for the same hit; left in, it becomes a
phantom `species="unknown"` false positive and double-counts at genus/family,
so only `tax_level == 1` is kept. A pre-resolved `taxpath.tsv` (full lineage,
no `tax_level`) passes through unchanged.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from ngs_benchmark.io.taxonomy import enrich_with_taxonomy
from ngs_benchmark.schema import REQUIRED_COLS, TAXRANKS

logger = logging.getLogger(__name__)


def load(path: str | Path) -> pd.DataFrame:
    # CZ ID Taxon Reports are comma-separated (`*_taxon_report.csv`); older
    # taxpath exports are tab-separated. Sniff the delimiter from the extension,
    # defaulting to comma.
    sep = "\t" if str(path).endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(
        path,
        sep=sep,
        header=0,
        keep_default_na=True,
        na_values=["", "NA", "NaN", "None"],
    )
    if "taxid" not in df.columns and "tax_id" in df.columns:
        df = df.rename(columns={"tax_id": "taxid"})
    if "taxid" in df.columns:
        df["taxid"] = pd.to_numeric(df["taxid"], errors="coerce").astype("Int64")
    if "tax_level" in df.columns:
        # keep species-level rows; genus aggregates would be phantom "unknown" FPs
        tax_level = pd.to_numeric(df["tax_level"], errors="coerce")
        df = df[tax_level.isna() | (tax_level == 1)].copy()

    if "is_phage" in df.columns:
        is_phage_mask = df["is_phage"].astype("string").str.lower().isin(["true", "yes", "1"])
        df = df[~is_phage_mask].copy()
        df["bphage"] = "no"

    nt = pd.to_numeric(df.get("nt_count"), errors="coerce") if "nt_count" in df.columns else None
    nr = pd.to_numeric(df.get("nr_count"), errors="coerce") if "nr_count" in df.columns else None
    if nt is not None and nr is not None:
        readn = pd.concat([nt.fillna(0.0), nr.fillna(0.0)], axis=1).max(axis=1)
    elif nt is not None:
        readn = nt.fillna(0.0)
    elif nr is not None:
        readn = nr.fillna(0.0)
    elif "readn" in df.columns:
        readn = pd.to_numeric(df["readn"], errors="coerce").fillna(0.0)
    else:
        raise ValueError(
            f"{path}: no nt_count/nr_count/readn column found. "
            f"Found columns: {list(df.columns)}"
        )
    df["readn"] = readn

    if "division" not in df.columns:
        for alt in ("superkingdom", "phylum"):
            if alt in df.columns:
                df = df.rename(columns={alt: "division"})
                break

    missing_ranks = [r for r in TAXRANKS if r not in df.columns]
    if missing_ranks:
        if "taxid" not in df.columns:
            raise ValueError(
                f"{path}: CZ ID input has no lineage column(s) {missing_ranks} "
                f"and no 'taxid'/'tax_id' to resolve them from. Found: {list(df.columns)}"
            )
        enrich_with_taxonomy(df, taxid_col="taxid", overwrite=False)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path}: CZ ID adapter missing required column(s) after normalisation: "
            f"{missing}. Found: {list(df.columns)}"
        )

    for rank in TAXRANKS:
        df[rank] = df[rank].astype("string")
    return df.reset_index(drop=True)
