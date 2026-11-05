"""Adapter for metaMix species-assignments output.

metaMix (Morfopoulou & Plagnol 2015) emits a TSV summary with one row per
candidate species, containing at minimum:

    Species         scientific name
    Reads           assigned read count (sometimes `assignedReads`)
    MarginalLogLike Bayesian marginal log-likelihood (optional, kept as `prob_raw`)

The exact column names vary by metaMix version (`Species_name`, `species`,
`Taxon`, ...). We look for the most common spellings.

If the file already carries a `taxid` / `TaxId` column, we use it; otherwise
we back-resolve from the species name.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ngs_benchmark.io.adapters.megan_lr import _lookup_taxid_by_name
from ngs_benchmark.io.taxonomy import _resolve_db_path, enrich_with_taxonomy
from ngs_benchmark.schema import REQUIRED_COLS, TAXRANKS, UNKNOWN


_NAME_COLS = ("Species", "species", "Species_name", "Taxon", "taxon", "name", "organism")
_READ_COLS = ("Reads", "reads", "assignedReads", "AssignedReads", "n_reads", "count")
_TAXID_COLS = ("taxid", "TaxId", "taxId", "tax_id", "NCBI_taxid")
_PROB_COLS = ("MarginalLogLike", "marginalLogLike", "posterior", "prob")


def _first_present(df: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def load(path: str | Path) -> pd.DataFrame:
    # metaMix is sometimes CSV, sometimes TSV — sniff.
    with open(path, encoding="utf-8") as fh:
        first = fh.readline()
    sep = "\t" if "\t" in first else ","
    df = pd.read_csv(path, sep=sep, header=0, keep_default_na=True)

    read_col = _first_present(df, _READ_COLS)
    name_col = _first_present(df, _NAME_COLS)
    taxid_col = _first_present(df, _TAXID_COLS)
    prob_col = _first_present(df, _PROB_COLS)
    if read_col is None:
        raise ValueError(
            f"{path}: metaMix adapter could not find a read-count column. "
            f"Tried {_READ_COLS}. Found: {list(df.columns)}"
        )
    if name_col is None and taxid_col is None:
        raise ValueError(
            f"{path}: metaMix adapter needs either a species name column "
            f"({_NAME_COLS}) or a taxid column ({_TAXID_COLS})"
        )

    out = pd.DataFrame()
    out["readn"] = pd.to_numeric(df[read_col], errors="coerce").fillna(0.0)
    if prob_col is not None:
        out["prob"] = pd.to_numeric(df[prob_col], errors="coerce")

    if taxid_col is not None:
        out["taxid"] = pd.to_numeric(df[taxid_col], errors="coerce").astype("Int64")
    else:
        names = df[name_col].astype("string").fillna("").tolist()
        name_map = _lookup_taxid_by_name(names, _resolve_db_path(None))
        out["taxid"] = pd.Series([name_map.get(n) for n in names], dtype="Int64")
        out["species"] = df[name_col].astype("string")

    if out["taxid"].notna().any():
        enrich_with_taxonomy(out, taxid_col="taxid", overwrite=(taxid_col is not None))
    else:
        for rank in TAXRANKS:
            if rank not in out.columns:
                out[rank] = UNKNOWN

    missing = [c for c in REQUIRED_COLS if c not in out.columns]
    if missing:
        raise ValueError(f"{path}: metamix adapter missing column(s) {missing}")
    for rank in TAXRANKS:
        out[rank] = out[rank].astype("string")
    return out.reset_index(drop=True)
