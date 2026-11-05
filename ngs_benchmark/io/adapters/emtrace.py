"""Adapter for EMtrace / Lazypipe `abund_table.tsv` output.

Canonical format. See the project notes for the schema.
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

from ngs_benchmark.schema import REQUIRED_COLS, TAXRANKS


def load(path: str | Path) -> pd.DataFrame:
    """Load a EMtrace/Lazypipe abundance table.

    The file is a TSV with header. Required columns: readn, taxid, and the
    six taxonomic ranks (species, genus, family, order, class, division).
    Optional columns (passed through): contign, csum, prob, pvalue, bphage.

    `division` may appear as `superkingdom` or `phylum` — both are renamed
    to `division` to match the canonical schema.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=0,
        dtype={"taxid": "Int64"},
        keep_default_na=True,
        na_values=["", "NA", "NaN", "None"],
    )
    # Normalise top-rank column name.
    if "division" not in df.columns:
        for alt in ("superkingdom", "phylum"):
            if alt in df.columns:
                df = df.rename(columns={alt: "division"})
                break

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path}: missing required column(s) {missing}. Expected canonical "
            f"schema with {REQUIRED_COLS}. Found columns: {list(df.columns)}"
        )

    # Coerce types defensively.
    df["readn"] = pd.to_numeric(df["readn"], errors="coerce").fillna(0.0)
    for rank in TAXRANKS:
        df[rank] = df[rank].astype("string")
    return df
