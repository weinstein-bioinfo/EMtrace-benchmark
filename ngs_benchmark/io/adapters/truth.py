"""Truth-table loader.

Two flavours of truth file are supported (per the project notes):

1. Membership truth — `taxid, species, genus, ..., division [, readn] [, gc_per_ml]`.
   Taxa listed are present; everything else is absent.
2. Dilution / Buddle truth — adds `sample_id`, `gc_per_ml`, `nc_sample_id`, and
   optional `matrix_id` columns. One row per (virus × sample).

Auto-detection rule: presence of `nc_sample_id` ⇒ dilution flavour.
"""

from __future__ import annotations

from pathlib import Path
import pandas as pd

from ngs_benchmark.schema import TAXRANKS


REQUIRED_TRUTH_COLS = ("taxid", *TAXRANKS)


def load(path: str | Path) -> pd.DataFrame:
    """Load a truth table.

    Returns a DataFrame with at least `REQUIRED_TRUTH_COLS`. `readn` is added
    if missing (filled with 1.0 to indicate presence). The flavour is recorded
    in the returned DataFrame's `.attrs["truth_flavour"]`.
    """
    df = pd.read_csv(
        path,
        sep="\t",
        header=0,
        dtype={"taxid": "Int64"},
        keep_default_na=True,
        na_values=["", "NA", "NaN", "None"],
    )
    if "division" not in df.columns:
        for alt in ("superkingdom", "phylum"):
            if alt in df.columns:
                df = df.rename(columns={alt: "division"})
                break

    missing = [c for c in REQUIRED_TRUTH_COLS if c not in df.columns]
    if missing:
        raise ValueError(
            f"{path}: truth file missing required column(s) {missing}. "
            f"Found columns: {list(df.columns)}"
        )

    if "readn" not in df.columns:
        df["readn"] = 1.0
    df["readn"] = pd.to_numeric(df["readn"], errors="coerce").fillna(1.0)
    for rank in TAXRANKS:
        df[rank] = df[rank].astype("string")

    flavour = "dilution" if "nc_sample_id" in df.columns else "membership"
    df.attrs["truth_flavour"] = flavour
    return df
