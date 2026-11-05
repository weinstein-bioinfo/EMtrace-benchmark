#!/usr/bin/env python3
"""Framework 16S — convert Spaghetti output to a table the scorer can read.

Spaghetti writes a phyloseq pair keyed by SILVA 138 reference ID:

    otu_table.csv          #OTU ID, <sample>          primary-alignment read counts
    phyloseq_taxonomy.csv  #OTU ID, Domain ... Species QIIME-prefixed labels

This joins them and sums reads per lineage, writing one row per
genus x species:

    division  genus  species  readn

`readn` carries the real read counts so the scorer's readn_min sweep means the
same thing as for Emu and EMtrace.

The converter only reformats. It does NOT fix labels: SILVA 138 predates the
2020 Lactobacillus split and files L. fermentum under `g__Lactobacillus`, and
that rename is a scoped rule in data/truth/zymobiomics_16S_collapse.tsv,
applied by the scorer where it is visible. The `species` column is kept so
that rule can match; Spaghetti is still scored at genus only, because the
manifest marks it `taxonomy=silva` and SILVA 138 carries no species keys.

The one thing dropped here is unassigned mass (genus empty, `uncultured`,
`unidentified`), which must not be counted as a named genus
(the project notes).

Usage:
    python3 jobs/framework_16S/framework_16S.spaghetti2table.py \
        results/framework_16S/spaghetti data/comparators/spaghetti.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

UNASSIGNED = re.compile(r"^(|uncultured|unidentified|unknown)$", re.IGNORECASE)


def strip_prefix(v: object) -> str:
    """`g__Bacillus` -> `Bacillus`; `s__Lactobacillus_fermentum` -> `Lactobacillus fermentum`."""
    if not isinstance(v, str):
        return ""
    return re.sub(r"^[a-z]__", "", v.strip()).replace("_", " ")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir", help="Spaghetti output directory (otu_table.csv, phyloseq_taxonomy.csv)")
    ap.add_argument("out", help="table to write")
    args = ap.parse_args()

    rundir = Path(args.rundir)
    otu = pd.read_csv(rundir / "otu_table.csv")
    if otu.shape[1] != 2:
        print(f"ERROR: expected one sample column, found {otu.shape[1] - 1}", file=sys.stderr)
        return 1
    otu.columns = ["id", "reads"]
    tax = pd.read_csv(rundir / "phyloseq_taxonomy.csv").rename(columns={"#OTU ID": "id"})
    df = otu.merge(tax, on="id", how="left", validate="one_to_one")
    if df["Genus"].isna().all():
        print("ERROR: no taxonomy joined to the OTU table", file=sys.stderr)
        return 1

    df["division"] = df["Domain"].map(strip_prefix)
    df["genus"] = df["Genus"].map(strip_prefix)
    df["species"] = df["Species"].map(strip_prefix)

    unassigned = df["genus"].str.match(UNASSIGNED)
    out = (
        df[~unassigned]
        .groupby(["division", "genus", "species"], as_index=False)["reads"].sum()
        .rename(columns={"reads": "readn"})
        .sort_values("readn", ascending=False)
    )
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    print(f"reads assigned        {int(df['reads'].sum()):,}")
    print(f"unassigned genus      {int(df.loc[unassigned, 'reads'].sum()):,} (excluded)")
    print(f"reads written         {int(out['readn'].sum()):,} in {len(out):,} genus x species rows, "
          f"{out['genus'].nunique():,} genera")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
