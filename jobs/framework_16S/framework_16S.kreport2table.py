#!/usr/bin/env python3
"""Kraken 2 / Bracken kreport -> the scorer's `table` format.

Emits one row per taxon with `division, species, genus, family, readn`, which is what
jobs/framework_16S/framework_16S.score.py's `table` loader expects.

Why parse the kreport rather than Bracken's own TSV: the kreport carries the
whole lineage as indentation, so species, genus and family come from one file
with no taxonomy lookup. Bracken's `-w` writes a kreport with its corrected
abundances, so Kraken 2 and Bracken go through this same path and differ only
in which file is passed.

Read counts use the CLADE count at species level and the taxon-only count at
genus/family level. A read assigned no deeper than `Escherichia` is real mass
that belongs in the genus total but names no species, so it is emitted with an
empty species rather than dropped -- dropping it would silently inflate the
species proportions of everything else.

Usage:
    python3 framework_16S.kreport2table.py <kreport> <out.tsv>
"""
from __future__ import annotations

import sys
from pathlib import Path

# D (domain) becomes `division`: ngs_benchmark.detection.set_metrics filters
# predictions on that column, and a table without it contributes no TP and no
# FP at all -- it silently scores as if nothing had been predicted.
RANK_COL = {"D": "division", "S": "species", "G": "genus", "F": "family"}


def parse(kreport: Path):
    rows, stack = [], []
    for line in open(kreport, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        clade, taxon, rank, name_raw = int(f[1]), int(f[2]), f[3], f[5]
        depth = (len(name_raw) - len(name_raw.lstrip())) // 2
        name = name_raw.strip()

        # The indentation is the tree: drop anything at or below this depth.
        stack[depth:] = [(rank, name)]

        lin = {}
        for r, n in stack:
            # S1/S2 are strain sub-levels; fold them onto the species above.
            base = r[0] if r and r[0] in RANK_COL else None
            if base and base not in lin:
                lin[RANK_COL[base]] = n

        if rank.startswith("S"):
            if rank != "S":
                continue                      # strain rows are inside the species clade
            readn = clade                     # species clade: includes its strains
        elif rank in ("G", "F"):
            readn = taxon                     # only reads stopping at this rank
        else:
            continue
        if readn <= 0:
            continue
        rows.append({
            "division": lin.get("division", ""),
            "species": lin.get("species", ""),
            "genus": lin.get("genus", ""),
            "family": lin.get("family", ""),
            "readn": readn,
        })
    return rows


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    src, dst = Path(sys.argv[1]), Path(sys.argv[2])
    rows = parse(src)
    if not rows:
        print(f"{src}: no S/G/F rows found — is this a kreport?", file=sys.stderr)
        return 1
    dst.parent.mkdir(parents=True, exist_ok=True)
    with open(dst, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("division\tspecies\tgenus\tfamily\treadn\n")
        for r in rows:
            fh.write(f"{r['division']}\t{r['species']}\t{r['genus']}\t"
                     f"{r['family']}\t{r['readn']}\n")
    named = sum(1 for r in rows if r["species"])
    print(f"{dst}: {len(rows)} rows ({named} with a species), "
          f"{sum(r['readn'] for r in rows):,} reads")
    return 0


if __name__ == "__main__":
    sys.exit(main())
