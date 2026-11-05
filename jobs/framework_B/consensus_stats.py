#!/usr/bin/env python3
"""Summarise a consensus FASTA — or compare two of them.

Used by framework_B.emtrace.assemblies.bash to verify that the 2026-09-23
assembly changes actually reached consensus.fa.  The signal to look for is N
content: min_coverage 5 -> 20 masks every position under 20x depth.

    consensus_stats.py NEW            # summarise
    consensus_stats.py OLD NEW        # compare, per record
"""
from __future__ import annotations

import sys
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks)
            name = line[1:].split()[0] if len(line) > 1 else "unnamed"
            chunks = []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks)
    return records


def summarise(seq: str) -> tuple[int, int, float]:
    n = seq.upper().count("N")
    length = len(seq)
    return length, n, (100.0 * n / length if length else 0.0)


def main(argv: list[str]) -> int:
    if not 1 <= len(argv) <= 2:
        print(__doc__, file=sys.stderr)
        return 2

    if len(argv) == 1:
        recs = read_fasta(Path(argv[0]))
        print(f"  {'record':<20} {'len':>8} {'N':>8} {'%N':>7}  {'called':>7}")
        for name, seq in recs.items():
            length, n, pct = summarise(seq)
            print(f"  {name:<20} {length:>8} {n:>8} {pct:>6.2f}% {100 - pct:>6.2f}%")
        return 0

    old = read_fasta(Path(argv[0]))
    new = read_fasta(Path(argv[1]))
    print(f"  {'record':<20} {'len':>8} {'N before':>9} {'N after':>9} "
          f"{'%N before':>10} {'%N after':>9} {'delta %N':>9}")
    for name in sorted(set(old) | set(new)):
        o, nw = old.get(name, ""), new.get(name, "")
        ol, on, opct = summarise(o)
        nl, nn, npct = summarise(nw)
        length = nl or ol
        print(f"  {name:<20} {length:>8} {on:>9} {nn:>9} "
              f"{opct:>9.2f}% {npct:>8.2f}% {npct - opct:>+8.2f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
