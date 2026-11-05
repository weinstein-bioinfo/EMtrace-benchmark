"""Per-category fraction-removed metrics for background filtering (Module D).

Reproduces the LazypipeX "Table 1" benchmark: given the set of reads a host
filter *removed* (one FASTQ per filter setting), report the fraction of reads
removed in each true-origin category. A good host filter removes ~100% of
human/bacterial reads while retaining viral reads.

The removed FASTQ is streamed once; only the small truth id-sets are held in
memory (see :mod:`ngs_benchmark.bgfilter.truth`).
"""

from __future__ import annotations

import gzip
from pathlib import Path

from .truth import BgfilterTruth, DivisionTruth

#: Internal category key → human-readable column label (and emission order).
CATEGORY_LABELS: dict[str, str] = {
    "total": "Total",
    "homo_sapiens": "Homo sapiens",
    "bacteria": "Bacteria",
    "viruses": "Viruses",
    "herv": "HERV",
    "viruses_excl_retro": "Viruses excl. Retroviridae",
}

#: Column order for the output table (leading ``filter`` label then categories).
COLUMN_ORDER: list[str] = ["filter", *CATEGORY_LABELS.keys()]


def _strip_mate_suffix(token: str) -> str:
    """Drop a trailing ``/1`` / ``/2`` (mate) suffix so paired reads match a
    per-fragment truth id (e.g. ``gi|..-35000/1`` → ``gi|..-35000``).

    Only a trailing ``/<digits>`` is removed; ids without one are returned as-is.
    """
    base, sep, tail = token.rpartition("/")
    return base if (sep and tail.isdigit()) else token


def iter_fastq_ids(path: str | Path, strip_mate: bool = False):
    """Yield the read id (first whitespace token of each header) of a FASTQ.

    Handles plain or gzipped FASTQ. Assumes the canonical 4-lines-per-record
    layout (sequence not wrapped), so every 4th line is a header. When
    *strip_mate* is True a trailing ``/1`` / ``/2`` is removed (paired reads
    whose truth is keyed per fragment).
    """
    p = str(path)
    opener = gzip.open(p, "rt") if p.endswith(".gz") else open(p)
    with opener as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            if not line:
                continue
            head = line[1:] if line[0] == "@" else line
            token = head.split(None, 1)[0] if head.strip() else ""
            if token:
                yield _strip_mate_suffix(token) if strip_mate else token


def count_removed(
    removed_fastq: str | Path, truth: BgfilterTruth, strip_mate: bool = False
) -> dict[str, int]:
    """Count removed reads per category by streaming the removed FASTQ once.

    Uses the set containments ``vi ⊆ nohuman`` and ``vi_nophix_noretro, phix ⊆ vi``
    so each read needs at most three set lookups. *strip_mate* drops a trailing
    ``/1`` / ``/2`` so paired reads match per-fragment truth ids.
    """
    nohuman = truth.nohuman
    vi = truth.vi
    vinpnr = truth.vi_nophix_noretro
    phix = truth.phix

    r_total = r_nohuman = r_vi = r_vinpnr = r_phix = 0
    for rid in iter_fastq_ids(removed_fastq, strip_mate=strip_mate):
        r_total += 1
        if rid in nohuman:
            r_nohuman += 1
            if rid in vi:
                r_vi += 1
                if rid in vinpnr:
                    r_vinpnr += 1
                elif rid in phix:
                    r_phix += 1

    return {
        "total": r_total,
        "homo_sapiens": r_total - r_nohuman,
        "bacteria": r_nohuman - r_vi,
        "viruses": r_vi,
        "viruses_excl_retro": r_vinpnr + r_phix,
        "herv": r_vi - r_vinpnr - r_phix,
    }


def fractions_removed(
    removed: dict[str, int], denominators: dict[str, int]
) -> dict[str, float | None]:
    """Per-category percentage removed = removed/denominator × 100 (None if denom 0)."""
    out: dict[str, float | None] = {}
    for cat in CATEGORY_LABELS:
        d = denominators.get(cat, 0)
        out[cat] = round(removed[cat] / d * 100, 3) if d else None
    return out


def classify_removed(
    removed_fastq: str | Path, truth: BgfilterTruth, label: str, strip_mate: bool = False
) -> dict[str, float | None | str]:
    """One table row: ``filter`` label + per-category fraction removed (%)."""
    removed = count_removed(removed_fastq, truth, strip_mate=strip_mate)
    fracs = fractions_removed(removed, truth.denominators)
    row: dict[str, float | None | str] = {"filter": label}
    row.update(fracs)
    return row


# --- division-labeled truth (simulated-Nanopore datasets) --------------------

def count_removed_division(
    removed_fastq: str | Path, truth: DivisionTruth, strip_mate: bool = False
) -> dict[str, int]:
    """Count removed reads per division by streaming the removed FASTQ once.

    Returns ``{"total": N, <division>: count, ...}``. Reads not found in the truth
    table are counted in ``total`` only (should not happen — removed reads are a
    subset of the input).
    """
    div_of = truth.div_of
    removed: dict[str, int] = {"total": 0}
    for rid in iter_fastq_ids(removed_fastq, strip_mate=strip_mate):
        removed["total"] += 1
        div = div_of.get(rid)
        if div is not None:
            removed[div] = removed.get(div, 0) + 1
    return removed


def division_columns(truth: DivisionTruth) -> list[str]:
    """Output column order for a division table: ``filter, total, <divisions…>``."""
    return ["filter", "total", *truth.divisions]


def classify_removed_division(
    removed_fastq: str | Path, truth: DivisionTruth, label: str, strip_mate: bool = False
) -> dict[str, float | None | str]:
    """One row: ``filter`` + total + per-division fraction removed (%).

    A good host filter maximises Eukaryota/Bacteria removal while minimising
    Viruses removal — i.e. lets you pick the optimal ``-s``.
    """
    removed = count_removed_division(removed_fastq, truth, strip_mate=strip_mate)
    denom = truth.denominators
    row: dict[str, float | None | str] = {"filter": label}
    for cat in ("total", *truth.divisions):
        d = denom.get(cat, 0)
        row[cat] = round(removed.get(cat, 0) / d * 100, 3) if d else None
    return row
