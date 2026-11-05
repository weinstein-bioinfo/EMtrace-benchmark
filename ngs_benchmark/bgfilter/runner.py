"""Compose Module D into runnable single / batch evaluations.

Single:  one removed-reads FASTQ (one filter setting) → one table row.
Batch:   a manifest of (filter, removed_fastq) pairs → one table, one row per
         setting (e.g. AS0..AS125 × {default, k11w5}), reproducing LazypipeX
         Table 1. The truth sets are loaded once and shared across all rows.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from ngs_benchmark.io.writers import write_tsv

from .metrics import (
    COLUMN_ORDER,
    classify_removed,
    classify_removed_division,
    division_columns,
)
from .truth import BgfilterTruth, DivisionTruth


@dataclass
class BgfilterConfig:
    """One filter setting to score."""

    #: Reads the host filter removed (emtrace ``reads/hflt.<host>.reads.fq.gz``).
    removed_fastq: str
    #: Filter label for the table row (e.g. ``AS100`` or ``k11w5.AS100``).
    label: str


def evaluate(cfg: BgfilterConfig, truth: BgfilterTruth, strip_mate: bool = False) -> dict:
    """Per-category fraction-removed row for one filter setting — no files written."""
    return classify_removed(cfg.removed_fastq, truth, cfg.label, strip_mate=strip_mate)


def run(cfg: BgfilterConfig, truth: BgfilterTruth, outdir: str, strip_mate: bool = False) -> dict:
    """Score one setting and write a single-row ``table.tsv`` to *outdir*."""
    row = evaluate(cfg, truth, strip_mate=strip_mate)
    df = pd.DataFrame([row], columns=COLUMN_ORDER)
    write_tsv(df, Path(outdir) / "table.tsv")
    return row


def run_batch(
    manifest: pd.DataFrame, truth: BgfilterTruth, outdir: str, strip_mate: bool = False
) -> pd.DataFrame:
    """Score many settings from a manifest, writing one combined ``table.tsv``.

    Manifest columns: ``filter``, ``removed_fastq``. Rows are emitted in manifest
    order so the sweep reads top-to-bottom (AS0 → AS125, per config). *strip_mate*
    drops a trailing ``/1`` / ``/2`` so paired reads match per-fragment truth ids.
    """
    rows: list[dict] = []
    for _, r in manifest.iterrows():
        cfg = BgfilterConfig(
            removed_fastq=str(r["removed_fastq"]),
            label=str(r["filter"]),
        )
        rows.append(evaluate(cfg, truth, strip_mate=strip_mate))

    df = pd.DataFrame(rows, columns=COLUMN_ORDER)
    write_tsv(df, Path(outdir) / "table.tsv")
    return df


def run_division_batch(
    manifest: pd.DataFrame, truth: DivisionTruth, outdir: str, strip_mate: bool = False
) -> pd.DataFrame:
    """Score many settings against a per-read division truth table.

    For the simulated-Nanopore datasets: one row per filter setting (e.g. an
    ``-s`` sweep × seeding config), columns = fraction of reads removed per
    division (Viruses / Bacteria / Eukaryota). Manifest columns: ``filter``,
    ``removed_fastq``. Writes ``<outdir>/table.tsv``.
    """
    rows: list[dict] = []
    for _, r in manifest.iterrows():
        rows.append(classify_removed_division(
            str(r["removed_fastq"]), truth, str(r["filter"]), strip_mate=strip_mate))

    df = pd.DataFrame(rows, columns=division_columns(truth))
    write_tsv(df, Path(outdir) / "table.tsv")
    return df
