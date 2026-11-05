"""Compose Module B into an end-to-end run.

Single run:

    align consensus→reference (MAFFT) → per-segment coverage/identity table
        → summary row (+ optional EMtrace viral_reads) → write outputs

Outputs (per ROADMAP §"Module B"):
    coverage.tsv   per-segment: segment, ref_len, covered_bp, aligned_bp,
                   matching_bp, coverage_pct, identity_pct
    summary.tsv    one row: [sample,] n_segments, total_ref_len, mean/weighted
                   coverage_pct + identity_pct, [viral_reads]
    config.json    provenance (config + summary)

Robustness: a missing/empty consensus yields NA coverage/identity (the assembly
either was not produced, or — with an explicit reference — is all-N). The
reference must exist (segment lengths are read from it).
"""

from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field
from pathlib import Path

import pandas as pd

from ngs_benchmark.io.writers import write_json, write_tsv

from .metrics import coverage_table, empty_coverage_table, summarise


@dataclass
class AssemblyConfig:
    consensus: str
    reference: str
    outdir: str
    sample: str | None = None
    #: EMtrace results/reports/stats_div.csv. The `read_count` for
    #: division=="Viruses" is reported as `viral_reads`; absent file → "no hits".
    stats_div: str | None = None
    extra: dict = field(default_factory=dict)


def count_viral_reads(stats_div: str | None) -> int | str | None:
    """`viral_reads` from a EMtrace results/reports/stats_div.csv.

    - path not given      → None      (column omitted / NA)
    - path given, missing → "no hits" (no viral reads reported)
    - path given, present → read_count for division=="Viruses" (0 if row absent)
    """
    if stats_div is None:
        return None
    if not os.path.exists(stats_div):
        return "no hits"
    import csv
    with open(stats_div, newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("division", "").strip() == "Viruses":
                try:
                    return int(row["read_count"])
                except (KeyError, ValueError):
                    return "no hits"
    return 0


def _consensus_usable(path: str | None) -> bool:
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def evaluate(cfg: AssemblyConfig) -> tuple[pd.DataFrame, dict]:
    """Compute the per-segment coverage table and summary row — no files written.

    Shared by `run` (single-sample, writes per-sample files) and `run_batch`
    (writes only combined TSVs).
    """
    if not os.path.exists(cfg.reference):
        raise FileNotFoundError(f"reference not found: {cfg.reference}")

    if _consensus_usable(cfg.consensus):
        cov_df = coverage_table(cfg.consensus, cfg.reference)
    else:
        cov_df = empty_coverage_table(cfg.reference)

    summary = summarise(cov_df)
    viral_reads = count_viral_reads(cfg.stats_div)
    if viral_reads is not None:
        summary["viral_reads"] = viral_reads
        cov_df = cov_df.copy()
        cov_df["viral_reads"] = viral_reads
    if cfg.sample is not None:
        summary = {"sample": cfg.sample, **summary}
    return cov_df, summary


def run(cfg: AssemblyConfig) -> dict:
    """Run one assembly-quality evaluation. Returns the summary row dict.

    Writes per-sample ``coverage.tsv``, ``summary.tsv`` and ``config.json`` to
    ``cfg.outdir`` (the ROADMAP Module B single-sample layout).
    """
    cov_df, summary = evaluate(cfg)

    outdir = Path(cfg.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(cov_df, outdir / "coverage.tsv")
    write_tsv(pd.DataFrame([summary]), outdir / "summary.tsv")
    write_json(
        {"config": {k: v for k, v in asdict(cfg).items() if k != "extra"},
         "summary": summary},
        outdir / "config.json",
    )
    return summary


def run_batch(manifest: pd.DataFrame, outdir: str) -> pd.DataFrame:
    """Evaluate many consensus/reference pairs from a manifest.

    Manifest columns: ``sample_id``, ``consensus``, ``reference``, optional
    ``candidates``. Writes exactly two TSVs at ``<outdir>/`` — no per-sample
    subdirectories:

        summary.tsv    one row per sample (leading ``sample`` column)
        coverage.tsv   one row per reference segment per sample (leading
                       ``sample`` column)

    Missing consensus paths are tolerated (NA metrics). Returns the summary
    DataFrame.
    """
    outroot = Path(outdir)
    outroot.mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict] = []
    coverage_frames: list[pd.DataFrame] = []
    for _, r in manifest.iterrows():
        sample_id = str(r["sample_id"])
        stats_div = r["stats_div"] if "stats_div" in manifest.columns and pd.notna(r.get("stats_div")) else None
        cfg = AssemblyConfig(
            consensus=str(r["consensus"]),
            reference=str(r["reference"]),
            outdir=str(outroot),          # unused: run_batch never writes per-sample files
            sample=sample_id,
            stats_div=str(stats_div) if stats_div is not None else None,
        )
        cov_df, summary = evaluate(cfg)
        summary_rows.append(summary)
        cov_df = cov_df.copy()
        cov_df.insert(0, "sample", sample_id)
        coverage_frames.append(cov_df)

    combined_summary = pd.DataFrame(summary_rows)
    combined_coverage = (
        pd.concat(coverage_frames, ignore_index=True)
        if coverage_frames else pd.DataFrame()
    )
    write_tsv(combined_summary, outroot / "summary.tsv")
    write_tsv(combined_coverage, outroot / "coverage.tsv")
    return combined_summary
