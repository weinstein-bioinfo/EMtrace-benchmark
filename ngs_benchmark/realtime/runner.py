"""Compose Module C into a runnable replay: drive a pipeline, score each batch.

The pipeline is kept behind a small protocol so this package does not import
EMtrace — the same separation the ``io.adapters`` use. A driver has to do two
things: process whatever is in the watch directory, and hand back the cumulative
abundance table once a batch is done.

    replay.py   releases batch files and measures both clocks
    runner.py   (here) drives the pipeline and scores each snapshot
    metrics.py  turns the snapshots into the three output tables

Every batch's cumulative abundance table is also kept, as
``abund/batch_NNNN.tsv`` beside the tables, so a finished replay can be scored
again under different settings (:func:`rescore`) without replaying it.
"""

from __future__ import annotations

import os
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import pandas as pd

from ngs_benchmark.io.writers import write_tsv
from ngs_benchmark.tailtrim import resolve_count_col, tail_trim

# Import the names, not the module: the package __init__ re-exports a *function*
# called `replay`, so `from . import replay` would bind that function.
from .replay import contiguity_gaps, replay as replay_batches, scan_library
from .metrics import (
    PANEL_OUTCOMES,
    Snapshot,
    accuracy,
    first_call,
    first_detection,
    latency,
    panel,
    panel_summary,
    set_scores,
    timeseries,
)


#: Truth-table abundance column; rows at 0 are absent taxa (e.g. yeasts in a 16S
#: truth). Matches abundance_metrics' default truth_weight.
TRUTH_WEIGHT = "pct_16S_only"


class PipelineDriver(Protocol):
    """What a real-time pipeline must expose to be replayed against."""

    def start(self, watch_dir: Path, outdir: Path) -> None:
        """Begin watching *watch_dir*; return once ready for the first batch."""

    def await_batch(self, index: int) -> None:
        """Block until batch *index* has been processed and reported."""

    def abundance(self) -> pd.DataFrame | None:
        """Current cumulative abundance table, or ``None`` if not written yet."""

    def stop(self) -> None:
        """Shut the pipeline down."""


@dataclass
class RealtimeConfig:
    """One library to replay."""

    #: MinKNOW ``fastq_pass`` directory; per-read ``start_time`` headers required.
    fastq_dir: str
    #: Name used in every output table.
    library: str
    #: Expected dominant taxon, for the first-correct-call metric.
    truth_taxon: str | None = None
    #: Multi-taxon truth table (profiling runs) — enables the L1 column.
    truth_table: str | None = None
    #: Rank the truth table and the abundance table are compared at.
    rank: str = "species"
    #: Column holding the EM posterior in the pipeline's abundance table.
    prob_col: str = "prob"
    #: Column holding per-taxon read counts.
    score_col: str = "readn"
    #: Taxa below this many reads are ignored when scoring — without it every
    #: one-read noise taxon counts as a false positive at every batch. 10 matches
    #: the Framework A / 16S ``readn_min`` and Emu's own minimum.
    read_min: float = 10.0
    #: Cumulative-share tail trim as a FRACTION, applied as a ROW filter to
    #: the abundance table before the per-rank sum below. 0.0 = off, which is
    #: what every completed Module C result was produced with.
    tail: float = 0.0
    #: Count column the tail trim ranks on; falls back to `score_col`.
    tail_col: str = "readn"
    #: Curve that defines stabilisation: ``L1`` for profiling, ``F1`` otherwise.
    stabilise_on: str = "L1"
    #: Tolerance band for stabilisation.
    stabilise_tol: float = 0.02
    #: Refuse libraries whose sequencing time jumps more than this (combined runs).
    max_gap_h: float = 24.0
    #: Hard-link batch files into the watch dir instead of copying.
    hardlink: bool = False
    #: Library has no ``start_time`` headers (SRA deposit): replay on read count
    #: only. ``seq_time_s`` and the throughput ratio are left empty.
    reads_only: bool = False
    #: Field-library panel (see :func:`load_panel`): taxa the source study
    #: reports, scored detected yes/no. Replaces ``truth_taxon`` / ``truth_table``
    #: — a field library has no complete truth, so no Pr/Rc/FP are computed.
    panel: str | None = None


PANEL_REQUIRED = ("taxon", "evidence")


def load_panel(path: str | Path, rank: str) -> pd.DataFrame:
    """Read a field-library panel TSV.

    Columns: ``taxon`` and ``evidence`` (``pcr_pos``, ``pcr_neg`` or
    ``reported``) are required; ``source`` (the label reported with the result,
    e.g. "Ergunay 2023, Table S4 (BLAST)") and ``rank`` (default: the run's
    rank) are optional. Lines starting with ``#`` are comments.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#").fillna("")
    missing = set(PANEL_REQUIRED) - set(df.columns)
    if missing:
        raise ValueError(f"{path}: panel is missing column(s) {', '.join(sorted(missing))}")
    df["taxon"] = df["taxon"].str.strip()
    df["evidence"] = df["evidence"].str.strip()
    bad = sorted(set(df["evidence"]) - set(PANEL_OUTCOMES))
    if bad:
        raise ValueError(f"{path}: unknown evidence {bad}; use {sorted(PANEL_OUTCOMES)}")
    dup = df.loc[df["taxon"].duplicated(), "taxon"].tolist()
    if dup:
        raise ValueError(f"{path}: taxon listed twice: {dup}")
    if "source" not in df.columns:
        df["source"] = ""
    if "rank" not in df.columns:
        df["rank"] = ""
    df["rank"] = df["rank"].where(df["rank"] != "", rank)
    return df[["taxon", "rank", "evidence", "source"]].reset_index(drop=True)


class EMtraceDriver:
    """Drive ``emtrace watch`` as a subprocess.

    The command mirrors the ``watch`` CLI: ``--dir`` is the directory batch
    files are released into, ``-r`` the results root, ``-s`` the sample name.
    ``watch`` writes everything under ``<results root>/<sample>/``: per-batch
    work in ``batches/batch_NNNN/`` (numbered from 1) and the cumulative
    ``results/abund_table.tsv``. ``--config``, ``--anns`` and ``-t`` are passed
    through so a replay uses the benchmark's settings — the 16S profile needs its
    own chopper window and databases.

    Start-up is a handshake. ``watch`` only reacts to files created after its
    observer starts, and ignores files already present unless told otherwise,
    while a replay releases batch 1 as soon as ``start`` returns. So the command
    carries ``-v`` (the watcher logs ``Watching <dir> for FASTQ files`` when it
    begins) and ``--scan-existing`` (a file landing during start-up is still
    picked up), and ``start`` blocks until that line appears.

    A batch counts as done when ``watch`` prints ``--- Batch N ---``: its batch
    callback, which runs only after EM and the cumulative tables are written.
    The batch's annotation table appears earlier, before EM, so polling for it
    could score the previous batch's abundances. A batch whose pipeline fails
    logs ``Batch N pipeline failed`` and never prints the line, so that is
    raised at once rather than after the timeout.
    """

    def __init__(
        self,
        sample: str = "replay",
        poll_s: float = 1.0,
        timeout_s: float = 3600.0,
        config: str | Path | None = None,
        anns: str | None = None,
        threads: int | None = None,
        executable: str = "emtrace",
    ):
        self.sample = sample
        self.poll_s = poll_s
        self.timeout_s = timeout_s
        self.config = config
        self.anns = anns
        self.threads = threads
        self.executable = executable
        self._proc: subprocess.Popen | None = None
        self._outdir: Path | None = None

    def command(self, watch_dir: Path, outdir: Path) -> list[str]:
        """The ``emtrace watch`` command line for this replay."""
        cmd = [self.executable, "watch", "--dir", str(watch_dir),
               "-r", str(outdir), "-s", self.sample,
               "--scan-existing", "-v"]
        if self.config:
            cmd += ["--config", str(self.config)]
        if self.anns:
            cmd += ["--anns", self.anns]
        if self.threads:
            cmd += ["-t", str(self.threads)]
        return cmd

    @property
    def run_dir(self) -> Path:
        """Where ``watch`` writes: ``<results root>/<sample>/``."""
        assert self._outdir is not None
        return self._outdir / self.sample

    @property
    def log_path(self) -> Path:
        assert self._outdir is not None
        return self._outdir / "watch.log"

    def start(self, watch_dir: Path, outdir: Path) -> None:
        self._outdir = outdir
        outdir.mkdir(parents=True, exist_ok=True)
        watch_dir.mkdir(parents=True, exist_ok=True)
        env = dict(os.environ, PYTHONUNBUFFERED="1")
        self._proc = subprocess.Popen(
            self.command(watch_dir, outdir),
            stdout=self.log_path.open("w"),
            stderr=subprocess.STDOUT,
            env=env,
        )
        self._await_ready()

    def _await_ready(self, timeout_s: float = 300.0) -> None:
        """Block until the watcher logs that it is watching the directory."""
        deadline = time.monotonic() + timeout_s
        while time.monotonic() < deadline:
            if "for FASTQ files" in self._log_text():
                time.sleep(1.0)  # let the observer thread start after the log line
                return
            if self._proc is not None and self._proc.poll() is not None:
                raise RuntimeError(
                    f"emtrace watch exited (code {self._proc.returncode}) "
                    f"during start-up — see {self.log_path}"
                )
            time.sleep(self.poll_s)
        raise TimeoutError(f"emtrace watch not ready within {timeout_s}s — see {self.log_path}")

    def _batch_dir(self, index: int) -> Path:
        return self.run_dir / "batches" / f"batch_{index:04d}"

    def _log_text(self) -> str:
        try:
            return self.log_path.read_text(errors="replace")
        except FileNotFoundError:
            return ""

    def await_batch(self, index: int) -> None:
        deadline = time.monotonic() + self.timeout_s
        done = re.compile(rf"^--- Batch {index} ---$", re.MULTILINE)
        failed = f"Batch {index} pipeline failed"
        while time.monotonic() < deadline:
            log = self._log_text()
            if done.search(log):
                return
            if failed in log:
                raise RuntimeError(
                    f"emtrace watch: batch {index} pipeline failed — see {self.log_path}"
                )
            if self._proc is not None and self._proc.poll() is not None:
                raise RuntimeError(
                    f"emtrace watch exited (code {self._proc.returncode}) "
                    f"while waiting for batch {index} — see {self.log_path}"
                )
            time.sleep(self.poll_s)
        raise TimeoutError(f"batch {index} not reported within {self.timeout_s}s")

    def abundance(self) -> pd.DataFrame | None:
        for candidate in (
            self.run_dir / "results" / "abund_table.tsv",
            self.run_dir / "results" / "reports" / "abund_table.tsv",
        ):
            if candidate.exists() and candidate.stat().st_size:
                return pd.read_csv(candidate, sep="\t", dtype={"taxid": str})
        return None

    def stop(self) -> None:
        if self._proc is not None and self._proc.poll() is None:
            self._proc.terminate()
            try:
                self._proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                self._proc.kill()


def score_snapshot(
    abund: pd.DataFrame | None,
    cfg: RealtimeConfig,
    truth: pd.DataFrame | None,
    targets: pd.DataFrame | None = None,
) -> dict:
    """Top taxon, EM probability, L1 and set metrics for one batch.

    Set metrics are computed over taxa at or above ``cfg.read_min``: the top
    taxon is read off the unthresholded ranking (it is the pipeline's answer
    either way), but precision would be meaningless if every one-read taxon
    counted as a false positive. With a panel (*targets*) no set metrics are
    computed; ``panel_reads`` holds the reads of each panel taxon instead.
    """
    out: dict = {
        "top_taxon": None, "em_prob": None, "L1": None, "n_truth_taxa": None,
        "n_pred": None, "n_true": None, "TP": None, "FP": None, "FN": None,
        "Pr": None, "Rc": None, "F1": None, "panel_reads": {},
    }
    if abund is None or abund.empty or cfg.rank not in abund.columns:
        return out

    df = abund.copy()
    df[cfg.score_col] = pd.to_numeric(df.get(cfg.score_col), errors="coerce").fillna(0)
    # Optional tail trim, off by default. Deliberately ADDITIVE: it drops rows
    # before the per-rank sum, and `read_min` below still thresholds the sums
    # afterwards. The two filters are different questions -- "is this taxon
    # part of the library's bulk?" and "does this rank have enough reads to
    # report?" -- and answer them at different granularities.
    if cfg.tail > 0:
        df = tail_trim(
            df, cfg.tail,
            resolve_count_col(df, prefer=cfg.tail_col, fallback=cfg.score_col),
        )
    ranked = df.groupby(cfg.rank, dropna=True)[cfg.score_col].sum().sort_values(ascending=False)
    ranked = ranked[ranked.index.astype(str).str.len() > 0]
    if not ranked.empty:
        out["top_taxon"] = str(ranked.index[0])
        if cfg.prob_col in df.columns:
            prob = pd.to_numeric(df[cfg.prob_col], errors="coerce")
            top_rows = df[df[cfg.rank].astype(str) == out["top_taxon"]]
            if not top_rows.empty:
                out["em_prob"] = float(
                    pd.to_numeric(top_rows[cfg.prob_col], errors="coerce").max()
                ) if prob.notna().any() else None

    if truth is not None and not truth.empty:
        from ngs_benchmark.detection.abundance import abundance_metrics

        try:
            m = abundance_metrics(df, truth, ranks=(cfg.rank,), pred_weight=cfg.score_col)
            if not m.empty:
                out["L1"] = float(m.iloc[0]["L1"])
        except (ValueError, KeyError):
            pass
        # Only taxa actually present count as truth: the 16S table carries the two
        # yeasts at 0 % on the pct_16S_only basis, which would otherwise score as
        # permanent false negatives (L1 above already ignores them).
        present = truth
        if TRUTH_WEIGHT in truth.columns:
            present = truth[pd.to_numeric(truth[TRUTH_WEIGHT], errors="coerce").fillna(0) > 0]
        truth_taxa = set(present[cfg.rank].dropna().astype(str)) if cfg.rank in present.columns else set()
        out["n_truth_taxa"] = len(truth_taxa & set(ranked.index.astype(str)))
    elif cfg.truth_taxon:
        truth_taxa = {cfg.truth_taxon}
    else:
        truth_taxa = set()

    if truth_taxa:
        pred_taxa = set(ranked[ranked >= cfg.read_min].index.astype(str))
        out.update(set_scores(pred_taxa, truth_taxa))

    if targets is not None:
        for rank, grp in targets.groupby("rank"):
            if rank not in df.columns:
                continue
            sums = df.groupby(rank, dropna=True)[cfg.score_col].sum()
            for taxon in grp["taxon"]:
                out["panel_reads"][taxon] = float(sums.get(taxon, 0.0))

    return out


def _check_scoring(cfg: RealtimeConfig) -> pd.DataFrame | None:
    """Load the panel, refusing it beside a truth: the two schemes exclude each other."""
    if not cfg.panel:
        return None
    if cfg.truth_taxon or cfg.truth_table:
        raise ValueError(
            "--panel scores a field library without a truth; drop --truth-taxon / --truth"
        )
    return load_panel(cfg.panel, cfg.rank)


def _tables(
    snapshots: list[Snapshot], cfg: RealtimeConfig, targets: pd.DataFrame | None
) -> dict[str, pd.DataFrame]:
    tables = {
        "timeseries": timeseries(snapshots, cfg.library),
        "latency": latency(snapshots, cfg.library),
        "accuracy": accuracy(
            snapshots, cfg.library, truth_taxon=cfg.truth_taxon,
            stabilise_on=cfg.stabilise_on, tol=cfg.stabilise_tol,
        ),
    }
    if cfg.truth_taxon:
        tables["first_call"] = first_call(snapshots, cfg.truth_taxon, cfg.library)
        tables["first_detection"] = first_detection(snapshots, cfg.library)
    if targets is not None:
        tables["panel"] = panel(snapshots, targets, cfg.library, cfg.read_min)
        tables["panel_summary"] = panel_summary(tables["panel"], cfg.library)
    return tables


def run(cfg: RealtimeConfig, driver: PipelineDriver, outdir: str | Path) -> dict[str, pd.DataFrame]:
    """Replay one library and write the three Module C tables into *outdir*.

    Returns the tables as well, so a caller scoring several libraries can
    concatenate them without re-reading from disk.
    """
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    batches = scan_library(cfg.fastq_dir, reads_only=cfg.reads_only)
    gaps = contiguity_gaps(batches, max_gap_h=cfg.max_gap_h)
    if gaps:
        first, hours = gaps[0]
        raise ValueError(
            f"{cfg.library}: sequencing time jumps {hours:.1f} h at batch {first} — "
            "this looks like two runs combined into one directory, which cannot be "
            "replayed as a single run. Split it or choose a contiguous library."
        )

    targets = _check_scoring(cfg)
    truth = None
    if cfg.truth_table:
        truth = pd.read_csv(cfg.truth_table, sep="\t")

    watch_dir = out / "watch"
    pipeline_out = out / "pipeline"
    abund_dir = out / "abund"
    abund_dir.mkdir(exist_ok=True)
    driver.start(watch_dir, pipeline_out)

    snapshots: list[Snapshot] = []
    try:
        for rel in replay_batches(
            batches, watch_dir, driver.await_batch, copy=not cfg.hardlink
        ):
            abund = driver.abundance()
            if abund is not None:
                write_tsv(abund, abund_dir / f"batch_{rel.batch.index:04d}.tsv")
            scored = score_snapshot(abund, cfg, truth, targets)
            snapshots.append(
                Snapshot(
                    batch=rel.batch.index,
                    seq_time_s=rel.batch.seq_time_s,
                    wall_s=round(rel.wall_s, 3),
                    reads_cum=rel.reads_cum,
                    latency_s=round(rel.latency_s, 3),
                    **scored,
                )
            )
    finally:
        driver.stop()

    tables = _tables(snapshots, cfg, targets)
    for name, df in tables.items():
        write_tsv(df, out / f"{name}.tsv")
    return tables


def rescore(cfg: RealtimeConfig, outdir: str | Path) -> dict[str, pd.DataFrame]:
    """Score a finished replay again, from its own outputs, without replaying.

    The clocks and read counts come from the replay's ``timeseries.tsv``, the
    per-batch abundances from ``abund/batch_NNNN.tsv`` (written by :func:`run`,
    or back-filled for replays made before it kept them). A batch without an
    abundance table scores as an empty one, as it did live. The tables are
    rewritten in place; tables the new scoring does not produce are removed, so
    the directory never mixes two schemes.
    """
    out = Path(outdir)
    targets = _check_scoring(cfg)
    truth = pd.read_csv(cfg.truth_table, sep="\t") if cfg.truth_table else None
    abund_dir = out / "abund"
    if not abund_dir.is_dir():
        raise ValueError(f"{out}: no per-batch abundance tables (abund/) to rescore from")

    ts = pd.read_csv(out / "timeseries.tsv", sep="\t")
    snapshots: list[Snapshot] = []
    for r in ts.itertuples(index=False):
        path = abund_dir / f"batch_{int(r.batch):04d}.tsv"
        abund = (pd.read_csv(path, sep="\t", dtype={"taxid": str})
                 if path.is_file() and path.stat().st_size else None)
        snapshots.append(
            Snapshot(
                batch=int(r.batch),
                seq_time_s=None if pd.isna(r.seq_time_s) else float(r.seq_time_s),
                wall_s=float(r.wall_s),
                reads_cum=int(r.reads_cum),
                latency_s=float(r.latency_s),
                **score_snapshot(abund, cfg, truth, targets),
            )
        )

    tables = _tables(snapshots, cfg, targets)
    for name in ("first_call", "first_detection", "panel", "panel_summary"):
        if name not in tables:
            (out / f"{name}.tsv").unlink(missing_ok=True)
    for name, df in tables.items():
        write_tsv(df, out / f"{name}.tsv")
    return tables
