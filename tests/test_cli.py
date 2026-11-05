"""End-to-end CLI tests using click.testing.CliRunner."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from ngs_benchmark.cli import cli


def test_cli_detection_e2e(minimal_pred_tsv: Path, minimal_truth_tsv: Path, tmp_path: Path) -> None:
    outdir = tmp_path / "out"
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "detection",
            "--pred", str(minimal_pred_tsv),
            "--truth", str(minimal_truth_tsv),
            "--outdir", str(outdir),
            "--ranks", "species",
            "--taxgroup", "Viruses",
            "--tail-ba", "0",
        ],
    )
    assert result.exit_code == 0, result.output

    # All output files exist.
    assert (outdir / "metrics.tsv").exists()
    assert (outdir / "prc.tsv").exists()
    assert (outdir / "per_taxon.tsv").exists()
    assert (outdir / "summary.json").exists()

    metrics = pd.read_csv(outdir / "metrics.tsv", sep="\t")
    row = metrics.iloc[0]
    assert row["TP"] == 2
    assert row["FP"] == 2
    assert row["FN"] == 1

    summary = json.loads((outdir / "summary.json").read_text())
    assert "AUPR" in summary
    assert summary["primary_taxgroup"] == "Viruses"


def test_cli_buddle_requires_nc_pred(minimal_pred_tsv: Path, minimal_truth_tsv: Path, tmp_path: Path) -> None:
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "detection",
            "--pred", str(minimal_pred_tsv),
            "--truth", str(minimal_truth_tsv),
            "--outdir", str(tmp_path / "out"),
            "--buddle",
        ],
    )
    assert result.exit_code != 0
    msg = str(result.exception).lower() if result.exception else result.output.lower()
    assert "nc-pred" in msg or "nc_pred" in msg


def test_cli_tn_universe_not_implemented(minimal_pred_tsv: Path, minimal_truth_tsv: Path, tmp_path: Path) -> None:
    tn_file = tmp_path / "tn.tsv"
    tn_file.write_text("species\n")
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "detection",
            "--pred", str(minimal_pred_tsv),
            "--truth", str(minimal_truth_tsv),
            "--outdir", str(tmp_path / "out"),
            "--tn-universe", str(tn_file),
        ],
    )
    assert result.exit_code != 0
    msg = str(result.exception).lower() if result.exception else result.output.lower()
    assert "deferred" in msg or "beta" in msg


def test_cli_unknown_adapter_errors(minimal_pred_tsv: Path, minimal_truth_tsv: Path, tmp_path: Path) -> None:
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "detection",
            "--pred", str(minimal_pred_tsv),
            "--truth", str(minimal_truth_tsv),
            "--outdir", str(tmp_path / "out"),
            "--pred-fmt", "not-a-real-adapter",
        ],
    )
    assert result.exit_code != 0
    msg = str(result.exception).lower() if result.exception else result.output.lower()
    assert "unknown adapter" in msg
