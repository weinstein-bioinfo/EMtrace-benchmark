"""Tests for the legacy-R `summary.xlsx` writer."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from openpyxl import load_workbook

from ngs_benchmark.detection.metrics import set_metrics
from ngs_benchmark.io.writers import write_xlsx_summary


def test_xlsx_summary_has_expected_tabs(pred_simple, truth_simple, tmp_path: Path) -> None:
    metrics_df = set_metrics(
        pred_simple, truth_simple, ranks=("species",), taxgroups=("Viruses", "Bacteria")
    )
    out = tmp_path / "summary.xlsx"
    write_xlsx_summary(
        metrics_df,
        pred_simple,
        truth_simple,
        path=out,
        taxgroups=("Viruses", "Bacteria"),
    )
    wb = load_workbook(out)
    assert wb.sheetnames[0] == "SUMMARY"
    # Each taxgroup gets 5 detail tabs: true, pred, TP, FN, FP.
    for short in ("vi", "ba"):
        for tab in ("true", "pred", "TP", "FN", "FP"):
            assert f"{short}_{tab}" in wb.sheetnames


def test_xlsx_summary_tab_contents(pred_simple, truth_simple, tmp_path: Path) -> None:
    metrics_df = set_metrics(
        pred_simple, truth_simple, ranks=("species",), taxgroups=("Viruses",)
    )
    out = tmp_path / "summary.xlsx"
    write_xlsx_summary(
        metrics_df, pred_simple, truth_simple, path=out, taxgroups=("Viruses",)
    )
    wb = load_workbook(out)

    # SUMMARY tab header should include TP, FP, FN columns.
    ws = wb["SUMMARY"]
    header = [c.value for c in ws[1]]
    for col in ("TP", "FP", "FN", "Pr", "Rc", "F1", "mean_accuracy"):
        assert col in header

    # vi_TP must list the two true-positive viruses (vA, vB).
    ws_tp = wb["vi_TP"]
    species_col = [c.value for c in ws_tp[1]].index("species")
    species_in_tp = {
        ws_tp.cell(row=r, column=species_col + 1).value for r in range(2, ws_tp.max_row + 1)
    }
    assert species_in_tp == {"vA", "vB"}

    # vi_FP must list the two off-target viruses (vE, vF).
    ws_fp = wb["vi_FP"]
    species_col = [c.value for c in ws_fp[1]].index("species")
    species_in_fp = {
        ws_fp.cell(row=r, column=species_col + 1).value for r in range(2, ws_fp.max_row + 1)
    }
    assert species_in_fp == {"vE", "vF"}

    # vi_FN must list the missed virus vC.
    ws_fn = wb["vi_FN"]
    species_col = [c.value for c in ws_fn[1]].index("species")
    species_in_fn = {
        ws_fn.cell(row=r, column=species_col + 1).value for r in range(2, ws_fn.max_row + 1)
    }
    assert species_in_fn == {"vC"}


def test_xlsx_fp_rows_red_filled(pred_simple, truth_simple, tmp_path: Path) -> None:
    metrics_df = set_metrics(
        pred_simple, truth_simple, ranks=("species",), taxgroups=("Viruses",)
    )
    out = tmp_path / "summary.xlsx"
    write_xlsx_summary(
        metrics_df, pred_simple, truth_simple, path=out, taxgroups=("Viruses",)
    )
    wb = load_workbook(out)
    ws_pred = wb["vi_pred"]
    species_col = [c.value for c in ws_pred[1]].index("species") + 1
    fill_status: dict[str, str | None] = {}
    for r in range(2, ws_pred.max_row + 1):
        species = ws_pred.cell(row=r, column=species_col).value
        fill = ws_pred.cell(row=r, column=species_col).fill
        fill_status[species] = fill.fgColor.rgb if fill and fill.fgColor else None
    # vA, vB are TPs ⇒ unfilled. vE, vF are FPs ⇒ red.
    assert fill_status["vA"] in (None, "00000000")
    assert fill_status["vB"] in (None, "00000000")
    assert fill_status["vE"] == "00FF5050"
    assert fill_status["vF"] == "00FF5050"


def test_xlsx_fn_rows_red_filled_in_true_tab(pred_simple, truth_simple, tmp_path: Path) -> None:
    metrics_df = set_metrics(
        pred_simple, truth_simple, ranks=("species",), taxgroups=("Viruses",)
    )
    out = tmp_path / "summary.xlsx"
    write_xlsx_summary(
        metrics_df, pred_simple, truth_simple, path=out, taxgroups=("Viruses",)
    )
    wb = load_workbook(out)
    ws_true = wb["vi_true"]
    species_col = [c.value for c in ws_true[1]].index("species") + 1
    for r in range(2, ws_true.max_row + 1):
        species = ws_true.cell(row=r, column=species_col).value
        fill = ws_true.cell(row=r, column=species_col).fill
        rgb = fill.fgColor.rgb if fill and fill.fgColor else None
        if species == "vC":
            assert rgb == "00FF5050"
        else:
            assert rgb in (None, "00000000")


def test_xlsx_percent_format_on_summary(pred_simple, truth_simple, tmp_path: Path) -> None:
    metrics_df = set_metrics(
        pred_simple, truth_simple, ranks=("species",), taxgroups=("Viruses",)
    )
    out = tmp_path / "summary.xlsx"
    write_xlsx_summary(
        metrics_df, pred_simple, truth_simple, path=out, taxgroups=("Viruses",)
    )
    wb = load_workbook(out)
    ws = wb["SUMMARY"]
    header = [c.value for c in ws[1]]
    pr_idx = header.index("Pr") + 1
    # Row 2 (first data row) should have percent format on Pr.
    assert ws.cell(row=2, column=pr_idx).number_format == "0.0%"


def test_xlsx_empty_pred(truth_simple, tmp_path: Path) -> None:
    empty = truth_simple.iloc[0:0].copy()
    metrics_df = set_metrics(
        empty, truth_simple, ranks=("species",), taxgroups=("Viruses",)
    )
    out = tmp_path / "summary.xlsx"
    # Should write without crashing even when pred is empty.
    write_xlsx_summary(
        metrics_df, empty, truth_simple, path=out, taxgroups=("Viruses",)
    )
    wb = load_workbook(out)
    assert "SUMMARY" in wb.sheetnames
    assert "vi_pred" in wb.sheetnames
