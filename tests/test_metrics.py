"""Unit tests for §3.2 set metrics."""

from __future__ import annotations

import pandas as pd

from ngs_benchmark.detection.metrics import per_taxon_table, set_metrics


def test_set_metrics_simple(pred_simple, truth_simple) -> None:
    df = set_metrics(pred_simple, truth_simple, ranks=("species",), taxgroups=("Viruses",))
    row = df.iloc[0]
    # Viruses: pred {vA, vB, vE, vF}, truth {vA, vB, vC} → TP=2, FP=2, FN=1.
    assert row["TP"] == 2
    assert row["FP"] == 2
    assert row["FN"] == 1
    assert row["Pr"] == 0.5      # 2/(2+2)
    assert row["Rc"] == 2 / 3    # 2/(2+1)
    assert row["F1"] == 2 * 0.5 * (2 / 3) / (0.5 + 2 / 3)
    assert row["mean_accuracy"] == (0.5 + 2 / 3) / 2


def test_set_metrics_bacteria_perfect(pred_simple, truth_simple) -> None:
    df = set_metrics(pred_simple, truth_simple, ranks=("species",), taxgroups=("Bacteria",))
    row = df.iloc[0]
    assert row["TP"] == 1
    assert row["FP"] == 0
    assert row["FN"] == 0
    assert row["Pr"] == 1.0
    assert row["Rc"] == 1.0
    assert row["F1"] == 1.0


def test_set_metrics_all_mismatch() -> None:
    pred = pd.DataFrame(
        {"readn": [10.0], "taxid": [1], "species": ["x"], "genus": ["g"], "family": ["f"],
         "order": ["o"], "class": ["c"], "division": ["Viruses"]}
    )
    truth = pd.DataFrame(
        {"readn": [10.0], "taxid": [2], "species": ["y"], "genus": ["g"], "family": ["f"],
         "order": ["o"], "class": ["c"], "division": ["Viruses"]}
    )
    df = set_metrics(pred, truth, ranks=("species",), taxgroups=("Viruses",))
    assert df.iloc[0]["TP"] == 0
    assert df.iloc[0]["FP"] == 1
    assert df.iloc[0]["FN"] == 1
    assert df.iloc[0]["Pr"] == 0.0
    assert df.iloc[0]["Rc"] == 0.0
    assert df.iloc[0]["F1"] == 0.0


def test_set_metrics_empty_pred(truth_simple) -> None:
    empty = truth_simple.iloc[0:0].copy()
    df = set_metrics(empty, truth_simple, ranks=("species",), taxgroups=("Viruses",))
    assert df.iloc[0]["TP"] == 0
    assert df.iloc[0]["FP"] == 0
    assert df.iloc[0]["FN"] == 3
    assert df.iloc[0]["Rc"] == 0.0
    # Pr = 0/0 → 0 by convention.
    assert df.iloc[0]["Pr"] == 0.0


def test_per_taxon_table_statuses(pred_simple, truth_simple) -> None:
    df = per_taxon_table(pred_simple, truth_simple, ranks=("species",), taxgroups=("Viruses",))
    statuses = dict(zip(df["taxon"], df["status"]))
    assert statuses["vA"] == "TP"
    assert statuses["vB"] == "TP"
    assert statuses["vC"] == "FN"
    assert statuses["vE"] == "FP"
    assert statuses["vF"] == "FP"
