"""Unit tests for §3.6 dilution series."""

from __future__ import annotations

import pandas as pd
import pytest

from ngs_benchmark.detection.dilution import SamplePred, dilution_metrics


def _viral_df(species: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "readn": [10.0] * len(species),
            "taxid": list(range(1, len(species) + 1)),
            "species": species,
            "genus": ["g"] * len(species),
            "family": ["f"] * len(species),
            "order": ["o"] * len(species),
            "class": ["c"] * len(species),
            "division": ["Viruses"] * len(species),
        }
    )


def _truth() -> pd.DataFrame:
    # Six mock-community viruses (Buddle-like).
    return _viral_df(["v1", "v2", "v3", "v4", "v5", "v6"])


def test_dilution_sensitivity_curve() -> None:
    # 60 gc/mL: detect 1; 600 gc/mL: detect 3; 6000 gc/mL: detect all 6.
    samples = [
        SamplePred("s60", 60.0, None, _viral_df(["v1"])),
        SamplePred("s600", 600.0, None, _viral_df(["v1", "v2", "v3"])),
        SamplePred("s6000", 6000.0, None, _viral_df(["v1", "v2", "v3", "v4", "v5", "v6"])),
    ]
    df, summary = dilution_metrics(samples, _truth(), rank="species", taxgroup="Viruses")
    rows = {r["gc_per_ml"]: r for r in df.to_dict(orient="records")}
    assert rows[60.0]["sensitivity"] == 1 / 6
    assert rows[600.0]["sensitivity"] == 3 / 6
    assert rows[6000.0]["sensitivity"] == 1.0
    assert summary["total_unique_fps"] == 0


def test_dilution_fp_counts() -> None:
    samples = [
        SamplePred("s60", 60.0, None, _viral_df(["v1", "off1"])),
        SamplePred("s600", 600.0, None, _viral_df(["v1", "v2", "off1", "off2"])),
    ]
    df, summary = dilution_metrics(samples, _truth(), rank="species", taxgroup="Viruses")
    rows = {r["gc_per_ml"]: r for r in df.to_dict(orient="records")}
    assert rows[60.0]["FP"] == 1
    assert rows[60.0]["fp_unique_taxa"] == 1
    assert rows[600.0]["FP"] == 2
    assert rows[600.0]["fp_unique_taxa"] == 2
    # off1 appears at both levels; off2 only at one; total distinct = 2.
    assert summary["total_unique_fps"] == 2


def test_fp_persistence_with_matrix_id() -> None:
    # Both samples share matrix M; off1 appears at both → persistent.
    # off2 appears at one only → not persistent.
    samples = [
        SamplePred("s60", 60.0, "M", _viral_df(["v1", "off1"])),
        SamplePred("s600", 600.0, "M", _viral_df(["v1", "off1", "off2"])),
    ]
    df, _ = dilution_metrics(samples, _truth(), rank="species", taxgroup="Viruses")
    rows = {r["gc_per_ml"]: r for r in df.to_dict(orient="records")}
    assert rows[60.0]["fp_persistent_taxa"] == 1   # off1
    assert rows[600.0]["fp_persistent_taxa"] == 1  # off1


def test_fp_persistence_na_when_no_matrix_id() -> None:
    samples = [
        SamplePred("s60", 60.0, None, _viral_df(["v1", "off1"])),
    ]
    df, summary = dilution_metrics(samples, _truth(), rank="species", taxgroup="Viruses")
    assert df["fp_persistent_taxa"].iloc[0] is None
    assert "fp_persistent_taxa_warning" in summary


def test_low_n_flag() -> None:
    samples = [
        SamplePred("s60a", 60.0, None, _viral_df(["v1"])),
        SamplePred("s60b", 60.0, None, _viral_df(["v2"])),
    ]
    df, _ = dilution_metrics(samples, _truth(), rank="species", taxgroup="Viruses")
    assert df["low_n"].iloc[0]
