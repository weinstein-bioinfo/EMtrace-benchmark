"""Adapter / IO loader tests."""

from __future__ import annotations

from pathlib import Path
import pandas as pd
import pytest

from ngs_benchmark.io.adapters import get_adapter, available
from ngs_benchmark.io.adapters import truth as truth_io


def test_emtrace_adapter_roundtrips(minimal_pred_tsv: Path) -> None:
    df = get_adapter("emtrace")(minimal_pred_tsv)
    assert {"readn", "taxid", "species", "genus", "family", "division"}.issubset(df.columns)
    assert len(df) == 5


def test_emtrace_adapter_passes_prob_column(minimal_pred_tsv: Path) -> None:
    df = get_adapter("emtrace")(minimal_pred_tsv)
    assert "prob" in df.columns


def test_emtrace_adapter_renames_superkingdom(tmp_path: Path) -> None:
    p = tmp_path / "alt.tsv"
    pd.DataFrame(
        {
            "readn": [10.0],
            "taxid": [1],
            "species": ["s"],
            "genus": ["g"],
            "family": ["f"],
            "order": ["o"],
            "class": ["c"],
            "superkingdom": ["Viruses"],
        }
    ).to_csv(p, sep="\t", index=False)
    df = get_adapter("emtrace")(p)
    assert "division" in df.columns
    assert df["division"].iloc[0] == "Viruses"


def test_truth_loader_flags_membership(minimal_truth_tsv: Path) -> None:
    df = truth_io.load(minimal_truth_tsv)
    assert df.attrs["truth_flavour"] == "membership"


def test_truth_loader_flags_dilution(tmp_path: Path) -> None:
    p = tmp_path / "buddle_truth.tsv"
    pd.DataFrame(
        {
            "sample_id": ["s1"],
            "gc_per_ml": [60.0],
            "nc_sample_id": ["nc1"],
            "taxid": [1],
            "species": ["vA"],
            "genus": ["g"],
            "family": ["f"],
            "order": ["o"],
            "class": ["c"],
            "division": ["Viruses"],
        }
    ).to_csv(p, sep="\t", index=False)
    df = truth_io.load(p)
    assert df.attrs["truth_flavour"] == "dilution"


def test_unknown_adapter_errors() -> None:
    with pytest.raises(ValueError, match="Unknown adapter"):
        get_adapter("not-a-real-adapter")


def test_available_lists_all_adapters() -> None:
    names = set(available())
    assert {"emtrace", "czid", "epi2me", "kraken_bracken", "megan_lr", "metamix"}.issubset(names)
