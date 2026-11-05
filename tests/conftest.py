"""Shared fixtures: small in-memory DataFrames + on-disk minimal TSVs."""

from __future__ import annotations

from pathlib import Path
import pandas as pd
import pytest


# Six taxa: A, B, C real viruses; D real bacterium; E, F off-target.
# Used across unit tests.
@pytest.fixture
def truth_simple() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "readn": [100, 100, 100, 1000],
            "taxid": [1, 2, 3, 100],
            "species": ["vA", "vB", "vC", "bD"],
            "genus": ["GA", "GB", "GC", "BG"],
            "family": ["FA", "FB", "FC", "BF"],
            "order": ["OV", "OV", "OV", "OB"],
            "class": ["CV", "CV", "CV", "CB"],
            "division": ["Viruses", "Viruses", "Viruses", "Bacteria"],
        }
    )


@pytest.fixture
def pred_simple() -> pd.DataFrame:
    """vA, vB hit (TP); vE, vF off-target (FP); vC missed (FN); bD present (TP)."""
    return pd.DataFrame(
        {
            "readn": [50, 40, 30, 5, 500],
            "taxid": [1, 2, 5, 6, 100],
            "species": ["vA", "vB", "vE", "vF", "bD"],
            "genus": ["GA", "GB", "GE", "GF", "BG"],
            "family": ["FA", "FB", "FE", "FF", "BF"],
            "order": ["OV", "OV", "OV", "OV", "OB"],
            "class": ["CV", "CV", "CV", "CV", "CB"],
            "division": ["Viruses", "Viruses", "Viruses", "Viruses", "Bacteria"],
            "prob": [0.9, 0.8, 0.5, 0.1, 0.95],
        }
    )


@pytest.fixture
def minimal_pred_tsv(tmp_path: Path, pred_simple: pd.DataFrame) -> Path:
    p = tmp_path / "pred.tsv"
    pred_simple.to_csv(p, sep="\t", index=False)
    return p


@pytest.fixture
def minimal_truth_tsv(tmp_path: Path, truth_simple: pd.DataFrame) -> Path:
    p = tmp_path / "truth.tsv"
    truth_simple.to_csv(p, sep="\t", index=False)
    return p
