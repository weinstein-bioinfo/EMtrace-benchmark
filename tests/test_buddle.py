"""Unit tests for §3.5 Buddle PMR/RPMR rule."""

from __future__ import annotations

import math
import pandas as pd

from ngs_benchmark.detection.buddle import annotate, filter_by_rule


def _df(species: list[str], readn: list[float], division: list[str] | None = None) -> pd.DataFrame:
    n = len(species)
    return pd.DataFrame(
        {
            "readn": readn,
            "taxid": list(range(1, n + 1)),
            "species": species,
            "genus": ["g"] * n,
            "family": ["f"] * n,
            "order": ["o"] * n,
            "class": ["c"] * n,
            "division": division if division is not None else ["Viruses"] * n,
        }
    )


def test_pmr_threshold_excludes_low_share_virus() -> None:
    # Two viruses: vA dominates; vB at 0.005% of microbial reads.
    # PMR threshold for viruses is 0.01% — vB must fail.
    # readn vA=99995, vB=5 → total=100000 → vB PMR = 0.005%.
    sample = _df(["vA", "vB"], [99995.0, 5.0])
    nc = _df([], [])
    out = annotate(sample, nc, total_reads_sample=1e6, total_reads_nc=1e6)
    pass_map = dict(zip(out["species"], out["buddle_pass"]))
    assert pass_map["vA"]
    assert not pass_map["vB"]


def test_pmr_threshold_admits_moderate_share() -> None:
    # vA at 80 reads, vB at 0.02% (above 0.01% threshold) with empty NC ⇒ pass.
    sample = _df(["vA", "vB"], [80, 0.020 * 80 / (100 - 0.020) + 0.001])
    # Easier construction: directly give 0.02% share.
    total = 1000
    sample = pd.DataFrame(
        {
            "readn": [total * (1 - 0.0002), total * 0.0002],
            "taxid": [1, 2],
            "species": ["vA", "vB"],
            "genus": ["g", "g"],
            "family": ["f", "f"],
            "order": ["o", "o"],
            "class": ["c", "c"],
            "division": ["Viruses", "Viruses"],
        }
    )
    nc = _df([], [])
    out = annotate(sample, nc, total_reads_sample=1e6, total_reads_nc=1e6)
    pass_map = dict(zip(out["species"], out["buddle_pass"]))
    assert pass_map["vA"]
    assert pass_map["vB"]


def test_rpmr_inf_when_nc_zero() -> None:
    sample = _df(["vA"], [100])
    nc = _df([], [])
    out = annotate(sample, nc, total_reads_sample=1e6, total_reads_nc=1e6)
    assert math.isinf(out["RPMR"].iloc[0])
    assert out["buddle_pass"].iloc[0]


def test_rpmr_finite_when_nc_present() -> None:
    sample = _df(["vA"], [100])
    nc = _df(["vA"], [50])
    out = annotate(sample, nc, total_reads_sample=1e6, total_reads_nc=1e6)
    # RPM_sample/RPM_nc = 100/50 = 2 — fails virus threshold (>5).
    assert out["RPMR"].iloc[0] == 2.0
    # NC nonzero AND RPMR=2 < 5 ⇒ fails.
    assert not out["buddle_pass"].iloc[0]


def test_filter_by_rule_drops_failed() -> None:
    sample = _df(["vA", "vB"], [80, 0.001])
    nc = _df([], [])
    annotated = annotate(sample, nc, total_reads_sample=1e6, total_reads_nc=1e6)
    out = filter_by_rule(annotated)
    assert set(out["species"]) == {"vA"}


def test_other_division_uses_higher_threshold() -> None:
    # A bacterium at PMR = 0.5% → fails (other_pmr threshold = 1%).
    sample = pd.DataFrame(
        {
            "readn": [995.0, 5.0],
            "taxid": [1, 2],
            "species": ["bA", "bB"],
            "genus": ["g", "g"],
            "family": ["f", "f"],
            "order": ["o", "o"],
            "class": ["c", "c"],
            "division": ["Bacteria", "Bacteria"],
        }
    )
    nc = _df([], [], division=[])
    out = annotate(sample, nc, total_reads_sample=1e6, total_reads_nc=1e6)
    pass_map = dict(zip(out["species"], out["buddle_pass"]))
    # bA = 99.5% PMR, but RPMR = inf (NC empty) → 99.5%>1%, but rule requires RPMR>10
    # AND PMR>1% for non-viruses. RPMR=inf satisfies; PMR=99.5%>1% satisfies → pass.
    assert pass_map["bA"]
    # bB = 0.5% PMR < 1% → fails.
    assert not pass_map["bB"]
