"""Unit tests for §3.3 PR curve + AUPR."""

from __future__ import annotations

import pandas as pd
import pytest

from ngs_benchmark.detection.curves import pr_curve


def _viral_df(species: list[str], scores: list[float]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "readn": scores,
            "taxid": list(range(1, len(species) + 1)),
            "species": species,
            "genus": ["g"] * len(species),
            "family": ["f"] * len(species),
            "order": ["o"] * len(species),
            "class": ["c"] * len(species),
            "division": ["Viruses"] * len(species),
        }
    )


def test_perfect_ranking_aupr_approaches_one() -> None:
    # Perfect ranking with K=3 positives. The PR curve goes (0,0)→(1/3,1)→(2/3,1)→(1,1)
    # → then degrades as FPs are added. Trapezoidal AUPR = 1 - 1/(2K) = 5/6.
    # (sklearn's average_precision uses a right-step rule and would return 1.0;
    # the trapezoidal convention from ROADMAP §3.3 under-estimates by 1/(2K).)
    # F-max is still 1.0 because we hit (Pr=1, Rc=1) at threshold 8.
    pred = _viral_df(["t1", "t2", "t3", "fp1", "fp2"], [10, 9, 8, 1, 0.5])
    truth = _viral_df(["t1", "t2", "t3"], [1, 1, 1])
    res = pr_curve(pred, truth, rank="species", taxgroup="Viruses", score_col="readn")
    assert res.aupr == pytest.approx(5 / 6, abs=1e-9)
    assert res.f_max == pytest.approx(1.0, abs=1e-9)


def test_no_true_positives_gives_aupr_zero() -> None:
    pred = _viral_df(["fp1", "fp2"], [10, 5])
    truth = _viral_df(["t1"], [1])
    res = pr_curve(pred, truth, rank="species", taxgroup="Viruses", score_col="readn")
    assert res.aupr == pytest.approx(0.0, abs=1e-9)
    assert res.f_max == 0.0


def test_known_asymmetric_aupr() -> None:
    # truth = {t1, t2}. pred ordered: t1(10), fp1(8), t2(5), fp2(1)
    # Sweep operating points (score >= t):
    #  t=+inf:  P=0, R=0
    #  t=10:    pred={t1}            TP=1 FP=0 FN=1 → Pr=1.0  Rc=0.5
    #  t=8:     pred={t1,fp1}        TP=1 FP=1 FN=1 → Pr=0.5  Rc=0.5
    #  t=5:     pred={t1,fp1,t2}     TP=2 FP=1 FN=0 → Pr=2/3  Rc=1.0
    #  t=1:     pred={all}           TP=2 FP=2 FN=0 → Pr=0.5  Rc=1.0
    # Sorted by recall asc: (0,0), (0.5,1.0), (0.5,0.5), (1.0,2/3), (1.0,0.5)
    # trapezoidal area = sum of trapezoids between consecutive (Rc, Pr) points.
    # By hand: (R, P) sorted asc by R: [(0,0), (0.5,1.0), (0.5,0.5), (1.0,2/3), (1.0,0.5)]
    # Between (0,0)→(0.5,1.0):      width 0.5, mean P = 0.5  → 0.25
    # Between (0.5,1.0)→(0.5,0.5):  width 0   → 0
    # Between (0.5,0.5)→(1.0,2/3):  width 0.5, mean P = (0.5+2/3)/2 = 7/12 → 7/24
    # Between (1.0,2/3)→(1.0,0.5):  width 0  → 0
    # total = 0.25 + 7/24 = 6/24 + 7/24 = 13/24 ≈ 0.5417
    pred = _viral_df(["t1", "fp1", "t2", "fp2"], [10, 8, 5, 1])
    truth = _viral_df(["t1", "t2"], [1, 1])
    res = pr_curve(pred, truth, rank="species", taxgroup="Viruses", score_col="readn")
    assert res.aupr == pytest.approx(13 / 24, abs=1e-6)


def test_fmax_operating_point() -> None:
    # F1 max occurs when Pr=1, Rc=1 → at threshold t=5 with pred={t1,t2}.
    pred = _viral_df(["t1", "t2", "fp1"], [10, 5, 1])
    truth = _viral_df(["t1", "t2"], [1, 1])
    res = pr_curve(pred, truth, rank="species", taxgroup="Viruses", score_col="readn")
    assert res.f_max == pytest.approx(1.0)
    assert res.pr_fmax == pytest.approx(1.0)
    assert res.rc_fmax == pytest.approx(1.0)
    assert res.threshold_at_fmax == 5.0


def test_uses_score_col_kwarg() -> None:
    # The same set of predictions but with non-monotonic readn vs prob — using
    # prob as score should reorder.
    pred = _viral_df(["t1", "t2", "fp1"], [1, 10, 5])
    pred["prob"] = [0.99, 0.1, 0.5]
    truth = _viral_df(["t1", "t2"], [1, 1])
    res = pr_curve(pred, truth, rank="species", taxgroup="Viruses", score_col="prob")
    # By prob: t1(0.99), fp1(0.5), t2(0.1).
    # Threshold sweep gives  (Rc=0.5,Pr=1.0), (Rc=0.5,Pr=0.5), (Rc=1.0,Pr=2/3).
    # F-max = 2*1*0.5/(1.5) = 2/3 ... but at (Rc=1, Pr=2/3) F1 = 2*2/3/(5/3) = 4/5 = 0.8.
    assert res.f_max == pytest.approx(0.8, abs=1e-9)
