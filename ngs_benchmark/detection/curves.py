"""§3.3 Precision–Recall curve and AUPR.

For each distinct value `t` of the score column in `pred`, recompute
(Pr_t, Rc_t) using the predictions with score >= t. AUPR is trapezoidal
integration over recall.

Also returned:
    F_max               — maximum F1 across thresholds
    threshold_at_Fmax   — score value at F_max
    Pr_Fmax, Rc_Fmax    — operating-point coords

§3.4 ROC + TN-universe handling is deferred from the beta slice.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .metrics import _safe_div


@dataclass(frozen=True)
class CurveResult:
    prc: pd.DataFrame
    aupr: float
    f_max: float
    threshold_at_fmax: float
    pr_fmax: float
    rc_fmax: float


def pr_curve(
    pred: pd.DataFrame,
    truth: pd.DataFrame,
    rank: str = "species",
    taxgroup: str = "Viruses",
    score_col: str = "readn",
    top_rank: str = "division",
) -> CurveResult:
    """Compute the threshold-swept PR curve at one (taxgroup, rank).

    The pred DataFrame is restricted to rows where pred[top_rank] == taxgroup
    before the sweep. The score column must contain non-negative numbers.
    The truth set is the set of unique `rank` values where truth[top_rank] == taxgroup.
    """
    if score_col not in pred.columns:
        raise ValueError(f"Score column '{score_col}' missing from pred")

    pred_grp = pred[pred[top_rank] == taxgroup] if top_rank in pred.columns else pred.iloc[0:0]
    true_grp = truth[truth[top_rank] == taxgroup] if top_rank in truth.columns else truth.iloc[0:0]
    true_set = set(true_grp[rank].dropna().unique())

    # Aggregate per taxon (max score across rows for that taxon — robust to duplicates).
    if len(pred_grp) == 0:
        agg = pd.DataFrame({rank: [], score_col: []})
    else:
        agg = pred_grp.groupby(rank, dropna=True)[score_col].max().reset_index()

    # Threshold sweep: walk over unique scores in descending order, accumulating predictions.
    sweep = agg.sort_values(score_col, ascending=False, kind="mergesort").reset_index(drop=True)

    rows: list[dict] = []
    pred_so_far: set = set()
    if len(sweep) > 0:
        # Add a leading row at threshold = +inf (no predictions yet) so the curve starts at Rc=0.
        rows.append({"threshold": float("inf"), "TP": 0, "FP": 0, "FN": len(true_set), "Pr": 0.0, "Rc": 0.0, "F1": 0.0})
        # Walk thresholds top-down.
        for i in range(len(sweep)):
            taxon = sweep.iloc[i][rank]
            t = float(sweep.iloc[i][score_col])
            pred_so_far.add(taxon)
            # Emit one row per distinct threshold (skip until score changes).
            if i + 1 < len(sweep) and sweep.iloc[i + 1][score_col] == t:
                continue
            tp = len(pred_so_far & true_set)
            fp = len(pred_so_far - true_set)
            fn = len(true_set - pred_so_far)
            pr = _safe_div(tp, tp + fp)
            rc = _safe_div(tp, tp + fn)
            f1 = 2 * pr * rc / (pr + rc) if (pr + rc) > 0 else 0.0
            rows.append({"threshold": t, "TP": tp, "FP": fp, "FN": fn, "Pr": pr, "Rc": rc, "F1": f1})
    else:
        rows.append({"threshold": float("inf"), "TP": 0, "FP": 0, "FN": len(true_set), "Pr": 0.0, "Rc": 0.0, "F1": 0.0})

    prc_df = pd.DataFrame(rows)

    # AUPR: trapezoidal integration over recall (ascending).
    rc_arr = prc_df["Rc"].to_numpy()
    pr_arr = prc_df["Pr"].to_numpy()
    order = np.argsort(rc_arr, kind="mergesort")
    rc_sorted = rc_arr[order]
    pr_sorted = pr_arr[order]
    aupr = float(np.trapezoid(pr_sorted, rc_sorted)) if len(rc_sorted) >= 2 else 0.0

    # F-max operating point.
    if prc_df["F1"].max() > 0:
        idx = int(prc_df["F1"].idxmax())
        f_max = float(prc_df.loc[idx, "F1"])
        thr_fmax = float(prc_df.loc[idx, "threshold"])
        pr_fmax = float(prc_df.loc[idx, "Pr"])
        rc_fmax = float(prc_df.loc[idx, "Rc"])
    else:
        f_max = 0.0
        thr_fmax = float("nan")
        pr_fmax = 0.0
        rc_fmax = 0.0

    return CurveResult(
        prc=prc_df,
        aupr=aupr,
        f_max=f_max,
        threshold_at_fmax=thr_fmax,
        pr_fmax=pr_fmax,
        rc_fmax=rc_fmax,
    )
