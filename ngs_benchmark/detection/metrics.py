"""§3.2 Set metrics.

For each (taxgroup, rank) combination, compute TP, FP, FN from the set of
unique taxon labels in pred vs truth. Direct port of NGSlib.R:212-219.
"""

from __future__ import annotations

import warnings
import pandas as pd

from ngs_benchmark.schema import TAXRANKS


def _safe_div(num: float, den: float) -> float:
    """`num / den` with 0/0 := 0 — matches the R reference convention."""
    return num / den if den > 0 else 0.0


def set_metrics(
    pred: pd.DataFrame,
    truth: pd.DataFrame,
    ranks: tuple[str, ...] = ("species", "genus", "family"),
    taxgroups: tuple[str, ...] = ("Viruses", "Bacteria"),
    top_rank: str = "division",
) -> pd.DataFrame:
    """Compute set metrics per (taxgroup, rank).

    Both inputs must be the §3.1-normalised DataFrames. The result has one row
    per (taxgroup, rank) with columns:
        taxgroup, rank, n_pred, n_true, TP, FP, FN, Pr, Rc, F1, mean_accuracy
    """
    for rank in ranks:
        if rank not in TAXRANKS:
            raise ValueError(f"Unknown rank '{rank}'. Valid: {TAXRANKS}")
        if rank not in pred.columns or rank not in truth.columns:
            raise ValueError(f"Rank column '{rank}' missing in pred or truth")

    rows: list[dict] = []
    for group in taxgroups:
        pred_grp = pred[pred[top_rank] == group] if top_rank in pred.columns else pred.iloc[0:0]
        true_grp = truth[truth[top_rank] == group] if top_rank in truth.columns else truth.iloc[0:0]
        for rank in ranks:
            pred_set = set(pred_grp[rank].dropna().unique())
            true_set = set(true_grp[rank].dropna().unique())
            tp_set = pred_set & true_set
            fp_set = pred_set - true_set
            fn_set = true_set - pred_set
            tp, fp, fn = len(tp_set), len(fp_set), len(fn_set)
            pr = _safe_div(tp, tp + fp)
            rc = _safe_div(tp, tp + fn)
            if pr + rc == 0:
                warnings.warn(
                    f"Both Pr and Rc are zero for taxgroup={group}, rank={rank}; "
                    f"setting F1 = 0.",
                    stacklevel=2,
                )
                f1 = 0.0
            else:
                f1 = 2 * pr * rc / (pr + rc)
            mean_acc = (pr + rc) / 2
            rows.append(
                {
                    "taxgroup": group,
                    "rank": rank,
                    "n_pred": len(pred_set),
                    "n_true": len(true_set),
                    "TP": tp,
                    "FP": fp,
                    "FN": fn,
                    "Pr": pr,
                    "Rc": rc,
                    "F1": f1,
                    "mean_accuracy": mean_acc,
                }
            )
    return pd.DataFrame(rows)


def per_taxon_table(
    pred: pd.DataFrame,
    truth: pd.DataFrame,
    ranks: tuple[str, ...] = ("species", "genus", "family"),
    taxgroups: tuple[str, ...] = ("Viruses", "Bacteria"),
    top_rank: str = "division",
) -> pd.DataFrame:
    """One row per (taxgroup, rank, taxon) with status ∈ {TP, FP, FN}.

    Columns: taxgroup, rank, taxon, status, readn_pred, readn_true.
    """
    rows: list[dict] = []
    for group in taxgroups:
        pred_grp = pred[pred[top_rank] == group] if top_rank in pred.columns else pred.iloc[0:0]
        true_grp = truth[truth[top_rank] == group] if top_rank in truth.columns else truth.iloc[0:0]
        for rank in ranks:
            # Aggregate readn per taxon at the requested rank.
            pred_readn = (
                pred_grp.groupby(rank, dropna=True)["readn"].sum() if len(pred_grp) else pd.Series(dtype=float)
            )
            true_readn = (
                true_grp.groupby(rank, dropna=False)["readn"].sum() if len(true_grp) else pd.Series(dtype=float)
            )
            pred_set = set(pred_readn.index)
            true_set = set(true_readn.index)
            for taxon in sorted(pred_set | true_set, key=str):
                in_pred = taxon in pred_set
                in_true = taxon in true_set
                if in_pred and in_true:
                    status = "TP"
                elif in_pred:
                    status = "FP"
                else:
                    status = "FN"
                rows.append(
                    {
                        "taxgroup": group,
                        "rank": rank,
                        "taxon": taxon,
                        "status": status,
                        "readn_pred": float(pred_readn.get(taxon, 0.0)),
                        "readn_true": float(true_readn.get(taxon, 0.0)),
                    }
                )
    return pd.DataFrame(rows)
