"""§3.6 Dilution-series detection metrics.

Per `(pipeline, gc_per_ml)` level, compute:
    sensitivity = recall      (primary; manuscript Figure 2)
    TP, FP, FN, precision     (augmentary)
    fp_unique_taxa            (per-level deduplicated FP count)
    fp_persistent_taxa        (FPs called at every level of the same matrix)
    low_n                     (n_samples < 3)

Inputs are a list of per-sample prediction frames (already §3.1-normalised
and post-Buddle-filter if applicable) and a truth table containing
`sample_id`, `gc_per_ml`, and per-virus rows. An optional `matrix_id` column
in the truth table groups dilution series of the same biological matrix.
"""

from __future__ import annotations

from dataclasses import dataclass
import pandas as pd

from .metrics import _safe_div


@dataclass(frozen=True)
class SamplePred:
    """One sample's prediction frame, pre-keyed with metadata."""

    sample_id: str
    gc_per_ml: float
    matrix_id: str | None
    pred: pd.DataFrame  # §3.1-normalised


def dilution_metrics(
    samples: list[SamplePred],
    truth: pd.DataFrame,
    rank: str = "species",
    taxgroup: str = "Viruses",
    top_rank: str = "division",
) -> tuple[pd.DataFrame, dict]:
    """Compute the §3.6 long-format `sensitivity_vs_gcml.tsv`.

    Returns the per-level DataFrame and a `summary` dict containing
    `total_unique_fps` (deduplicated FP-taxon count across all samples).
    """
    truth_grp = truth[truth[top_rank] == taxgroup] if top_rank in truth.columns else truth.iloc[0:0]
    truth_set_global = set(truth_grp[rank].dropna().unique())

    # Bucket samples by gc_per_ml level.
    by_level: dict[float, list[SamplePred]] = {}
    for s in samples:
        by_level.setdefault(s.gc_per_ml, []).append(s)

    # First pass: per-(level, sample) FP sets, for persistence analysis.
    fps_by_sample: dict[tuple[float, str], set] = {}
    for s in samples:
        pred_grp = s.pred[s.pred[top_rank] == taxgroup] if top_rank in s.pred.columns else s.pred.iloc[0:0]
        pred_set = set(pred_grp[rank].dropna().unique())
        fps_by_sample[(s.gc_per_ml, s.sample_id)] = pred_set - truth_set_global

    # Persistence: an FP is "persistent" iff it appears in every level of its matrix.
    has_matrix = any(s.matrix_id is not None for s in samples)
    persistent_by_matrix: dict[str, set] = {}
    if has_matrix:
        # Group samples by matrix → for each matrix collect FP sets per level, intersect them.
        from collections import defaultdict

        by_matrix_level: dict[str, dict[float, set]] = defaultdict(lambda: defaultdict(set))
        for s in samples:
            if s.matrix_id is None:
                continue
            by_matrix_level[s.matrix_id][s.gc_per_ml] |= fps_by_sample[(s.gc_per_ml, s.sample_id)]
        for mid, level_sets in by_matrix_level.items():
            sets = list(level_sets.values())
            persistent_by_matrix[mid] = set.intersection(*sets) if sets else set()

    rows: list[dict] = []
    all_fps: set = set()
    for level, items in sorted(by_level.items()):
        # Pool predictions across samples at this level.
        pred_set: set = set()
        n_tp = 0
        for s in items:
            pred_grp = s.pred[s.pred[top_rank] == taxgroup] if top_rank in s.pred.columns else s.pred.iloc[0:0]
            pred_set |= set(pred_grp[rank].dropna().unique())

        # TP, FP, FN are evaluated against the same global truth set
        # (every dilution shares the same 6 viruses in Buddle's design).
        tp_set = pred_set & truth_set_global
        fp_set = pred_set - truth_set_global
        fn_set = truth_set_global - pred_set
        all_fps |= fp_set

        fp_persistent_taxa: int | None
        if has_matrix:
            matrices_at_level = {s.matrix_id for s in items if s.matrix_id is not None}
            persistent = set()
            for mid in matrices_at_level:
                persistent |= persistent_by_matrix.get(mid, set())
            fp_persistent_taxa = len(fp_set & persistent)
        else:
            fp_persistent_taxa = None  # surfaced as NA in TSV

        tp = len(tp_set)
        fp = len(fp_set)
        fn = len(fn_set)
        sensitivity = _safe_div(tp, tp + fn)
        precision = _safe_div(tp, tp + fp)

        rows.append(
            {
                "gc_per_ml": level,
                "n_samples": len(items),
                "sensitivity": sensitivity,
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "precision": precision,
                "fp_unique_taxa": fp,
                "fp_persistent_taxa": fp_persistent_taxa,
                "low_n": len(items) < 3,
            }
        )

    df = pd.DataFrame(rows)
    summary = {"total_unique_fps": len(all_fps)}
    if not has_matrix and len(samples) > 0:
        summary["fp_persistent_taxa_warning"] = (
            "Truth file has no `matrix_id` column — fp_persistent_taxa is NA. "
            "Add `matrix_id` to enable persistence analysis."
        )
    return df, summary
