#!/usr/bin/env python3
"""Framework C — rescore the `barcode79` replay at species as well as genus.

Module C scores one rank per replay, and the `barcode79` run was scored at
genus (`--rank genus --score-col readn_em`). The stabilisation figure needs
species-level L1 on the same batches, and the near-synonym collapse the 16S
benchmark applies to its headline numbers (data/truth/*_collapse.tsv)
is not something Module C knows about.

Neither is recoverable from the run directory: `watch` overwrites
`results/abund_table.tsv` after every batch, so only the final one survives.
What does survive is every batch's annotation table, and the cumulative
abundance is a pure function of those — `watch` re-runs EM over the
accumulated annotations after each batch (realtime/processor.py `_run_em`),
warm-starting from the previous F(t). This script repeats exactly that:

    for k in 1..26:  EM(concat(annot1 of batches 1..k), F_init=F(k-1))

`readn_em = prob x reads entering EM`, so a rank's score is its summed
posterior; taxid -> lineage is read off the final abundance table, which
carries every taxid the run ever saw. Scoring then follows
`ngs_benchmark.realtime.runner.score_snapshot`: L1 from `abundance_metrics`,
set metrics over taxa at or above `--read-min` reads.

The strict genus columns reproduce the recorded timeseries.tsv to within 1e-4
L1 (EM stops on a log-likelihood threshold, so a rerun is not bit-identical),
which is what makes the species columns trustworthy; the run asserts it.

Set metrics are carried three ways per batch: Module C's rule at each rank
(`*_strict`), the same with near-synonyms collapsed (`*_coll`), and the 16S
benchmark's row-filter-then-collapse cell with the 1 % tail trim (`*_tail`).

    metrics/framework_C/barcode79.ranks.tsv

Also `*_floor` per rank: collapsed, 10 reads plus the 0.01 % abundance floor of
Framework A and the static 16S tables.

Usage:
    python3 jobs/framework_C/framework_C.barcode79.ranks.py
    python3 jobs/framework_C/framework_C.barcode79.ranks.py --library gutmock_sup_r1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]   # jobs/<framework>/<script>
sys.path.insert(0, str(REPO / "ngs-benchmark"))

from ngs_benchmark.detection.abundance import (  # noqa: E402
    abundance_metrics,
    apply_collapse,
    load_collapse_map,
)
from ngs_benchmark.detection.metrics import set_metrics  # noqa: E402
from ngs_benchmark.realtime.metrics import set_scores  # noqa: E402
from ngs_benchmark.tailtrim import tail_trim  # noqa: E402
from emtrace.analysis.em import em_loop  # noqa: E402

# --library: run dir results/framework_C/<library>, output
# metrics/framework_C/<library>.ranks.tsv. The gut-mock replicates
# (gutmock_sup_r1..r4) share one truth and collapse map.
TRUTHS = {
    "barcode79": ("zymobiomics_16S_truth.tsv", "zymobiomics_16S_collapse.tsv"),
    "gutmock":   ("gutmock_16S_truth.tsv", "gutmock_16S_collapse.tsv"),
}

# emtrace.config.16S.frameC.N25.yaml, `em:` block — the replay's settings.
EM_TEMP = 10.0
EM_MAX_ITER = 1000
EM_LOGL_THRESHOLD = 0.01
EM_EXCLUDE_NO_SPECIES = False

TRUTH_WEIGHT = "pct_16S_only"
ANNOT_COLS = ["qseqid", "staxid", "bitscore", "species"]

# The 1 % tail trim of the 16S benchmark (framework_16S.score.py, axis 4). It is a
# scoring step, not a pipeline step, and it does NOT work like Framework C's cutoff:
# it drops abundance-table ROWS and only then collapses to a rank, where Module C
# sums to the rank first and thresholds after. Reproduced here so the tail-trimmed
# precision in Figure 3 is the same number the 16S tables report.
TAIL = 0.01
READN_MIN = 10.0
# The Framework A / 16S-table floor: a taxon must also hold >= 0.01 % of the rank's
# summed readn_em among taxa that pass READN_MIN (`*_floor` columns).
MIN_REL = 0.0001


def lineage(run_dir: Path, library: str) -> pd.DataFrame:
    """taxid -> species/genus/division, from the final (cumulative) abundance table."""
    abund = pd.read_csv(run_dir / f"pipeline/{library}/results/abund_table.tsv",
                        sep="\t", dtype={"taxid": str})
    return abund[["taxid", "species", "genus", "division"]].drop_duplicates("taxid")


def score_tail(abund: pd.DataFrame, truth: pd.DataFrame) -> dict:
    """Set metrics the 16S benchmark's way: row filter, then collapse to rank.

    `filter_cell(readn_min=10, tail=0.01)` on `readn_em`, then `set_metrics` over
    Bacteria — framework_16S.score.py, the cell the 16S tables report.

    The trim is `ngs_benchmark.tailtrim`, the same code the pipeline's
    `report.tail` and framework_16S.score.py's `filter_cell` both run.
    """
    counts = pd.to_numeric(abund["readn_em"], errors="coerce").fillna(0)
    cell = tail_trim(abund[counts >= READN_MIN], TAIL, "readn_em",
                     preserve_order=False)
    if cell.empty:
        return {}
    m = set_metrics(cell, truth, ranks=("genus", "species"),
                    taxgroups=("Bacteria",)).set_index("rank")
    return {f"{col}_{rank}_tail": m.loc[rank, col]
            for rank in ("genus", "species")
            for col in ("n_pred", "TP", "FP", "FN", "Pr", "Rc", "F1")}


def score(abund: pd.DataFrame, truth: pd.DataFrame, rank: str,
          read_min: float, min_rel: float = 0.0) -> dict:
    """score_snapshot's metrics for one rank (runner.py:238)."""
    ranked = abund.groupby(rank, dropna=True)["readn_em"].sum().sort_values(ascending=False)
    ranked = ranked[ranked.index.astype(str).str.len() > 0]
    m = abundance_metrics(abund, truth, ranks=(rank,), pred_weight="readn_em")
    present = truth[pd.to_numeric(truth[TRUTH_WEIGHT], errors="coerce").fillna(0) > 0]
    truth_taxa = set(present[rank].dropna().astype(str))
    passed = ranked[ranked >= read_min]
    if min_rel > 0 and passed.sum() > 0:
        passed = passed[passed / passed.sum() >= min_rel]
    pred_taxa = set(passed.index.astype(str))
    s = set_scores(pred_taxa, truth_taxa)
    return {
        "top": str(ranked.index[0]),
        "L1": float(m.iloc[0]["L1"]),
        "other_pct": float(m.iloc[0]["other_pct"]),
        "n_pred": s["n_pred"], "n_true": s["n_true"],
        "TP": s["TP"], "FP": s["FP"], "FN": s["FN"],
        "Pr": s["Pr"], "Rc": s["Rc"], "F1": s["F1"],
        "recovered": len(truth_taxa & set(ranked.index.astype(str))),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--library", default="barcode79",
                    help="replay library: barcode79 or gutmock_<model>_r<N>")
    ap.add_argument("--run", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--read-min", type=float, default=10.0)
    ap.add_argument("--no-check", action="store_true",
                    help="skip the strict-genus agreement check against timeseries.tsv")
    args = ap.parse_args()
    lib = args.library
    args.run = args.run or REPO / "results/framework_C" / lib
    args.out = args.out or REPO / "metrics/framework_C" / f"{lib}.ranks.tsv"
    truth_file, collapse_file = TRUTHS[lib.split("_", 1)[0]]

    truth = pd.read_csv(REPO / "data/truth" / truth_file, sep="\t")
    cmap = load_collapse_map(REPO / "data/truth" / collapse_file)
    tax = lineage(args.run, lib)
    ts = pd.read_csv(args.run / "timeseries.tsv", sep="\t")

    batches = sorted((args.run / f"pipeline/{lib}/batches").glob("batch_*"))
    accum: list[pd.DataFrame] = []
    F: dict[str, float] | None = None
    rows = []
    for k, bdir in enumerate(batches, start=1):
        ann = pd.read_csv(bdir / "annotation" / "annot1.tsv", sep="\t",
                          usecols=ANNOT_COLS, dtype={"staxid": str})
        accum.append(ann)
        full = pd.concat(accum, ignore_index=True)
        res = em_loop(full, temp=EM_TEMP, max_iter=EM_MAX_ITER,
                      log_L_threshold=EM_LOGL_THRESHOLD,
                      exclude_no_species=EM_EXCLUDE_NO_SPECIES, F_init=F)
        F = res.F
        abund = pd.DataFrame({"taxid": list(F), "prob": list(F.values())})
        abund["readn_em"] = abund["prob"] * res.n_queries
        abund = abund.merge(tax, on="taxid", how="left")

        row = {"batch": k, "reads_em": res.n_queries}
        collapsed = apply_collapse(abund, cmap)
        for label, table in (("strict", abund), ("coll", collapsed)):
            for rank in ("genus", "species"):
                row.update({f"{key}_{rank}_{label}": value
                            for key, value in score(table, truth, rank, args.read_min).items()})
        for rank in ("genus", "species"):   # collapsed, 10 reads + 0.01 % floor
            row.update({f"{key}_{rank}_floor": value
                        for key, value in score(collapsed, truth, rank, args.read_min,
                                                MIN_REL).items()})
        # Collapsed only: the tail-trimmed cell is the 16S benchmark's headline basis.
        row.update(score_tail(collapsed, truth))
        rows.append(row)
        print(f"[batch {k:>2}] EM {res.n_iterations} it, {res.n_queries} reads  "
              f"genus L1={row['L1_genus_strict']:.4f}  "
              f"species L1={row['L1_species_coll']:.4f} (collapsed)  "
              f"species Pr={row['Pr_species_coll']:.3f} / "
              f"{row.get('Pr_species_tail', float('nan')):.3f} (1 % tail)", flush=True)

    out = pd.DataFrame(rows)
    out.insert(0, "library", lib)
    out = out.merge(ts[["batch", "seq_time_s", "wall_s", "reads_cum"]], on="batch", how="left")

    if not args.no_check:
        delta = (out["L1_genus_strict"].to_numpy() - ts["L1"].to_numpy())
        worst = abs(delta).max()
        print(f"[check] strict genus L1 vs timeseries.tsv: max |delta| = {worst:.2e}")
        if worst > 1e-3:
            raise SystemExit("Reconstructed EM does not reproduce the recorded genus L1.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    out = args.out.resolve()
    shown = out.relative_to(REPO) if out.is_relative_to(REPO) else out
    print(f"Table -> {shown}")


if __name__ == "__main__":
    main()
