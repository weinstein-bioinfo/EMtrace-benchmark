#!/usr/bin/env python3
"""Framework C — collect Module C replay results into the metrics tables.

`ngs-benchmark realtime` already scores every batch of a replay and writes
timeseries / first_call / latency / accuracy.tsv beside the run, and, for a
panel-scored library (the tick pools), panel.tsv and panel_summary.tsv. This
script gathers the manuscript's replays (RUNS) under results/framework_C/
into one place:

    metrics/framework_C/summary.tsv     one row per replay (manuscript Table 2)
    metrics/framework_C/timeseries.tsv  per batch, all replays
    metrics/framework_C/accuracy.tsv    first_call / stabilised / final, all replays
    metrics/framework_C/panel.tsv       per panel taxon, panel-scored replays

The `run` column is the directory name (e.g. `ticks_58-37.refseq_abv1g_k11w5`),
so one library replayed against two databases stays two rows. Target-read
curves per batch are built by framework_C.curve.py, not here.

For a panel-scored replay the first call is replaced by the panel's first
detection: the batch from which every PCR-positive taxon has stayed at or above
the read threshold (`panel_first_batch`). These libraries have no read start
times, so their sequencing-time columns and throughput ratio are empty.

Usage:
    python3 jobs/framework_C/framework_C.score.py            # RUNS
    python3 jobs/framework_C/framework_C.score.py RUN [RUN ...]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]   # jobs/<framework>/<script>
RESULTS = REPO / "results/framework_C"
METRICS = REPO / "metrics/framework_C"

# The replays reported in the manuscript (Table 2, Fig. 4).
RUNS = (
    "Orbivirus_kevo.refseq_abv1g_k11w5",
    "Orbivirus_kevo.refseq_abv1g_k11w5_dmnd",
    "TBEV_FBC92293.refseq_abv1g",
    "ticks_58-37.refseq_abv1g_k11w5",
    "ticks_58-37.refseq_abv1g_k11w5_dmnd",
    "ticks_35-33.refseq_abv1g_k11w5",
    "ticks_35-33.refseq_abv1g_k11w5_dmnd",
    "barcode79",
    "gutmock_sup_r1",
    "gutmock_sup_r2",
    "gutmock_sup_r3",
    "gutmock_sup_r4",
)
LATENCY_TARGET_S_PER_10K = 120.0  # the project notes: < 2 min per 10 k reads


def read_provenance(path: Path) -> dict:
    prov = {}
    for line in path.read_text().splitlines():
        key, _, value = line.partition(":")
        prov[key.strip()] = value.strip()
    return prov


def truth_taxon(prov: dict) -> str | None:
    m = re.search(r"--truth-taxon (.+?) --", prov.get("scoring", ""))
    return m.group(1) if m else None


def scored_rank(prov: dict) -> str | None:
    m = re.search(r"--rank (\S+)", prov.get("scoring", ""))
    return m.group(1) if m else None


def panel_columns(run_dir: Path) -> tuple[dict, pd.DataFrame]:
    """Panel counts and first detection of all PCR targets; empty if not panel-scored."""
    if not (run_dir / "panel_summary.tsv").is_file():
        return {}, pd.DataFrame()
    ps = pd.read_csv(run_dir / "panel_summary.tsv", sep="\t").iloc[0]
    panel = pd.read_csv(run_dir / "panel.tsv", sep="\t")
    pos = panel[panel["evidence"] == "pcr_pos"]
    # Every PCR target must be detected; the panel is detected when the last one is.
    first = (int(pos["first_detection_batch"].max())
             if len(pos) and (pos["outcome"] == "detected").all() else pd.NA)
    cols: dict[str, object] = {k: int(ps[k]) for k in ("pcr_pos_detected", "pcr_pos_n", "pcr_neg_clear",
                                     "pcr_neg_n", "reported_detected", "reported_n")}
    cols["panel_first_batch"] = first
    return cols, panel


def score_run(run_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    prov = read_provenance(run_dir / "provenance.txt")
    library = prov["library"]
    run = run_dir.name
    ts = pd.read_csv(run_dir / "timeseries.tsv", sep="\t")
    acc = pd.read_csv(run_dir / "accuracy.tsv", sep="\t")
    # first_call.tsv exists only for single-taxon (--truth-taxon) replays that
    # reached a call; a profiling or --reads-only replay has none.
    fc_path = run_dir / "first_call.tsv"
    fc = pd.read_csv(fc_path, sep="\t") if fc_path.is_file() else pd.DataFrame()
    lat = pd.read_csv(run_dir / "latency.tsv", sep="\t").iloc[0]
    panel_cols, panel = panel_columns(run_dir)

    meta = {"run": run, "anns": prov.get("anns"), "threads": prov.get("threads"),
            "host": f"{prov.get('cpu')} / {prov.get('ram_gb')} GB"}
    taxon = truth_taxon(prov)
    ts.insert(0, "run", run)
    acc.insert(0, "run", run)
    for key in ("anns",):
        ts.insert(2, key, meta[key])
        acc.insert(2, key, meta[key])
    if len(panel):
        panel.insert(0, "run", run)

    last = ts.iloc[-1]
    stab = acc[acc["event"] == "stabilised"]
    batch_reads = ts["reads_cum"].diff().fillna(ts["reads_cum"])
    seq_s = last["seq_time_s"]
    summary = {
        **meta,
        "library": library,
        "truth": taxon or prov.get("scoring"),
        "rank": scored_rank(prov),
        "batches": int(lat["batches"]),
        "reads": int(last["reads_cum"]),
        "batch_reads_median": int(batch_reads.median()),
        "seq_span_h": round(seq_s / 3600, 2) if pd.notna(seq_s) else pd.NA,
        "first_call_batch": fc["batch"].iloc[0] if len(fc) else pd.NA,
        "first_call_seq_min": round(fc["seq_time_s"].iloc[0] / 60, 1) if len(fc) else pd.NA,
        "first_call_wall_min": round(fc["wall_s"].iloc[0] / 60, 1) if len(fc) else pd.NA,
        "first_call_reads": fc["reads"].iloc[0] if len(fc) else pd.NA,
        "first_call_em_prob": fc["em_prob"].iloc[0] if len(fc) else pd.NA,
        **panel_cols,
        "latency_median_s": lat["median_s"],
        "latency_p95_s": lat["p95_s"],
        "s_per_10k_reads": round(lat["s_per_1k_reads"] * 10, 1),
        "meets_2min_per_10k": lat["s_per_1k_reads"] * 10 < LATENCY_TARGET_S_PER_10K,
        "throughput_ratio": lat["ratio"],
        "stabilised_batch": stab["batch"].iloc[0] if len(stab) else pd.NA,
        "stabilised_reads": stab["reads_cum"].iloc[0] if len(stab) else pd.NA,
        "final_top_taxon": last["top_taxon"],
        "final_em_prob": last["em_prob"],
        "final_TP": last["TP"], "final_FP": last["FP"], "final_FN": last["FN"],
        "final_F1": last["F1"],
        "final_L1": last["L1"],
        "wall_total_min": round(int(prov.get("wall_s", 0)) / 60, 1),
    }
    return ts, acc, panel, summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("runs", nargs="*", default=list(RUNS),
                    help="replay directories under --results (default: the manuscript's)")
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--outdir", type=Path, default=METRICS)
    args = ap.parse_args()

    runs = [args.results / r for r in args.runs]
    missing = [r.name for r in runs
               if not ((r / "timeseries.tsv").is_file() and (r / "provenance.txt").is_file())]
    if missing:
        raise SystemExit(f"Not a completed replay under {args.results}: {', '.join(missing)}")

    ts_all, acc_all, panel_all, summaries = [], [], [], []
    for run_dir in runs:
        ts, acc, panel, summary = score_run(run_dir)
        ts_all.append(ts)
        acc_all.append(acc)
        panel_all.append(panel)
        summaries.append(summary)
        print(f"[scored] {run_dir.name}: {summary['batches']} batches, "
              f"final top {summary['final_top_taxon']}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(args.outdir / "summary.tsv", sep="\t", index=False)
    pd.concat(ts_all).to_csv(args.outdir / "timeseries.tsv", sep="\t", index=False)
    pd.concat(acc_all).to_csv(args.outdir / "accuracy.tsv", sep="\t", index=False)
    pd.concat(panel_all).to_csv(args.outdir / "panel.tsv", sep="\t", index=False)
    out = args.outdir.resolve()
    shown = out.relative_to(REPO) if out.is_relative_to(REPO) else out
    print(f"Tables → {shown}/{{summary,timeseries,accuracy,panel}}.tsv")


if __name__ == "__main__":
    main()
