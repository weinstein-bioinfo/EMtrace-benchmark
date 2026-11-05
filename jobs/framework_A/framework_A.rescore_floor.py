#!/usr/bin/env python3
"""Framework A — rescore detection with the 16S benchmark's abundance floor.

Three scoring cells, every pipeline, every Framework A dataset:

  base      read_min 10 on readn, no floor            (earlier scoring, for comparison)
  em        read_min 10 on readn_em, no floor         (count switch only)
  floor_em  read_min 10 on readn_em, then a 0.01 %    (the scoring reported in the
            relative-abundance floor on readn_em       manuscript; Emu's published
                                                       floor, as in framework_16S.score.py)

readn_em falls back to readn for tables without it (EPI2ME, CZ ID, MEGAN-LR), so only
EMtrace's count basis changes between base and em.

Reads a tool places above a rank (MEGAN-LR's LCA, exported as "unknown") are no
taxon call at that rank (--drop-unknown, 9 Oct 2026); they still count in the
floor's denominator.

The floor's denominator leaves out division Primates (--min-rel-exclude).
EPI2ME reports human reads as taxa (99.99 % of a Buddle library) where
EMtrace and CZ ID remove them first; with them in the total, 0.01 % of all
reads exceeds EPI2ME's entire viral signal on Buddle.

Inputs (all local):
  EMtrace  Buddle: results/framework_A/emtrace/<ds>/<acc>/  (Roihu runs)
              Petersen, Gauthier: Framework B v20260924_refauto runs — the same
              config, host filter and --ann1 minimap.refseq.abv.1g as Framework A,
              newer code; Framework A has no local EMtrace runs for them.
  EPI2ME      results/framework_A/epi2me/<ds>/<acc>/<acc>_epi2me.bracken.tsv
  CZ ID       results/framework_A/czid/<ds>/<acc>/<acc>*_taxon_report.csv
  MEGAN-LR    results/framework_A/megan_lr/<ds>/<acc>/<acc>_megan_taxcount.tsv
              (taxid, read count: MEGAN6 long-read mode, read-count export, 9 Oct 2026)

Outputs:
  metrics/framework_A_floor/<cell>/<ds>/   detection-batch tree per cell
  metrics/framework_A_floor/summary.tsv    micro-averaged species Pr/Rc/F1

Run from the repo root:
  python3 jobs/framework_A/framework_A.rescore_floor.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[2]
CONFIG = REPO / "jobs/config.yaml"
RES_A = REPO / "results/framework_A"
RES_B = REPO / "results/framework_b"
OUT = REPO / "metrics/framework_A_floor"

READ_MIN = 10
FLOOR = 0.0001   # 0.01 %
CELLS = {
    "base":     {"filter_col": "readn",    "min_rel": 0.0},
    "em":       {"filter_col": "readn_em", "min_rel": 0.0},
    "floor_em": {"filter_col": "readn_em", "min_rel": FLOOR},
}
# Framework B directory names spell the dataset "petersen".
EMTRACE_B_DIR = {
    "peterson_2026": "sra_petersen_2026_v20260924_refauto",
    "gauthier_2021": "sra_gauthier_2021_v20260924_refauto",
}
DATASETS = ("buddle_2024_dna", "buddle_2024_rna", "peterson_2026", "gauthier_2021")


def emtrace_path(ds: str, acc: str) -> Path:
    root = RES_B / EMTRACE_B_DIR[ds] if ds in EMTRACE_B_DIR else RES_A / "emtrace" / ds
    return root / acc / "results/reports/abund_table.tsv"


def czid_path(ds: str, acc: str) -> Path | None:
    hits = sorted((RES_A / "czid" / ds / acc).glob(f"{acc}*_taxon_report.csv"))
    return hits[0] if len(hits) == 1 else None


def samples(ds: str, cfg: dict) -> list[dict]:
    """Config samples; for sample-keyed truths, only the samples the truth scores
    (Gauthier's 16 positives, not the negative controls)."""
    spec = cfg["datasets"][ds]
    rows = [s for s in spec["samples"] if s.get("sra") and s["sra"] != "PLACEHOLDER"]
    truth = pd.read_csv(REPO / spec["truth_file"], sep="\t")
    if "sample_id" in truth.columns:
        keep = set(truth["sample_id"].astype(str))
        rows = [s for s in rows if s["sra"] in keep]
    return rows


def manifest(ds: str, cfg: dict) -> pd.DataFrame:
    rows = []
    for s in samples(ds, cfg):
        acc = s["sra"]
        preds = {
            "emtrace": emtrace_path(ds, acc),
            "epi2me": RES_A / "epi2me" / ds / acc / f"{acc}_epi2me.bracken.tsv",
            "czid": czid_path(ds, acc),
            "megan_lr": RES_A / "megan_lr" / ds / acc / f"{acc}_megan_taxcount.tsv",
        }
        for pipeline, path in preds.items():
            if path is None or not path.is_file():
                sys.exit(f"[FATAL] {pipeline} {ds}/{acc}: prediction missing ({path})")
            rows.append({
                "pipeline": pipeline, "sample_id": acc,
                "path": str(path.relative_to(REPO)), "fmt": pipeline,
                "gc_per_ml": s.get("gc_per_ml", ""),
            })
    return pd.DataFrame(rows)


def score(cell: str, ds: str, cfg: dict) -> Path:
    outdir = OUT / cell / ds
    outdir.mkdir(parents=True, exist_ok=True)
    man = manifest(ds, cfg)
    man_path = outdir / "manifest.tsv"
    man.to_csv(man_path, sep="\t", index=False)
    nb = cfg["ngs_benchmark"]
    cmd = [
        "ngs-benchmark", "detection-batch",
        "--pred-list", str(man_path),
        "--truth", str(REPO / cfg["datasets"][ds]["truth_file"]),
        "--outdir", str(outdir),
        "--ranks", nb["ranks"],
        "--taxgroup", nb["taxgroup"],
        "--score-col", nb["score_col"],
        "--read-min", str(READ_MIN),
        "--min-rel", str(CELLS[cell]["min_rel"]),
        "--filter-col", CELLS[cell]["filter_col"],
        "--min-rel-exclude", "Primates",
        "--tail", "0",
        "--drop-unknown",
    ]
    if ds.startswith("buddle"):
        cmd.append("--dilution-series")
    print(f"[score] {cell:8s} {ds}  ({len(man)} predictions)")
    subprocess.run(cmd, cwd=REPO, check=True, stdout=subprocess.DEVNULL)
    return outdir / "combined_metrics.tsv"


def main() -> None:
    cfg = yaml.safe_load(CONFIG.read_text())
    rows = []
    for cell in CELLS:
        for ds in DATASETS:
            m = pd.read_csv(score(cell, ds, cfg), sep="\t")
            m = m[(m["rank"] == "species") & (m["taxgroup"] == "Viruses")]
            for pipeline, g in m.groupby("pipeline"):
                tp, fp, fn = g["TP"].sum(), g["FP"].sum(), g["FN"].sum()
                pr = tp / (tp + fp) if tp + fp else float("nan")
                rc = tp / (tp + fn) if tp + fn else float("nan")
                f1 = 2 * pr * rc / (pr + rc) if pr + rc else float("nan")
                rows.append({"cell": cell, "dataset": ds, "pipeline": pipeline,
                             "samples": len(g), "TP": tp, "FP": fp, "FN": fn,
                             "Pr": round(pr, 3), "Rc": round(rc, 3), "F1": round(f1, 3)})
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "summary.tsv", sep="\t", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
