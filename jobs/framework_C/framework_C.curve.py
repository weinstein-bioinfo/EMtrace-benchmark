#!/usr/bin/env python3
"""Framework C — per-batch truth-taxon curve of the viral replays, for Fig. 5a.

Module C's `em_prob` is the top *species'* probability. For Orbivirus_kevo,
scored at genus, that is 0.42-0.44 (the reads split between relatives), and the
replays keep no per-batch abundance snapshots. This script rebuilds the curve
of the truth taxon (genus Orbivirus, or species TBEV) from each replay's
per-batch annotations:

  - the warm-started EM is replayed batch by batch on the accumulated
    annot1.tsv tables, with the EM settings of the run's own config, exactly as
    `emtrace watch` runs it (it reproduces Module C's top-species values);
    `truth_em` is the EM probability summed over the truth taxa;
  - `truth_reads_cum` counts reads with any hit to a truth taxon, once however
    many hits they have (all-hits annotation), from any search in the chain
    (minimap2 and DIAMOND); `truth_reads_dmnd_cum` is the DIAMOND part.

Truth taxa are those with <rank> == <taxon> in the run's final cumulative
abundance table. Reads released and sequencing time come from Module C's
timeseries.tsv beside the run.

    metrics/framework_C/<library>.curve.tsv   per run and batch
        libraries: orbi_kevo (genus Orbivirus), tbev (species TBEV),
                   ticks (genus Orthonairovirus; PRJNA878809, no seq_time_s)

`--reads` instead writes the reads curve of the manuscript's Figure (a)
(the project notes): per run, batch and taxon, the `readn` of the
cumulative abundance table that Module C kept for the batch
(`<run>/abund/batch_NNNN.tsv`), summed over the taxon at its rank. This is the
count the panel scoring uses, not `truth_reads_cum` above (reads with any hit).
Taxa: the main target of each library, plus every panel taxon of the tick
pools (jobs/framework_C/panels/).

    metrics/framework_C/realtime_reads.curve.tsv

Must run with the Python that has emtrace installed (it imports its EM):

    PY=$(head -1 "$(command -v emtrace)" | sed -E 's/^#! *//')
    $PY jobs/framework_C/framework_C.curve.py --library tbev
    $PY jobs/framework_C/framework_C.curve.py --reads
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml
from emtrace.analysis.em import em_loop
from emtrace.models.annotation import read_annotation_tsv

REPO = Path(__file__).resolve().parents[2]   # jobs/<framework>/<script>
RESULTS = REPO / "results/framework_C"
METRICS = REPO / "metrics/framework_C"
# library: (replays, rank column of abund_table.tsv, truth taxon)
LIBRARIES = {
    "orbi_kevo": (["Orbivirus_kevo.refseq_abv1g_k11w5", "Orbivirus_kevo.refseq_abv1g_k11w5_dmnd"],
                  "genus", "Orbivirus"),
    "tbev": (["TBEV_FBC92293.refseq_abv1g"], "species", "Orthoflavivirus encephalitidis"),
    "ticks": (["ticks_80-43.refseq_abv1g_k11w5", "ticks_35-33.refseq_abv1g_k11w5_dmnd",
               "ticks_58-37.refseq_abv1g_k11w5_dmnd"], "genus", "Orthonairovirus"),
}


PANELS = REPO / "jobs/framework_C/panels"
# Figure (a) runs: run -> (rank, main target, panel file or None)
READS_RUNS = {
    "Orbivirus_kevo.refseq_abv1g_k11w5": ("genus", "Orbivirus", None),
    "Orbivirus_kevo.refseq_abv1g_k11w5_dmnd": ("genus", "Orbivirus", None),
    "TBEV_FBC92293.refseq_abv1g": ("species", "Orthoflavivirus encephalitidis", None),
    **{f"ticks_{pool}.refseq_abv1g_k11w5{chain}": ("genus", "Orthonairovirus",
                                                   PANELS / f"ticks_{pool}.panel.tsv")
       for pool in ("35-33", "58-37") for chain in ("", "_dmnd")},
}


def reads_curve(run_dir: Path, rank: str, main: str, panel: Path | None) -> pd.DataFrame:
    """Per batch and taxon: readn of the batch's cumulative abundance table."""
    name = run_dir.name
    if panel is None:
        taxa = pd.DataFrame({"taxon": [main], "evidence": ["target"], "source": [""]})
    else:
        taxa = pd.read_csv(panel, sep="\t", comment="#")[["taxon", "evidence", "source"]]
    ts = pd.read_csv(run_dir / "timeseries.tsv", sep="\t")[["batch", "seq_time_s", "reads_cum"]]
    rows = []
    for batch in ts["batch"]:
        path = run_dir / "abund" / f"batch_{batch:04d}.tsv"
        ab = pd.read_csv(path, sep="\t") if path.is_file() else pd.DataFrame(columns=[rank, "readn"])
        readn = ab.groupby(rank)["readn"].sum()
        for t in taxa.itertuples():
            rows.append({"run": name, "library": name.split(".")[0],
                         "chain": "+ DIAMOND" if name.endswith("_dmnd") else "RefSeq only",
                         "batch": batch, "taxon": t.taxon, "rank": rank,
                         "evidence": t.evidence, "source": t.source,
                         "main": t.taxon == main, "readn": int(readn.get(t.taxon, 0))})
    return pd.DataFrame(rows).merge(ts, on="batch", how="left")


def write_reads_curve(results: Path, outdir: Path) -> None:
    curves = []
    for name, (rank, main, panel) in READS_RUNS.items():
        run_dir = results / name
        if not (run_dir / "abund").is_dir():
            raise SystemExit(f"{name}: no abund/ (run framework_C.abund_backfill.py)")
        c = reads_curve(run_dir, rank, main, panel)
        last = c[c["batch"] == c["batch"].max()]
        print(f"[reads] {name}: " + ", ".join(f"{r.taxon} {r.readn:,}" for r in last.itertuples()))
        curves.append(c)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "realtime_reads.curve.tsv"
    pd.concat(curves).to_csv(out, sep="\t", index=False)
    print(f"Table → {out.resolve().relative_to(REPO)}")


def read_provenance(path: Path) -> dict:
    prov = {}
    for line in path.read_text().splitlines():
        key, _, val = line.partition(":")
        prov[key.strip()] = val.strip()
    return prov


def run_curve(run_dir: Path, rank: str, taxon: str) -> pd.DataFrame:
    """Per batch: truth-taxon EM probability and cumulative truth-taxon reads."""
    library = run_dir.name
    pdir = run_dir / "pipeline" / library
    em_cfg = yaml.safe_load((REPO / read_provenance(run_dir / "provenance.txt")["config"])
                            .read_text())["em"]
    ab = pd.read_csv(pdir / "results" / "abund_table.tsv", sep="\t", dtype=str)
    taxa = set(ab.loc[ab[rank] == taxon, "taxid"])
    ts = pd.read_csv(run_dir / "timeseries.tsv", sep="\t")[["batch", "seq_time_s", "reads_cum"]]

    rows, acc, F = [], None, None
    seen, seen_dmnd = set(), set()
    for i, bdir in enumerate(sorted((pdir / "batches").glob("batch_*")), start=1):
        new = read_annotation_tsv(bdir / "annotation" / "annot1.tsv")
        acc = new if acc is None else pd.concat([acc, new], ignore_index=True)
        res = em_loop(acc, temp=em_cfg["temp"], max_iter=em_cfg["max_iter"],
                      log_L_threshold=em_cfg["log_L_threshold"],
                      exclude_no_species=em_cfg["exclude_no_species"], F_init=F)
        F = dict(res.F)
        # The EM leaves NaN for a few taxa whose share collapses to zero (minor
        # taxa, <= 34 reads); the abundance table reports them as 0, so do we.
        probs = [p for t, p in F.items() if str(t) in taxa and p == p]
        hits = new[new["staxid"].isin(taxa)]
        seen.update(hits["qseqid"])
        seen_dmnd.update(hits.loc[hits["search"].str.startswith("Diamond"), "qseqid"])
        rows.append({"run": library, "batch": i,
                     "truth_em": round(sum(probs), 6),
                     "top_species_em": round(max(probs, default=0.0), 6),
                     "truth_reads_cum": len(seen),
                     "truth_reads_dmnd_cum": len(seen_dmnd)})
    out = pd.DataFrame(rows).merge(ts, on="batch", how="left")
    out["truth_share"] = (out["truth_reads_cum"] / out["reads_cum"]).round(5)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--library", choices=LIBRARIES, default="orbi_kevo")
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--outdir", type=Path, default=METRICS)
    ap.add_argument("--runs", nargs="+", help="default: the library's replays")
    ap.add_argument("--reads", action="store_true",
                    help="write realtime_reads.curve.tsv (Figure a) instead of a library curve")
    args = ap.parse_args()
    if args.reads:
        write_reads_curve(args.results, args.outdir)
        return
    runs, rank, taxon = LIBRARIES[args.library]

    curves = []
    for name in args.runs or runs:
        run_dir = args.results / name
        if not (run_dir / "timeseries.tsv").is_file():
            print(f"[skip] {name}: no completed replay (timeseries.tsv)")
            continue
        c = run_curve(run_dir, rank, taxon)
        last = c.iloc[-1]
        print(f"[curve] {name}: {len(c)} batches, {rank} EM {c['truth_em'].min():.4f}-"
              f"{c['truth_em'].max():.4f}, {int(last['truth_reads_cum']):,} {taxon} reads")
        curves.append(c)
    if not curves:
        raise SystemExit(f"No completed {args.library} replays")

    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / f"{args.library}.curve.tsv"
    pd.concat(curves).to_csv(out, sep="\t", index=False)
    print(f"Table → {out.resolve().relative_to(REPO)}")


if __name__ == "__main__":
    main()
