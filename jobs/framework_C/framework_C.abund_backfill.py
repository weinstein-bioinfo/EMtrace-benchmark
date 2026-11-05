#!/usr/bin/env python3
"""Framework C — back-fill the per-batch abundance tables of finished replays.

Module C now keeps every batch's cumulative `abund_table.tsv` as
`<run>/abund/batch_NNNN.tsv`, so a replay can be scored again
(`ngs-benchmark realtime-rescore`) without replaying it. Replays made before
that keep only the final table. This script rebuilds the missing ones from the
per-batch annotations, with emtrace's own code, exactly as `emtrace watch`
built them live:

  - each batch's annot2.tsv (else annot1.tsv) is appended to the accumulated
    annotations, and the warm-started EM is run on them with the run's config;
  - the abundance table is then built by `BatchProcessor._write_abundance_tsv`
    (taxonomy linkage, exclude_no_species and the tail trim included).

The last rebuilt table is checked against the replay's own final
`abund_table.tsv` (taxa and readn), and the script fails on a mismatch.

Must run with the Python that has emtrace installed:

    PY=$(head -1 "$(command -v emtrace)" | sed -E 's/^#! *//')
    $PY jobs/framework_C/framework_C.abund_backfill.py \\
        ticks_80-43.refseq_abv1g_k11w5 ticks_35-33.refseq_abv1g_k11w5_dmnd

Then rescore, e.g.:

    ngs-benchmark realtime-rescore --outdir results/framework_C/<run> \\
        --library <run> --rank genus --read-min 10 --stabilise-on em_prob \\
        --stabilise-tol 0.01 --panel jobs/framework_C/panels/<pool>.panel.tsv
"""

from __future__ import annotations

import argparse
import shutil
import tempfile
from pathlib import Path

import pandas as pd
from emtrace.core.config import LazypipeConfig
from emtrace.models.annotation import read_annotation_tsv
from emtrace.realtime.processor import BatchProcessor

REPO = Path(__file__).resolve().parents[2]   # jobs/<framework>/<script>
RESULTS = REPO / "results/framework_C"


def read_provenance(path: Path) -> dict:
    prov = {}
    for line in path.read_text().splitlines():
        key, _, val = line.partition(":")
        prov[key.strip()] = val.strip()
    return prov


def batch_annotation(bdir: Path) -> pd.DataFrame:
    """The table `watch` merged for this batch: annot2.tsv, else annot1.tsv."""
    for name in ("annot2.tsv", "annot1.tsv"):
        path = bdir / "annotation" / name
        if path.is_file():
            return read_annotation_tsv(path)
    return pd.DataFrame()


def pipeline_dir(run_dir: Path) -> Path:
    """`pipeline/<job>`: named after the run, or a copied run's original job (a replay copied to a new name)."""
    named = run_dir / "pipeline" / run_dir.name
    if named.is_dir():
        return named
    found = [p for p in (run_dir / "pipeline").iterdir() if (p / "batches").is_dir()]
    if len(found) != 1:
        raise SystemExit(f"{run_dir.name}: expected one pipeline/<job>/batches, found {len(found)}")
    return found[0]


def same_table(a: pd.DataFrame, b: pd.DataFrame) -> bool:
    """Same taxa with the same read counts (row order and float noise aside)."""
    ka = a.set_index(a["taxid"].astype(str))["readn"].sort_index()
    kb = b.set_index(b["taxid"].astype(str))["readn"].sort_index()
    return ka.index.equals(kb.index) and (ka.values == kb.values).all()


def backfill(run_dir: Path, force: bool) -> None:
    library = run_dir.name
    out = run_dir / "abund"
    if out.is_dir() and any(out.iterdir()) and not force:
        print(f"[skip] {library}: abund/ already filled (use --force)")
        return
    if not (run_dir / "timeseries.tsv").is_file():
        print(f"[skip] {library}: no timeseries.tsv (not scored)")
        return
    pdir = pipeline_dir(run_dir)
    batches = sorted((pdir / "batches").glob("batch_*"))
    ts = pd.read_csv(run_dir / "timeseries.tsv", sep="\t")
    if len(batches) != len(ts):
        raise SystemExit(f"{library}: {len(batches)} batch dirs but {len(ts)} timeseries rows")

    cfg = LazypipeConfig.from_yaml(REPO / read_provenance(run_dir / "provenance.txt")["config"])
    out.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        proc = BatchProcessor(cfg, Path(tmp), job_name=pdir.name)
        table = Path(tmp) / "results" / "abund_table.tsv"
        for bdir in batches:
            index = int(bdir.name.split("_")[1])
            annot = batch_annotation(bdir)
            em = proc._run_em(annot)
            proc.state.update(annot, em, 0, 1)
            if proc.state.F:
                proc._write_em_tsv(proc.state.F, proc._em_prob_path, proc._em_n_queries())
            table.unlink(missing_ok=True)
            proc._write_abundance_tsv()
            if table.is_file():
                shutil.copyfile(table, out / f"batch_{index:04d}.tsv")
            print(f"  {library} batch {index}: "
                  f"{'%d taxa' % len(pd.read_csv(table, sep=chr(9))) if table.is_file() else 'no table'}",
                  flush=True)

    last = out / f"batch_{int(batches[-1].name.split('_')[1]):04d}.tsv"
    final = pd.read_csv(pdir / "results" / "abund_table.tsv", sep="\t", dtype={"taxid": str})
    if not same_table(pd.read_csv(last, sep="\t", dtype={"taxid": str}), final):
        raise SystemExit(f"{library}: rebuilt last batch differs from the replay's abund_table.tsv")
    print(f"[ok] {library}: {len(batches)} tables; last matches the replay's final table")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("runs", nargs="+", help="replay directory names under --results")
    ap.add_argument("--results", type=Path, default=RESULTS)
    ap.add_argument("--force", action="store_true", help="rebuild an existing abund/")
    args = ap.parse_args()
    for name in args.runs:
        backfill(args.results / name, args.force)


if __name__ == "__main__":
    main()
