#!/usr/bin/env python3
"""Extract Framework-B metrics from INSaFLU/TELEVIR combined workflow reports.

For every sample under results/insaflu_out/<set>/<sample>/, read the
`<sample>_cr_date_*.tsv` (Validation Mapping Results), find the row for the
sample's TARGET virus, and emit:

    coverage %  = `Cov (%)`
    identity %  = 100 - `Mismatch rate`   (Mismatch rate is a percentage)

Samples where TELEVIR produced no target mapping (empty report, or no target
row among background hits) are emitted as coverage=0 / no consensus, matching
the "no consensus when recovery < 50%" convention used for ViMOP.

Column layout differs between reports (the no-consensus TSVs even lack the
`Mismatch rate` column), so columns are located BY HEADER NAME, not position.
Cells may carry leading/trailing whitespace and are stripped.

This is the TIER-1 (reported) metric: the platform's self-reported statistics vs
TELEVIR's own auto-selected reference. The companion TIER-2 (consensus) metric —
consensus rebuilt from the BAM and scored against the original benchmark
reference — is produced by jobs/framework_B.televir.metrics.bash and written to
televir_<label>.* (no "reported" marker), a drop-in peer of vimop_<label>.*.

Outputs (metrics/framework_B/), labels matching ViMOP (pure, mixed,
peterson_2026, gauthier_2021) so tier-1 pairs with tier-2 by label:
    televir_reported.<label>.summary.tsv   one per benchmark set
    televir_reported.all.summary.tsv       combined

Run from the repo root:
    python3 jobs/framework_B/extract_televir_metrics.py
"""

from __future__ import annotations

import csv
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))          # .../jobs/framework_B
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
INSAFLU_DIR = os.path.join(REPO, "results", "insaflu_out")
OUT_DIR = os.path.join(REPO, "metrics", "framework_B")

SETS = ["pure_mutants", "mixed_mutants", "petersen2026", "gauthier2021"]

# benchmark-set dir -> ViMOP label (so tier-1 files pair with tier-2 by label)
SET_LABEL = {
    "pure_mutants": "pure",
    "mixed_mutants": "mixed",
    "petersen2026": "peterson_2026",
    "gauthier2021": "gauthier_2021",
}

# Target-virus matcher: (description regex, taxid set). A row matches the target
# if its Description matches the regex OR its Taxid is in the taxid set.
TARGETS = {
    "COVID":  (re.compile(r"coronavirus 2\b|SARS-CoV-2", re.I), {"2697049"}),
    "DENV":   (re.compile(r"dengue", re.I), {"12637", "11053"}),
    "LASV_L": (re.compile(r"lassa.*segment\s*l|segment\s*l.*lassa|lassa.*\bL\b", re.I), set()),
    "LASV_S": (re.compile(r"lassa.*segment\s*s|segment\s*s.*lassa|lassa.*\bS\b", re.I), set()),
    "Zika":   (re.compile(r"zika", re.I), set()),
    "WestNile": (re.compile(r"west\s*nile", re.I), set()),
    "SARS":   (re.compile(r"coronavirus 2\b|SARS-CoV-2", re.I), {"2697049"}),
}


def target_for(setname: str, sample: str) -> str:
    """Map a sample id to its target-virus key."""
    if setname in ("pure_mutants", "mixed_mutants"):
        base = sample.replace("_mixed", "")
        for key in ("LASV_L", "LASV_S", "COVID", "DENV"):
            if base.startswith(key):
                return key
    elif setname == "petersen2026":
        prefix = sample.split("_", 1)[0]
        return {"Zika": "Zika", "Dengue": "DENV",
                "WestNile": "WestNile", "SARS": "SARS"}[prefix]
    elif setname == "gauthier2021":
        return "SARS"
    raise ValueError(f"cannot map target for {setname}/{sample}")


def read_report(tsv_path: str):
    """Yield dict rows (stripped) from a cr_date TSV, keyed by stripped header."""
    with open(tsv_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = [h.strip() for h in next(reader)]
        except StopIteration:
            return
        for raw in reader:
            if not any(c.strip() for c in raw):
                continue
            cells = [c.strip() for c in raw]
            cells += [""] * (len(header) - len(cells))
            yield dict(zip(header, cells))


def pick_target_row(rows, target_key):
    """Return the target row with the highest Cov (%), or None."""
    rx, taxids = TARGETS[target_key]
    best, best_cov = None, -1.0
    for r in rows:
        desc, taxid = r.get("Description", ""), r.get("Taxid", "")
        if not (rx.search(desc) or taxid in taxids):
            continue
        try:
            cov = float(r.get("Cov (%)", "") or "nan")
        except ValueError:
            cov = float("nan")
        if cov == cov and cov > best_cov:      # skip NaN
            best, best_cov = r, cov
    return best


def process_sample(setname: str, sample: str, sdir: str) -> dict:
    target = target_for(setname, sample)
    tsvs = glob.glob(os.path.join(sdir, "*cr_date*.tsv"))
    has_bam = bool(glob.glob(os.path.join(sdir, "*.bam")))

    # sample_id joins with tier-2 / vimop / emtrace: accession for sra, dir for mutants
    if setname in ("petersen2026", "gauthier2021"):
        sample_id = sample.split("_", 1)[1]
    else:
        sample_id = sample

    rec = {"set": setname, "label": SET_LABEL[setname],
           "sample": sample, "sample_id": sample_id, "target": target,
           "accID": "", "mean_coverage_pct": 0.0, "mean_identity_pct": "",
           "mismatch_rate": "", "mapped_reads": "", "depth": "",
           "consensus": "no", "has_bam": "yes" if has_bam else "no"}

    if not tsvs:
        rec["note"] = "no report tsv"
        return rec

    rows = list(read_report(tsvs[0]))
    if not rows:
        rec["note"] = "empty report (no mapping)"
        return rec

    row = pick_target_row(rows, target)
    if row is None:
        rec["note"] = "target not detected"
        return rec

    cov = float(row["Cov (%)"])
    rec["accID"] = row.get("accID", "")
    rec["mean_coverage_pct"] = round(cov, 2)
    rec["mapped_reads"] = row.get("Mapped reads", "")
    rec["depth"] = row.get("Depth", "")
    mm = row.get("Mismatch rate", "")
    if mm:
        rec["mismatch_rate"] = mm
        rec["mean_identity_pct"] = round(100.0 - float(mm), 2)
    rec["consensus"] = "yes"
    rec["note"] = ""
    return rec


def main() -> int:
    os.makedirs(OUT_DIR, exist_ok=True)
    cols = ["set", "label", "sample", "sample_id", "target", "accID",
            "mean_coverage_pct", "mean_identity_pct", "mismatch_rate",
            "mapped_reads", "depth", "consensus", "has_bam", "note"]
    combined = []

    for setname in SETS:
        setdir = os.path.join(INSAFLU_DIR, setname)
        if not os.path.isdir(setdir):
            continue
        recs = []
        for sample in sorted(os.listdir(setdir)):
            sdir = os.path.join(setdir, sample)
            if not os.path.isdir(sdir):
                continue
            recs.append(process_sample(setname, sample, sdir))

        label = SET_LABEL[setname]
        out = os.path.join(OUT_DIR, f"televir_reported.{label}.summary.tsv")
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(recs)
        n_cons = sum(r["consensus"] == "yes" for r in recs)
        print(f"{setname:14s} {len(recs):2d} samples  ({n_cons} with consensus)  -> {out}")
        combined.extend(recs)

    allout = os.path.join(OUT_DIR, "televir_reported.all.summary.tsv")
    with open(allout, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(combined)
    print(f"{'ALL':14s} {len(combined):2d} samples -> {allout}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
