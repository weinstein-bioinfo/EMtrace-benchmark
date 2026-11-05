#!/usr/bin/env python3
"""
D_sim06_combine_labeled.py — Framework D, step 6: labeled Nanopore datasets.

Builds the FINAL Framework D datasets: for each mutation rate the four viral
targets (COVID + DENV + LASV_L + LASV_S, each at that rate) are combined with the
human+bacteria background into one Nanopore FASTQ, together with a per-read
ground-truth table.

Reuses the existing PBSIM3-simulated reads (no re-simulation):
  viral      data/petersen_2026/fastq/<vname>_<rate>/<vname>_<rate>.fastq.gz
  background data/petersen_2026/fastq/<bg>/<bg>.fastq.gz   (per-organism!)

Per-organism background sampling is essential: the pooled background_only_*.fastq.gz
loses organism identity, so it cannot be labeled human vs bacteria. Each organism
is sampled to its yaml read count (seqkit sample, fixed seed) and labeled as it is
added.

Outputs (data/petersen_2026/combined/):
  combined_<rate>.fastq.gz     COVID+DENV+LASV_L+LASV_S (@rate) + background
  combined_<rate>.truth.tsv    readid  division  species  taxid  source  mutation_rate
  combined.summary.tsv         rate  n_viral  n_bacteria  n_eukaryota  n_total

Labels come from assets/source_taxonomy.tsv (run D_sim00_taxonomy.py first).

Usage:
  python3 jobs/framework_D/D_sim06_combine_labeled.py
         [--rates 00 10 20 30] [--targets COVID DENV LASV_L LASV_S]
         [--seed 42] [-v]

Dependencies: seqkit, pigz|gzip, python3 + pyyaml.
"""

import argparse
import logging
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

OUT_DIR = Path("data/petersen_2026")
FASTQ_DIR = OUT_DIR / "fastq"
COMBINED_DIR = OUT_DIR / "combined"
TARGETS_YAML = Path("jobs/framework_B/sim_targets.yaml")
BG_YAML = Path("jobs/assets/background_refs.yaml")
TAXONOMY_TSV = Path("jobs/assets/source_taxonomy.tsv")

DEFAULT_RATES = ["00", "10", "20", "30"]
DEFAULT_TARGETS = ["COVID", "DENV", "LASV_L", "LASV_S"]


# ---------------------------------------------------------------------------
# Config / metadata
# ---------------------------------------------------------------------------

def load_taxonomy(log: logging.Logger) -> dict[str, dict]:
    """source → {division, species, taxid} from source_taxonomy.tsv."""
    if not TAXONOMY_TSV.exists():
        log.error("missing %s — run D_sim00_taxonomy.py first", TAXONOMY_TSV)
        sys.exit(1)
    tax: dict[str, dict] = {}
    with open(TAXONOMY_TSV, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tax[f[idx["source"]]] = {
                "division": f[idx["division"]],
                "species": f[idx["species"]],
                "taxid": f[idx["taxid"]],
            }
    return tax


def load_targets(log: logging.Logger) -> tuple[dict[str, str], int]:
    """Return (name → vimop_name, seed) from sim_targets.yaml."""
    with open(TARGETS_YAML, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    by_name = {t["name"]: t["vimop_name"] for t in cfg["targets"]}
    return by_name, int(cfg.get("seed", 42))


def load_background(log: logging.Logger) -> list[tuple[str, Path, int]]:
    """Return [(name, fastq_path, n_reads)] for background organisms present."""
    with open(BG_YAML, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    out: list[tuple[str, Path, int]] = []
    for entry in cfg["sequences"]:
        fq = entry.get("fastq")
        if not fq or not Path(fq).exists():
            log.warning("  background fastq missing, skipping: %s", entry["name"])
            continue
        out.append((entry["name"], Path(fq), int(entry["n_reads"])))
    return out


# ---------------------------------------------------------------------------
# Shell helpers (mirrors B_sim05)
# ---------------------------------------------------------------------------

def _run(cmd: list, log: logging.Logger, **kw) -> None:
    log.debug("  $ %s", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, **kw)


def _gzip_compress(src: Path, dst: Path, log: logging.Logger) -> None:
    cmd = (["pigz", "-p", "8", "-c", "--", str(src)] if shutil.which("pigz")
           else ["gzip", "-c", "--", str(src)])
    with open(dst, "wb") as fout:
        subprocess.run(cmd, stdout=fout, check=True)


def _read_ids(fastq: Path, log: logging.Logger) -> list[str]:
    """Read ids (first header token) of a FASTQ via `seqkit seq -n -i`."""
    res = subprocess.run(["seqkit", "seq", "-n", "-i", str(fastq)],
                         check=True, capture_output=True, text=True)
    return [ln for ln in res.stdout.splitlines() if ln]


def _stream_all(src: Path, fout, log: logging.Logger) -> None:
    """Append the whole FASTQ (gz or plain) to open binary handle fout."""
    cmd = ["gzip", "-dc", "--", str(src)] if str(src).endswith(".gz") else ["cat", "--", str(src)]
    _run(cmd, log, stdout=fout)


def _sample(src: Path, n: int, seed: int, dst: Path, log: logging.Logger) -> None:
    """Write n seqkit-sampled reads from src into plain FASTQ dst."""
    with open(dst, "wb") as fout:
        _run(["seqkit", "sample", "-2", "-n", str(n), "-s", str(seed), str(src)],
             log, stdout=fout)


# ---------------------------------------------------------------------------
# Build one combined dataset
# ---------------------------------------------------------------------------

def build_rate(rate: str, targets: list[str], by_name: dict[str, str],
               bg: list[tuple[str, Path, int]], tax: dict[str, dict],
               seed: int, log: logging.Logger) -> dict:
    out_fq = COMBINED_DIR / f"combined_{rate}.fastq.gz"
    out_truth = COMBINED_DIR / f"combined_{rate}.truth.tsv"
    if out_fq.exists() and out_truth.exists():
        log.info("  skip (exists): %s", out_fq.name)
        return {}

    tmp_raw = COMBINED_DIR / f"combined_{rate}.raw.fastq"
    truth_rows: list[tuple] = []          # (readid, division, species, taxid, source)
    counts = {"Viruses": 0, "Bacteria": 0, "Eukaryota": 0}

    def _label(source: str, ids: list[str]) -> None:
        meta = tax.get(source)
        if meta is None:
            raise KeyError(f"no taxonomy for source {source!r} (run D_sim00_taxonomy.py)")
        for rid in ids:
            truth_rows.append((rid, meta["division"], meta["species"], meta["taxid"], source))
        counts[meta["division"]] = counts.get(meta["division"], 0) + len(ids)

    with open(tmp_raw, "wb") as fout:
        # --- viral: all reads per target at this rate ---
        for name in targets:
            vname = by_name.get(name)
            if vname is None:
                log.warning("  unknown target %r — skipping", name)
                continue
            src = FASTQ_DIR / f"{vname}_{rate}" / f"{vname}_{rate}.fastq.gz"
            if not src.exists():
                log.warning("  viral FASTQ missing, skipping: %s", src)
                continue
            log.info("  + viral %s (%s) — all reads", name, src.name)
            ids = _read_ids(src, log)
            _stream_all(src, fout, log)
            _label(name, ids)

        # --- background: per-organism sampling (labelable) ---
        for bname, bpath, n in bg:
            log.info("  + background %s — %d reads", bname, n)
            tmp_s = COMBINED_DIR / f"_bg_{rate}_{bname}.fastq"
            _sample(bpath, n, seed, tmp_s, log)
            ids = _read_ids(tmp_s, log)
            with open(tmp_s, "rb") as fin:
                shutil.copyfileobj(fin, fout)
            tmp_s.unlink()
            _label(bname, ids)

    log.info("  compressing → %s", out_fq.name)
    _gzip_compress(tmp_raw, out_fq, log)
    tmp_raw.unlink()

    with open(out_truth, "w", encoding="utf-8") as fh:
        fh.write("readid\tdivision\tspecies\ttaxid\tsource\tmutation_rate\n")
        for rid, div, sp, txid, source in truth_rows:
            fh.write(f"{rid}\t{div}\t{sp}\t{txid}\t{source}\t{rate}\n")

    n_total = len(truth_rows)
    log.info("  done: %d reads (%s) → %s", n_total, dict(counts), out_fq.name)
    return {"rate": rate, "n_viral": counts.get("Viruses", 0),
            "n_bacteria": counts.get("Bacteria", 0),
            "n_eukaryota": counts.get("Eukaryota", 0), "n_total": n_total}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rates", nargs="+", default=DEFAULT_RATES)
    ap.add_argument("--targets", nargs="+", default=DEFAULT_TARGETS)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="[%(levelname)s] %(message)s")
    log = logging.getLogger("D_sim06")

    if not shutil.which("seqkit"):
        log.error("seqkit not found (brew install seqkit | conda install -c bioconda seqkit)")
        sys.exit(1)

    tax = load_taxonomy(log)
    by_name, default_seed = load_targets(log)
    bg = load_background(log)
    seed = args.seed if args.seed is not None else default_seed

    COMBINED_DIR.mkdir(parents=True, exist_ok=True)
    log.info("rates=%s targets=%s seed=%d background=%d organisms",
             args.rates, args.targets, seed, len(bg))

    summary: list[dict] = []
    for rate in args.rates:
        log.info("=== combined_%s ===", rate)
        row = build_rate(rate, args.targets, by_name, bg, tax, seed, log)
        if row:
            summary.append(row)

    if summary:
        sfile = COMBINED_DIR / "combined.summary.tsv"
        cols = ["rate", "n_viral", "n_bacteria", "n_eukaryota", "n_total"]
        write_header = not sfile.exists()
        with open(sfile, "a", encoding="utf-8") as fh:
            if write_header:
                fh.write("\t".join(cols) + "\n")
            for r in summary:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")
        log.info("summary → %s", sfile)

    log.info("=== D_sim06 complete: %d dataset(s) ===", len(summary))


if __name__ == "__main__":
    main()
