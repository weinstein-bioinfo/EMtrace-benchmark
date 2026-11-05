#!/usr/bin/env python3
"""
B_sim05_mix_background.py — Framework B genome simulation, step 5.

Two modes (--mode):

  background (or both)
    Build a background-only FASTQ by sampling from host + bacterial read
    files according to jobs/assets/background_refs.yaml.
    Read counts are taken directly from the yaml (n_reads field, Table S8).
    Output: data/petersen_2026/mixed/background_only_<N>.fastq.gz

  mixed (or both)
    For each virus × deviation, combine ALL available virus reads with
    the background-only pool (no downsampling of viral reads).
    Output: data/petersen_2026/mixed/<VIMOP_NAME>_<DEV>_mixed.fastq.gz

Ported from ViMOP's mix.sh (seqkit-based sampling) and mix_background.sh.
Does not call those scripts directly.

Run sim_bg_download.py first to populate the fasta paths in background_refs.yaml.

Usage:
  python3 jobs/framework_B/B_sim05_mix_background.py
         [--targets COVID DENV LASV_L LASV_S]
         [--mode background|mixed|both]
         [--seed 42]
         [-v]

Dependencies:
  seqkit  (brew install seqkit | conda install -c bioconda seqkit)
  python3 + pyyaml
"""

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import yaml

OUT_DIR = Path("data/petersen_2026")
FASTQ_DIR = OUT_DIR / "fastq"
MIXED_DIR = OUT_DIR / "mixed"
TARGETS_YAML = Path("jobs/framework_B/sim_targets.yaml")
BG_YAML = Path("jobs/assets/background_refs.yaml")


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def load_targets_config() -> tuple[dict, list, int]:
    """Return (vimop_names_by_name, deviations, seed) from sim_targets.yaml."""
    with open(TARGETS_YAML, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    by_name = {t["name"]: t["vimop_name"] for t in cfg["targets"]}
    deviations = cfg.get(
        "deviations", [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
    seed = cfg.get("seed", 42)
    return by_name, deviations, seed


def dev_str(deviation: float) -> str:
    """0.05 → '05', 0.30 → '30'."""
    return f"{int(round(deviation * 100)):02d}"


def load_bg_yaml(log: logging.Logger) -> tuple[list[Path], list[int]]:
    """Read background_refs.yaml; return (fastq_paths, n_reads) for present entries."""
    with open(BG_YAML, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    files: list[Path] = []
    counts: list[int] = []
    for entry in cfg["sequences"]:
        fastq = entry.get("fastq")
        if not fastq:
            log.warning("  no fastq path for %s — run B_sim04_simulate_reads.py --mode background first",
                        entry["name"])
            continue
        p = Path(fastq)
        if not p.exists():
            log.warning("  background fastq missing, skipping: %s", p)
            continue
        files.append(p)
        counts.append(int(entry["n_reads"]))
    return files, counts


# ---------------------------------------------------------------------------
# seqkit wrapper
# ---------------------------------------------------------------------------

def run(cmd: list, log: logging.Logger, **kw) -> None:
    log.debug("  $ %s", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, **kw)


def _gzip_compress(src: Path, dst: Path, log: logging.Logger) -> None:
    """Compress src → dst using pigz -p 8 if available, else gzip."""
    if shutil.which("pigz"):
        cmd = ["pigz", "-p", "8", "-c", "--", str(src)]
    else:
        cmd = ["gzip", "-c", "--", str(src)]
    log.debug("  $ %s > %s", " ".join(cmd), dst)
    with open(dst, "wb") as fout:
        subprocess.run(cmd, stdout=fout, check=True)


def sample_reads(src: Path, count: int | str, seed: int,
                 log: logging.Logger, fout) -> None:
    """
    Stream `count` reads from `src` (gz or plain) into open binary file `fout`.
    Pass count="all" to stream the entire file without seqkit.
    """
    if count == "all":
        if str(src).endswith(".gz"):
            run(["gzip", "-dc", "--", str(src)], log, stdout=fout)
        else:
            run(["cat", "--", str(src)], log, stdout=fout)
    else:
        run(["seqkit", "sample", "-2", "-n", str(count), "-s", str(seed), str(src)],
            log, stdout=fout)


# ---------------------------------------------------------------------------
# Step A: build background-only pool
# ---------------------------------------------------------------------------

def build_background(seed: int, log: logging.Logger) -> Path:
    """Sample reads from each background genome and write a pooled fastq.gz."""
    files, counts = load_bg_yaml(log)
    if not files:
        raise RuntimeError(
            f"No valid background files found in {BG_YAML}. "
            "Run sim_bg_download.py first."
        )

    bg_total = sum(counts)
    bg_out = MIXED_DIR / f"background_only_{bg_total}.fastq.gz"

    if bg_out.exists():
        log.info("  skip (exists): %s", bg_out)
        return bg_out

    log.info("--- building background-only pool (%d reads) ---", bg_total)
    log.info("  sources: %d files", len(files))

    tmp_raw = MIXED_DIR / f"background_only_{bg_total}.raw.fastq"
    with open(tmp_raw, "wb") as fout:
        for f, n in zip(files, counts):
            log.info("  sampling %d reads from %s", n, f.name)
            sample_reads(f, n, seed, log, fout)

    log.info("  compressing → %s ...", bg_out.name)
    _gzip_compress(tmp_raw, bg_out, log)
    tmp_raw.unlink()

    n_actual = _count_reads(bg_out, log)
    log.info("  done: %d reads → %s", n_actual, bg_out)
    return bg_out


# ---------------------------------------------------------------------------
# Step B: mix virus reads with background
# ---------------------------------------------------------------------------

def build_mixed(targets: list[str], vimop_names: dict, deviations: list[float],
                seed: int, bg_out: Path,
                log: logging.Logger) -> None:
    if not bg_out.exists():
        raise FileNotFoundError(
            f"Background pool not found: {bg_out}  "
            f"(run with --mode background first)"
        )

    for target in targets:
        if target not in vimop_names:
            log.warning("Unknown target %r — skipping", target)
            continue
        vname = vimop_names[target]
        for dev in deviations:
            ds = dev_str(dev)
            virus_fq = FASTQ_DIR / f"{vname}_{ds}" / f"{vname}_{ds}.fastq.gz"
            out = MIXED_DIR / f"{vname}_{ds}_mixed.fastq.gz"

            if not virus_fq.exists():
                log.warning("  virus FASTQ missing, skipping: %s", virus_fq)
                continue
            if out.exists():
                log.info("  skip (exists): %s", out.name)
                continue

            log.info("--- mixing %s_%s: all virus reads + background ---", vname, ds)

            tmp_raw = MIXED_DIR / f"{vname}_{ds}_mixed.raw.fastq"
            with open(tmp_raw, "wb") as fout:
                sample_reads(virus_fq, "all", seed, log, fout)
                sample_reads(bg_out, "all", seed, log, fout)

            _gzip_compress(tmp_raw, out, log)
            tmp_raw.unlink()

            n_reads = _count_reads(out, log)
            log.info("  done: %d reads → %s", n_reads, out)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _count_reads(fq_gz: Path, log: logging.Logger) -> int:
    """Count reads in a gzipped FASTQ by counting '@' header lines."""
    n = 0
    try:
        with gzip.open(fq_gz, "rt") as fh:
            for line in fh:
                if line.startswith("@"):
                    n += 1
    except Exception as e:
        log.warning("  could not count reads in %s: %s", fq_gz.name, e)
    return n


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Main entry point"""
    ap = argparse.ArgumentParser(
        description=textwrap.dedent(__doc__ or ""),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--targets", nargs="+",
                    default=["COVID", "DENV", "LASV_L", "LASV_S"])
    ap.add_argument("--mode", choices=["background", "mixed", "both"],
                    default="both")
    ap.add_argument("--seed", type=int, default=None,
                    help="Random seed (default: from sim_targets.yaml).")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"sim_05_{time.strftime('%Y%m%d_%H%M%S')}.log"

    fmt = logging.Formatter(
        "[%(asctime)s %(levelname)s] %(message)s", datefmt="%H:%M:%S")
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)

    sh = logging.StreamHandler()
    sh.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    sh.setFormatter(fmt)
    root.addHandler(sh)

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    root.addHandler(fh)

    log = logging.getLogger("sim05")
    log.info("log -> %s", log_file)

    if not shutil.which("seqkit"):
        log.error(
            "seqkit not found (brew install seqkit or conda install -c bioconda seqkit)")
        sys.exit(1)
    if not BG_YAML.exists():
        log.error("background refs YAML not found: %s", BG_YAML)
        sys.exit(1)

    vimop_names, default_devs, default_seed = load_targets_config()
    seed = args.seed if args.seed is not None else default_seed
    deviations = default_devs

    log.info("mode:         %s", args.mode)
    log.info("seed:         %d", seed)
    log.info("targets:      %s", args.targets)

    MIXED_DIR.mkdir(parents=True, exist_ok=True)

    bg_out = MIXED_DIR / "background_only_500000.fastq.gz"  # updated after build

    if args.mode in ("background", "both"):
        bg_out = build_background(seed, log)

    if args.mode in ("mixed", "both"):
        build_mixed(args.targets, vimop_names, deviations,
                    seed, bg_out, log)

    n_mixed = sum(1 for _ in MIXED_DIR.glob("*_mixed.fastq.gz"))
    log.info("=== sim_05 complete ===")
    log.info("  background: %s", bg_out)
    log.info("  mixed files: %d", n_mixed)
    log.info("  log: %s", log_file)


if __name__ == "__main__":
    main()
