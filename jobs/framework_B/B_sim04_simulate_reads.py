#!/usr/bin/env python3
"""
B_sim04_simulate_reads.py — Framework B genome simulation, step 4.

Simulates Nanopore reads using PBSIM3 (SISPA profile) for two sets of inputs:

  viral (--mode viral or both)
    One FASTQ per mutant genome produced by sim_03.
    Output: data/petersen_2026/fastq/<VIMOP_NAME>_<DEV>/<VIMOP_NAME>_<DEV>.fastq.gz

  background (--mode background or both)
    One FASTQ per organism in background_refs.yaml (run sim_bg_download.py first).
    Depth is computed per genome so that at least n_reads × 2 reads are produced;
    sim_05 then samples the exact count with seqkit.
    Output: data/petersen_2026/fastq/<NAME>/<NAME>.fastq.gz
    Side-effect: writes the fastq path back into background_refs.yaml (fastq field).

MinION header conversion (gzopen, read_fastq, minionize) ported from
minionize.py in the ViMOP benchmark repository.

Usage:
  python3 jobs/framework_B/B_sim04_simulate_reads.py
         [--targets COVID DENV LASV_L LASV_S]
         [--mode viral|background|both]
         [--seed 42] [--depth 30] [--threads 8] [-v]

Dependencies:
  pbsim3         (conda install -c bioconda pbsim3)
  python3 + biopython + pyyaml
  PBSIM3 profile files in jobs/assets/pbsim_profile/
  background_refs.yaml fasta fields populated by sim_bg_download.py (background mode)
"""

import argparse
import gzip
import logging
import math
import random
import shutil
import subprocess
import sys
import tempfile
import textwrap
import time
import uuid
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import IO, cast

import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

OUT_DIR = Path("data/petersen_2026")
MUTANTS_DIR = OUT_DIR / "mutants"
FASTQ_DIR = OUT_DIR / "fastq"
TARGETS_YAML = Path("jobs/framework_B/sim_targets.yaml")
BG_YAML = Path("jobs/assets/background_refs.yaml")
ASSETS_DIR = Path("jobs/assets")
PROFILE_DIR = ASSETS_DIR / "pbsim_profile"

PROFILE_ID = "sispatrimmed"
MODEL_VERSION = "dna_r10.4.1_e8.2_400bps_hac@v3.5.2"
BARCODE = "barcode01"
MEAN_READ_LEN = 666   # SISPA profile mean read length (bp); from sample_profile_sispatrimmed.stats


# ---------------------------------------------------------------------------
# MinION header conversion — ported from ViMOP minionize.py
# ---------------------------------------------------------------------------

def gzopen(fname: str | Path, mode: str = "r") -> IO[str]:
    fname = str(fname)
    if fname.endswith(".gz"):
        return cast(IO[str], gzip.open(fname, mode + "t", encoding="utf-8"))
    return open(fname, mode, encoding="utf-8")


def read_fastq(file_path: Path):
    """Yield (header, sequence, comment, quality) tuples from a FASTQ file."""
    with gzopen(file_path) as fh:
        while True:
            header = fh.readline().strip()
            if not header:
                break
            seq = fh.readline().strip()
            comment = fh.readline().strip()
            quality = fh.readline().strip()
            if not (header.startswith("@") and comment.startswith("+")):
                raise ValueError(
                    f"Invalid FASTQ format at read starting: {header!r}")
            yield header[1:], seq, comment[1:], quality


def minionize(input_path: Path, output_path: Path,
              model_version_id: str, barcode: str) -> int:
    """
    Rewrite FASTQ headers to MinION-style UUIDs with run metadata.
    Returns the number of reads written.
    """
    runid = uuid.uuid4().hex
    start_time = datetime.now(timezone.utc)
    n = 0

    with gzopen(output_path, "w") as fout:
        for i, (_, seq, _, quality) in enumerate(read_fastq(input_path), 1):
            read_uuid = str(uuid.uuid4())
            ch = random.randint(100, 512)
            ts = start_time + timedelta(seconds=i)
            desc = (f"runid={runid} sampleid=no_sample read={i} ch={ch} "
                    f"start_time={ts.isoformat()}Z "
                    f"model_version_id={model_version_id} barcode={barcode}")
            record = SeqRecord(
                Seq(seq), id=read_uuid, name=read_uuid, description=desc,
                letter_annotations={"phred_quality": [
                    ord(q) - 33 for q in quality]},
            )
            SeqIO.write(record, fout, "fastq")
            n += 1
    return n


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def load_targets_config() -> tuple[dict, list, int]:
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


# ---------------------------------------------------------------------------
# PBSIM3 wrapper
# ---------------------------------------------------------------------------

def run(cmd: list, log: logging.Logger, **kw) -> subprocess.CompletedProcess:
    log.debug("  $ %s", " ".join(str(c) for c in cmd))
    return subprocess.run(cmd, check=True, **kw)


MAX_PBSIM_LEN = 900_000_000  # PBSIM3 hard limit is 1 Gb per reference record


def _chunk_fasta(src: Path, tmp: Path) -> list[Path]:
    """
    Read a (possibly multi-contig) FASTA and write it as one or more chunk FASTAs,
    each containing a single N-joined record ≤ MAX_PBSIM_LEN bp.

    Fixes two PBSIM3 limitations:
      - SIGPIPE: PBSIM3 opens one gzip pipe per contig; many-contig assemblies
        exhaust internal pipe limits.
      - Hard 1 Gb reference limit: genomes > 1 Gb (e.g. human) must be split.

    Returns list of chunk FASTA paths (usually one for small genomes).
    """
    contigs: list[str] = []
    current: list[str] = []
    with open(src) as fh:
        for line in fh:
            if line.startswith(">"):
                if current:
                    contigs.append("".join(current))
                current = []
            else:
                current.append(line.strip())
    if current:
        contigs.append("".join(current))

    N_JOIN = "N" * 100
    chunks: list[list[str]] = []
    chunk: list[str] = []
    chunk_len = 0
    for seq in contigs:
        extra = len(seq) + (len(N_JOIN) if chunk else 0)
        if chunk and chunk_len + extra > MAX_PBSIM_LEN:
            chunks.append(chunk)
            chunk = [seq]
            chunk_len = len(seq)
        else:
            chunk.append(seq)
            chunk_len += extra
    if chunk:
        chunks.append(chunk)

    paths: list[Path] = []
    for i, seqs in enumerate(chunks):
        joined = N_JOIN.join(seqs)
        p = tmp / f"chunk_{i:04d}.fna"
        with open(p, "w") as fh:
            fh.write(f">chunk_{i}\n")
            for j in range(0, len(joined), 60):
                fh.write(joined[j:j + 60] + "\n")
        paths.append(p)
    return paths


def _pbsim_run_and_collect(chunk: Path, prefix: Path, depth: int, seed: int,
                           log: logging.Logger) -> list[Path]:
    """Run PBSIM3 on a single chunk FASTA; return verified .fq.gz output paths."""
    run([
        "pbsim",
        "--strategy",         "wgs",
        "--method",           "sample",
        "--depth",            str(depth),
        "--seed",             str(seed),
        "--genome",           str(chunk.resolve()),
        "--difference-ratio", "39:24:36",
        "--sample-profile-id", PROFILE_ID,
        "--prefix",           str(prefix),
    ], log, cwd=str(PROFILE_DIR))

    raw_files = sorted(prefix.parent.glob(f"{prefix.name}_*.fq.gz"))
    if not raw_files:
        raise RuntimeError(f"PBSIM3 produced no output for chunk {chunk.name}")

    # PBSIM3 uses fclose() instead of pclose() on its internal gzip pipe,
    # so the gzip child may still be flushing when pbsim exits.  Poll each
    # output file until gzip -t confirms it is complete (max 10 s).
    for rf in raw_files:
        for _ in range(100):
            r = subprocess.run(["gzip", "-t", str(rf)], capture_output=True)
            if r.returncode == 0:
                break
            time.sleep(0.1)
        else:
            raise RuntimeError(f"PBSIM3 output incomplete after 10s: {rf}")
    return raw_files


def simulate_one(genome_fa: Path, out_fq: Path, depth: int, seed: int,
                 log: logging.Logger) -> None:
    """Run PBSIM3 on genome_fa and write MinION-formatted reads to out_fq."""
    if out_fq.exists():
        log.info("  skip (exists): %s", out_fq.name)
        return

    if not genome_fa.exists():
        log.warning("  genome not found, skipping: %s", genome_fa)
        return

    out_fq.parent.mkdir(parents=True, exist_ok=True)
    log.info("  simulating reads: %s (depth=%dx) ...", genome_fa.name, depth)

    with tempfile.TemporaryDirectory(prefix="sim04_") as tmp:
        tmp = Path(tmp)

        # PBSIM3 segfaults on gzipped input — decompress first.
        if str(genome_fa).endswith(".gz"):
            genome_plain = tmp / genome_fa.stem
            log.debug("  decompressing %s ...", genome_fa.name)
            with gzip.open(genome_fa, "rb") as fin, open(genome_plain, "wb") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            genome_plain = genome_fa

        # Split into ≤ 900 Mb chunks (fixes SIGPIPE on many-contig assemblies
        # and the 1 Gb hard limit for large genomes like human).
        chunks = _chunk_fasta(genome_plain, tmp)
        if len(chunks) > 1:
            log.debug("  split into %d chunks (≤ 900 Mb each)", len(chunks))

        all_raw: list[Path] = []
        for i, chunk in enumerate(chunks):
            prefix = tmp / f"raw_c{i:04d}"
            raw_files = _pbsim_run_and_collect(chunk, prefix, depth, seed, log)
            all_raw.extend(raw_files)

        concat = tmp / "concat.fastq.gz"
        with open(concat, "wb") as fout:
            for rf in all_raw:
                fout.write(rf.read_bytes())

        # Convert to MinION headers
        n = minionize(concat, out_fq, MODEL_VERSION, BARCODE)

    log.info("  done: %d reads → %s", n, out_fq)


# ---------------------------------------------------------------------------
# Background read simulation
# ---------------------------------------------------------------------------

def genome_length(fasta: Path) -> int:
    """Sum of all contig lengths in a FASTA file (gz or plain)."""
    total = 0
    opener = gzip.open if str(fasta).endswith(".gz") else open
    with opener(fasta, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def simulate_background(seed: int, log: logging.Logger) -> None:
    """Simulate reads from each background reference genome in background_refs.yaml.

    For each entry the sequencing depth is chosen so that PBSIM3 produces at
    least n_reads * 2 reads; sim_05 then samples the exact n_reads count with
    seqkit.  The fastq output path is written back into background_refs.yaml
    under the fastq field so sim_05 can find it.
    """
    if not BG_YAML.exists():
        log.error("background_refs.yaml not found: %s — run sim_bg_download.py first", BG_YAML)
        return

    with open(BG_YAML, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    changed = False
    for entry in cfg["sequences"]:
        name    = entry["name"]
        fasta   = entry.get("fasta")
        n_reads = int(entry["n_reads"])
        out_fq  = FASTQ_DIR / name / f"{name}.fastq.gz"

        if not fasta:
            log.warning("  %s: no fasta path — run sim_bg_download.py first", name)
            continue
        fasta_path = Path(fasta)
        if not fasta_path.exists():
            log.warning("  %s: fasta missing: %s", name, fasta_path)
            continue

        if entry.get("fastq") and Path(entry["fastq"]).exists():
            log.info("  %s: skip (fastq exists)", name)
        else:
            log.info("=== background %s — target %d reads ===", name, n_reads)
            log.info("  measuring genome length from %s ...", fasta_path.name)
            glen  = genome_length(fasta_path)
            depth = max(1, math.ceil(n_reads * MEAN_READ_LEN * 2.0 / glen))
            log.info("  genome: %d bp → depth=%dx (~%d reads)",
                     glen, depth, int(depth * glen / MEAN_READ_LEN))
            simulate_one(fasta_path, out_fq, depth, seed, log)

        if str(out_fq) != entry.get("fastq"):
            entry["fastq"] = str(out_fq)
            changed = True

    if changed:
        with open(BG_YAML, "w", encoding="utf-8") as fh:
            yaml.dump(cfg, fh, allow_unicode=True, sort_keys=False)
        log.info("updated %s with fastq paths", BG_YAML)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    "Main entry point"
    ap = argparse.ArgumentParser(
        description=textwrap.dedent(__doc__ or ""),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--targets", nargs="+",
                    default=["COVID", "DENV", "LASV_L", "LASV_S"],
                    help="Viral targets to simulate (default: all).")
    ap.add_argument("--mode", choices=["viral", "background", "both"],
                    default="both",
                    help="What to simulate: viral mutants, background genomes, or both (default: both).")
    ap.add_argument("--seed",    type=int, default=None,
                    help="Random seed (default: from sim_targets.yaml).")
    ap.add_argument("--depth",   type=int, default=30,
                    help="Sequencing depth for viral reads (default: 30). Background depth is computed per genome.")
    ap.add_argument("--threads", type=int, default=8,
                    help="PBSIM3 threads (default: 8).")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"sim_04_{time.strftime('%Y%m%d_%H%M%S')}.log"

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

    log = logging.getLogger("sim04")
    log.info("log -> %s", log_file)

    # Dependency checks
    if not shutil.which("pbsim"):
        log.error("pbsim not found (conda install -c bioconda pbsim3)")
        sys.exit(1)
    if not PROFILE_DIR.exists():
        log.error("PBSIM3 profile directory not found: %s", PROFILE_DIR)
        sys.exit(1)

    vimop_names, default_devs, default_seed = load_targets_config()
    deviations = default_devs
    seed = args.seed if args.seed is not None else default_seed

    log.info("mode:       %s",   args.mode)
    log.info("depth:      %dx",  args.depth)
    log.info("seed:       %d",   seed)
    log.info("targets:    %s",   args.targets)
    log.info("deviations: %s",   deviations)

    FASTQ_DIR.mkdir(parents=True, exist_ok=True)

    if args.mode in ("viral", "both"):
        for target in args.targets:
            if target not in vimop_names:
                log.warning("Unknown target %r — skipping", target)
                continue
            vname = vimop_names[target]
            log.info("=== viral %s ===", vname)
            for dev in deviations:
                ds = dev_str(dev)
                genome = MUTANTS_DIR / vname / f"{vname}_{ds}.fasta"
                out_fq = FASTQ_DIR / f"{vname}_{ds}" / f"{vname}_{ds}.fastq.gz"
                simulate_one(genome, out_fq, args.depth, seed, log)

    if args.mode in ("background", "both"):
        simulate_background(seed, log)

    n_fq = sum(1 for _ in FASTQ_DIR.rglob("*.fastq.gz"))
    log.info("done — %d fastq.gz files in %s", n_fq, FASTQ_DIR)
    log.info("log -> %s", log_file)


if __name__ == "__main__":
    main()
