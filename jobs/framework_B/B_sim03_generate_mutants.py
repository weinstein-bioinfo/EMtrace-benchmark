#!/usr/bin/env python3
"""
B_sim03_generate_mutants.py — Framework B genome simulation, step 3.

For each of 4 target viruses × 7 deviation fractions (0.00–0.30):
  1. Load the reference sequence (template genome).
  2. Load positional variability weights produced by sim_02.
  3. Determine how many sites to mutate:
       N = ceil(len_ref × deviation)
  4. Sample N positions without replacement, weighted by variability weights
     (higher weight → higher probability of selection — variable sites preferred).
  5. At each selected position, replace the nucleotide with one of the three
     alternative bases chosen uniformly at random.
  6. Write the mutated genome to:
       data/petersen_2026/mutants/<NAME>/<NAME>_<DEV>.fasta
     e.g.  mutants/COVID/COVID_00.fasta … mutants/LASV_S/LASV_S_30.fasta

Deviation 0.00 produces an unmutated copy of the reference (baseline).

Outputs (per target × deviation):
  mutants/<NAME>/<NAME>_<DEV>.fasta
Summary TSV written to mutants/summary.tsv.

Naming convention matches ViMOP's sim_virus.sh for downstream compatibility.

Dependencies: python3 + biopython, numpy, pyyaml.
"""

import argparse
import logging
import textwrap
import time
from pathlib import Path

import numpy as np
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

OUT_DIR = Path("data/petersen_2026")
REF_DIR = OUT_DIR / "refs"
WGT_DIR = OUT_DIR / "weights"
MUT_DIR = OUT_DIR / "mutants"

TARGETS_YAML = Path("jobs/framework_B/sim_targets.yaml")

BASES = list("ACGT")
ALT   = {b: [x for x in BASES if x != b] for b in BASES}

DEFAULT_DEVIATIONS = [0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def load_targets_config() -> tuple[dict, list, int]:
    """Load sim_targets.yaml; return (by_name dict, deviations, seed)."""
    with open(TARGETS_YAML, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    by_name = {}
    for t in cfg["targets"]:
        by_name[t["name"]] = {
            "ref_acc": t["ref_acc"],
            "vimop_name": t["vimop_name"],
        }
    deviations = cfg.get("deviations", DEFAULT_DEVIATIONS)
    seed       = cfg.get("seed", 42)
    return by_name, deviations, seed


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_reference(_name: str, ref_acc: str, log: logging.Logger) -> str:
    """Return the ungapped reference sequence as a plain uppercase string."""
    fa = REF_DIR / f"{ref_acc}.fa"
    if not fa.exists():
        raise FileNotFoundError(f"Reference not found: {fa}  (run sim_01 first)")
    rec = next(SeqIO.parse(fa, "fasta"))
    seq = str(rec.seq).upper().replace("-", "").replace("N", "")
    log.info("  reference %s: %d nt", ref_acc, len(seq))
    return seq


def load_weights(name: str, log: logging.Logger) -> np.ndarray:
    """Load per-position variability weights (numpy array, one value per ref nt)."""
    npy = WGT_DIR / f"{name}.weights.npy"
    if not npy.exists():
        raise FileNotFoundError(f"Weights not found: {npy}  (run sim_02 first)")
    w = np.load(str(npy)).astype(np.float64)
    log.info("  weights: %d positions, max=%.0f, variable=%d",
             len(w), w.max(), int((w > 0).sum()))
    return w


def mutate(seq: str, deviation: float, weights: np.ndarray,
           rng: np.random.Generator,
           log: logging.Logger) -> tuple[str, list[int]]:
    """
    Return (mutated_seq, list_of_mutated_positions).

    n_sites = ceil(len_ref × deviation)  — matches ViMOP supplementary intent.
    Positions sampled without replacement weighted by variability weights.
    """
    n_sites = int(np.ceil(len(seq) * deviation))
    if n_sites == 0:
        return seq, []

    w = weights[:len(seq)].copy()
    if w.sum() == 0:
        log.warning("  all weights zero — falling back to uniform sampling")
        w = np.ones(len(seq), dtype=np.float64)
    probs = w / w.sum()

    chosen = rng.choice(len(seq), size=n_sites, replace=False, p=probs)
    chosen_sorted = sorted(chosen.tolist())

    seq_list = list(seq)
    for pos in chosen_sorted:
        orig = seq_list[pos]
        alts = ALT.get(orig, [b for b in BASES if b != orig])
        if alts:
            seq_list[pos] = rng.choice(alts)

    return "".join(seq_list), chosen_sorted


def write_fasta(seq: str, ref_acc: str, vimop_name: str,
                deviation: float, out_path: Path,
                n_mutations: int) -> None:
    dev_str = f"{int(round(deviation * 100)):02d}"
    desc = (f"artificial genome derived from {ref_acc}, "
            f"deviation={deviation:.2f}, n_mutations={n_mutations}")
    record = SeqRecord(Seq(seq),
                       id=f"{vimop_name}_{dev_str}",
                       description=desc)
    SeqIO.write([record], out_path, "fasta")


# ---------------------------------------------------------------------------
# Per-target processing
# ---------------------------------------------------------------------------

def process_target(name: str, cfg: dict, deviations: list,
                   rng: np.random.Generator,
                   log: logging.Logger) -> list[dict]:
    log.info("=== %s ===", name.upper())

    ref_acc    = cfg["ref_acc"]
    vimop_name = cfg["vimop_name"]

    ref_seq = load_reference(name, ref_acc, log)
    weights = load_weights(name, log)

    if len(weights) != len(ref_seq):
        log.warning(
            "  weight length (%d) != ref length (%d); truncating/padding",
            len(weights), len(ref_seq)
        )
        if len(weights) < len(ref_seq):
            weights = np.pad(weights, (0, len(ref_seq) - len(weights)))
        else:
            weights = weights[:len(ref_seq)]

    target_dir = MUT_DIR / vimop_name
    target_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for deviation in deviations:
        dev_str = f"{int(round(deviation * 100)):02d}"   # 0.05 → "05"
        out_fa  = target_dir / f"{vimop_name}_{dev_str}.fasta"

        mut_seq, mut_pos = mutate(ref_seq, deviation, weights, rng, log)
        write_fasta(mut_seq, ref_acc, vimop_name, deviation, out_fa,
                    len(mut_pos))

        n_requested = int(np.ceil(len(ref_seq) * deviation))
        actual      = sum(1 for a, b in zip(ref_seq, mut_seq) if a != b)
        log.info("  dev=%.2f → requested %d, actual %d → %s",
                 deviation, n_requested, actual, out_fa.name)

        rows.append({
            "target":      name,
            "vimop_name":  vimop_name,
            "deviation":   f"{deviation:.2f}",
            "ref_len":     len(ref_seq),
            "n_requested": n_requested,
            "n_actual":    actual,
            "output":      str(out_fa),
        })
    return rows


def write_summary(rows: list[dict], log: logging.Logger) -> None:
    tsv = MUT_DIR / "summary.tsv"
    header = ["target", "vimop_name", "deviation", "ref_len",
              "n_requested", "n_actual", "output"]
    with open(tsv, "w", encoding="utf-8") as fout:
        fout.write("\t".join(header) + "\n")
        for row in rows:
            fout.write("\t".join(str(row[h]) for h in header) + "\n")
    log.info("summary -> %s", tsv)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Main entry point"""
    ap = argparse.ArgumentParser(
        description=textwrap.dedent(__doc__ or ""),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    try:
        targets_cfg, default_devs, default_seed = load_targets_config()
        target_choices = list(targets_cfg.keys())
    except Exception:
        targets_cfg, default_devs, default_seed = {}, DEFAULT_DEVIATIONS, 42
        target_choices = ["COVID", "DENV", "LASV_L", "LASV_S"]

    ap.add_argument("--targets", nargs="+", choices=target_choices,
                    default=target_choices,
                    help="Targets to simulate (default: all).")
    ap.add_argument("--deviations", nargs="+", type=float, default=None,
                    help="Deviation fractions to generate "
                         "(default: from sim_targets.yaml).")
    ap.add_argument("--seed", type=int, default=None,
                    help="Random seed (default: from sim_targets.yaml).")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    log_dir  = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"sim_03_{time.strftime('%Y%m%d_%H%M%S')}.log"

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

    log = logging.getLogger("sim03")
    log.info("log -> %s", log_file)

    targets_cfg, default_devs, default_seed = load_targets_config()
    deviations = args.deviations if args.deviations is not None else default_devs
    seed       = args.seed if args.seed is not None else default_seed
    log.info("deviations: %s", deviations)
    log.info("seed: %d", seed)

    MUT_DIR.mkdir(parents=True, exist_ok=True)

    rng  = np.random.default_rng(seed)
    rows = []
    for name in args.targets:
        rows.extend(process_target(name, targets_cfg[name], deviations, rng, log))

    write_summary(rows, log)
    log.info("done — %d artificial genomes in %s", len(rows), MUT_DIR)


if __name__ == "__main__":
    main()
