#!/usr/bin/env python3
"""
B_sim02_align_and_weights.py — Framework B genome simulation, step 2.

Two source modes (--source):

  vimop (default)
    Read the pre-aligned multifasta directly from ViMOP's curated
    data/virus_references/*.aligned.fasta.gz.  MAFFT is skipped.
    Alignment paths are resolved from sim_targets.yaml.

  local
    Read the multifasta produced by sim_01 and align with MAFFT
    (--addfragments --keeplength against the downloaded reference).
    Fallback for viruses not covered by ViMOP.

In both modes the weight computation is identical:
  weight[col] = number of sequences that differ from the per-column
                consensus at that column (gaps ignored).
  Highly variable sites are preferentially chosen for mutagenesis.

Outputs (per target):
  aligned/<target>.aln.fa           — alignment written only in 'local' mode
  weights/<target>.weights.tsv      — reference_pos, ref_base, weight
  weights/<target>.weights.npy      — numpy array indexed by reference position

Dependencies: mafft (local mode only), python3 + biopython + numpy + pyyaml.
"""

import argparse
import gzip
import logging
import subprocess
import textwrap
import time
from collections import Counter
from pathlib import Path

import numpy as np
import yaml
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

OUT_DIR = Path("data/petersen_2026")
REF_DIR = OUT_DIR / "refs"
MFA_DIR = OUT_DIR / "multifastas"
ALN_DIR = OUT_DIR / "aligned"
WGT_DIR = OUT_DIR / "weights"

TARGETS_YAML = Path("jobs/framework_B/sim_targets.yaml")

TARGETS = ["COVID", "DENV", "LASV_L", "LASV_S"]

REF_ACCS = {
    "COVID":  "OX637002.1",
    "DENV":   "NC_001474.2",
    "LASV_L": "MG812630.1",
    "LASV_S": "MG812631.1",
}


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------

def load_targets_config() -> dict:
    """Load sim_targets.yaml; return dict keyed by target name."""
    with open(TARGETS_YAML, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    vimop_root = Path(cfg["vimop_root"])
    by_name = {}
    for t in cfg["targets"]:
        by_name[t["name"]] = {
            "ref_acc":       t["ref_acc"],
            "vimop_name":    t["vimop_name"],
            "aligned_gz":    vimop_root / t["vimop_aligned"],
        }
    return by_name


# ---------------------------------------------------------------------------
# Alignment helpers
# ---------------------------------------------------------------------------

def run(cmd: list[str], log: logging.Logger, **kw) -> subprocess.CompletedProcess:
    log.debug("  $ %s", " ".join(str(c) for c in cmd))
    return subprocess.run(cmd, check=True, **kw)


def load_vimop_alignment(gz_path: Path, log: logging.Logger) -> MultipleSeqAlignment:
    """Read a gzip-compressed pre-aligned FASTA directly."""
    log.info("  loading ViMOP alignment from %s ...", gz_path.name)
    with gzip.open(gz_path, "rt") as fh:
        aln = AlignIO.read(fh, "fasta")
    log.info("  loaded: %d sequences × %d columns",
             len(aln), aln.get_alignment_length())
    return aln


def align_mafft(fa_in: Path, ref_fa: Path, aln_out: Path, threads: int,
                log: logging.Logger) -> Path:
    """
    Align fa_in against ref_fa using --addfragments --keeplength.
    The reference defines the alignment columns; insertions in query sequences
    are discarded so every alignment column maps 1-to-1 to a reference position.
    """
    if aln_out.exists():
        log.info("  alignment already exists — skipping MAFFT")
        return aln_out
    n = sum(1 for _ in SeqIO.parse(fa_in, "fasta"))
    log.info("  MAFFT: adding %d sequences to reference %s (--keeplength) ...",
             n, ref_fa.name)
    with open(aln_out, "w", encoding="utf-8") as fout:
        run([
            "mafft",
            "--addfragments", str(fa_in),
            "--keeplength",
            "--thread", str(threads),
            "--quiet",
            str(ref_fa),
        ], log, stdout=fout)
    log.info("  alignment saved to %s", aln_out.name)
    return aln_out


# ---------------------------------------------------------------------------
# Weight computation
# ---------------------------------------------------------------------------

def consensus_and_weights(aln: MultipleSeqAlignment,
                          log: logging.Logger) -> tuple[str, np.ndarray]:
    """
    Return (consensus_str, weights_array) over alignment columns.

    weights_array[col] = count of sequences differing from the column consensus
                         (gaps and N excluded from both numerator and denominator).
    """
    n_seqs = len(aln)
    n_cols = aln.get_alignment_length()
    log.info("  alignment: %d sequences × %d columns", n_seqs, n_cols)

    consensus = []
    weights = np.zeros(n_cols, dtype=np.int32)

    for col in range(n_cols):
        bases = [str(aln[i, col]).upper() for i in range(n_seqs)]
        non_gap = [b for b in bases if b not in ("-", "N", "?")]
        if not non_gap:
            consensus.append("N")
            weights[col] = 0
            continue
        cons_base = Counter(non_gap).most_common(1)[0][0]
        consensus.append(cons_base)
        weights[col] = sum(1 for b in non_gap if b != cons_base)

    return "".join(consensus), weights


def map_to_reference(aln: MultipleSeqAlignment, ref_acc: str,
                     weights: np.ndarray,
                     log: logging.Logger) -> tuple[np.ndarray, str]:
    """
    Find the reference row in the alignment; return (ref_weights, ref_seq).

    In vimop mode the alignment is already reference-length so the reference
    row has no gaps.  In local (--keeplength) mode the same holds.
    Any residual gaps are stripped and weights truncated to match.
    """
    ref_row = None
    for rec in aln:
        if rec.id.split()[0] == ref_acc or rec.id == ref_acc:
            ref_row = rec
            break
    if ref_row is None:
        raise ValueError(
            f"Reference accession {ref_acc!r} not found in alignment. "
            f"IDs present: {[r.id for r in aln][:5]}..."
        )

    ref_seq = str(ref_row.seq).upper().replace("-", "")
    if len(ref_seq) != len(weights):
        log.warning(
            "  reference row length (%d) != weight vector length (%d) "
            "— gaps present in reference; weights truncated to ref length",
            len(ref_seq), len(weights),
        )
    log.info("  reference ungapped length: %d nt", len(ref_seq))
    return weights[:len(ref_seq)].astype(np.float64), ref_seq


# ---------------------------------------------------------------------------
# Per-target processing
# ---------------------------------------------------------------------------

def process_target(name: str, source: str, targets_cfg: dict,
                   threads: int, log: logging.Logger) -> None:
    log.info("=== %s (source=%s) ===", name.upper(), source)

    tsv_out = WGT_DIR / f"{name}.weights.tsv"
    npy_out = WGT_DIR / f"{name}.weights.npy"
    ref_acc = REF_ACCS[name]

    if source == "vimop":
        cfg      = targets_cfg[name]
        gz_path  = cfg["aligned_gz"]
        if not gz_path.exists():
            log.error("  ViMOP alignment not found: %s", gz_path)
            return
        aln = load_vimop_alignment(gz_path, log)

    else:  # local
        fa_in   = MFA_DIR / f"{name}.fa"
        aln_out = ALN_DIR / f"{name}.aln.fa"
        ref_fa  = REF_DIR / f"{ref_acc}.fa"
        if not fa_in.exists():
            log.error("  multifasta not found: %s — run sim_01 first", fa_in)
            return
        if not ref_fa.exists():
            log.error("  reference FASTA not found: %s — run sim_01 first", ref_fa)
            return
        align_mafft(fa_in, ref_fa, aln_out, threads, log)
        aln = AlignIO.read(str(aln_out), "fasta")

    # Consensus + column weights
    _, col_weights = consensus_and_weights(aln, log)

    # Map to reference positions
    ref_weights, ref_seq = map_to_reference(aln, ref_acc, col_weights, log)

    # Save weights TSV
    with open(tsv_out, "w", encoding='utf-8') as fout:
        fout.write("ref_pos\tref_base\tweight\n")
        for pos, (base, w) in enumerate(zip(ref_seq, ref_weights)):
            fout.write(f"{pos}\t{base}\t{w:.0f}\n")
    log.info("  weights saved to %s", tsv_out.name)

    # Save numpy array
    np.save(str(npy_out), ref_weights)
    log.info("  numpy weights saved to %s", npy_out.name)

    nonzero = int((ref_weights > 0).sum())
    log.info("  variable positions: %d / %d (%.1f%%)",
             nonzero, len(ref_weights),
             100 * nonzero / len(ref_weights) if ref_weights.size else 0)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """Main entry point"""
    ap = argparse.ArgumentParser(
        description=textwrap.dedent(__doc__ or ""),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--targets", nargs="+", choices=TARGETS,
                    default=TARGETS, help="Targets to process (default: all).")
    ap.add_argument("--source", choices=["vimop", "local"], default="vimop",
                    help="Alignment source: 'vimop' reads ViMOP's pre-aligned "
                         "FASTAs (default); 'local' runs MAFFT on sim_01 output.")
    ap.add_argument("--threads", type=int, default=8,
                    help="MAFFT threads (local mode only).")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    log_dir  = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"sim_02_{time.strftime('%Y%m%d_%H%M%S')}.log"

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

    log = logging.getLogger("sim02")
    log.info("log -> %s", log_file)
    log.info("source = %s", args.source)

    ALN_DIR.mkdir(parents=True, exist_ok=True)
    WGT_DIR.mkdir(parents=True, exist_ok=True)

    targets_cfg = load_targets_config() if args.source == "vimop" else {}

    for name in args.targets:
        process_target(name, args.source, targets_cfg, args.threads, log)

    log.info("done — weights in %s", WGT_DIR)


if __name__ == "__main__":
    main()
