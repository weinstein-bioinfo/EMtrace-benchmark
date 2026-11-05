#!/usr/bin/env python3
"""
B_sim01_build_multifastas.py — Framework B genome simulation, step 1.

For each of 4 target viruses (SARS-CoV-2, DENV, Lassa-L, Lassa-S):
  1. Fetch the reference FASTA from NCBI Entrez.
  2. Extract matching sequences from the local GenBank viral database
     (minimap.genbank.vi.fa.gz) using assembly_summary.txt species_taxid
     → leaf taxids → acc2taxid → seqkit grep.
  3. For Lassa: separate L (>5 kb) from S (≤5 kb) segments by length.
  4. Concatenate reference + pool into a raw multifasta.
  5. Deduplicate by exact sequence (seqkit rmdup -s).
  6. SARS-CoV-2 only (optional): cluster at 99.9% identity with MMseqs2
     easy-cluster to reduce 12 k sequences to a manageable alignment input.
     Disabled by default — set "cluster": True in TARGETS to activate.
  7. Write final multifastas to data/petersen_2026/multifastas/.

Outputs (per target):
  multifastas/<target>.fa           — final deduplicated / clustered multifasta
  multifastas/<target>.raw.fa       — pre-dedup concatenation (kept for audit)

Dependencies: python3 + biopython, seqkit, mmseqs2 (SARS-CoV-2 clustering, optional).
  brew install mmseqs2               # or: conda install -c bioconda mmseqs2
"""

import argparse
import logging
import os
import shutil
import subprocess
import textwrap
import time
from pathlib import Path

from Bio import Entrez, SeqIO

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

ENTREZ_EMAIL = "your.email@example.org"

DB_DIR = Path(os.environ.get("GENBANK_VI_DIR", "databases/genbank.vi.20260507"))  # GenBank viral DB
MERGED_FA          = DB_DIR / "minimap.genbank.vi.fa.gz"
ACC2TAXID          = DB_DIR / "minimap.genbank.vi.fa.acc2taxid"
ASSEMBLY_SUMMARY   = DB_DIR / "genbank.vi.assembly_summary.txt"

OUT_DIR = Path("data/petersen_2026")
REF_DIR = OUT_DIR / "refs"
MFA_DIR = OUT_DIR / "multifastas"

TARGETS = {
    "COVID": {
        "ref_acc":         "OX637002.1",
        "species_taxid":   3418604,
        # 2697049 is a strain-level taxon; filter to high-quality assemblies only
        "assembly_level_filter": {"Complete Genome", "Chromosome"},
        "segment":         None,
        "cluster":         False,  # set True to enable MMseqs2 clustering at 99.9 %
    },
    "DENV": {
        "ref_acc":         "NC_001474.2",
        "species_taxid":   3052464,   # covers DENV-1/2/3/4
        "assembly_level_filter": None,
        "segment":         None,
        "cluster":         False,
    },
    "LASV_L": {
        "ref_acc":         "MG812630.1",
        "species_taxid":   3052310,   # Mammarenavirus lassaense
        "assembly_level_filter": None,
        "segment":         "L",    # length > 5000 nt
        "cluster":         False,
    },
    "LASV_S": {
        "ref_acc":         "MG812631.1",
        "species_taxid":   3052310,
        "assembly_level_filter": None,
        "segment":         "S",    # length ≤ 5000 nt
        "cluster":         False,
    },
}

LASSA_L_MIN_LEN = 5000   # threshold separating L and S segments


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def run(cmd: list[str], log: logging.Logger, **kw) -> subprocess.CompletedProcess:
    log.debug("  $ %s", " ".join(str(c) for c in cmd))
    result = subprocess.run(cmd, check=True, **kw)
    return result


def fetch_reference(acc: str, dest: Path, log: logging.Logger) -> Path:
    """Download a single accession from NCBI nuccore in FASTA format."""
    out = dest / f"{acc}.fa"
    if out.exists():
        log.info("  ref %s already on disk — skipping fetch", acc)
        return out
    log.info("  fetching %s from NCBI ...", acc)
    Entrez.email = ENTREZ_EMAIL
    for attempt in range(3):
        try:
            handle = Entrez.efetch(db="nuccore", id=acc,
                                   rettype="fasta", retmode="text")
            seq = handle.read()
            handle.close()
            break
        except Exception as exc:
            log.warning("  attempt %d failed: %s", attempt + 1, exc)
            time.sleep(5)
    else:
        raise RuntimeError(f"Failed to fetch {acc} after 3 attempts")
    out.write_text(seq, encoding="utf-8")
    log.info("  saved reference to %s", out)
    return out


def extract_taxid_accessions(cfg: dict, log: logging.Logger) -> list[str]:
    """
    Return all sequence accessions for the target by:
      1. Parsing assembly_summary.txt, filtering by species_taxid
         (and assembly_level if cfg['assembly_level_filter'] is set).
      2. Collecting leaf taxids from matching assemblies.
      3. Matching those leaf taxids against acc2taxid to get sequence accessions.
    """
    species_taxid  = cfg["species_taxid"]
    level_filter   = cfg.get("assembly_level_filter")

    # --- step 1 & 2: assembly_summary → leaf taxids ---
    leaf_taxids: set[int] = set()
    n_assemblies = 0
    with open(ASSEMBLY_SUMMARY, encoding="utf-8") as fh:
        raw_header = fh.readline().lstrip("# ").rstrip("\n")
        header = raw_header.split("\t")
        try:
            idx_sp  = header.index("species_taxid")
            idx_tax = header.index("taxid")
            idx_lvl = header.index("assembly_level")
        except ValueError as exc:
            raise RuntimeError(
                f"assembly_summary.txt missing expected column: {exc}\n"
                f"Found columns: {header[:10]}"
            )
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_sp, idx_tax, idx_lvl):
                continue
            if int(parts[idx_sp]) != species_taxid:
                continue
            if level_filter and parts[idx_lvl] not in level_filter:
                continue
            leaf_taxids.add(int(parts[idx_tax]))
            n_assemblies += 1

    log.info("  %d assemblies matched in assembly_summary (species_taxid=%d%s)",
             n_assemblies, species_taxid,
             f", level_filter={sorted(level_filter)}" if level_filter else "")

    # --- step 3: acc2taxid → sequence accessions ---
    matching: list[str] = []
    with open(ACC2TAXID, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if int(parts[1]) in leaf_taxids:
                matching.append(parts[0])
    log.info("  %d sequence accessions found for %d leaf taxids",
             len(matching), len(leaf_taxids))
    return matching


def seqkit_extract(accessions: list[str], out_fa: Path,
                   log: logging.Logger) -> Path:
    """Extract sequences from the merged FASTA by accession list using seqkit."""
    acc_file = out_fa.with_suffix(".accs.txt")
    acc_file.write_text("\n".join(accessions) + "\n", encoding="utf-8")
    run(["seqkit", "grep", "-f", str(acc_file),
         str(MERGED_FA), "-o", str(out_fa)], log)
    n = count_seqs(out_fa)
    log.info("  extracted %d sequences -> %s", n, out_fa.name)
    return out_fa


def filter_by_length(fa_in: Path, fa_out: Path, min_len: int,
                     mode: str, log: logging.Logger) -> Path:
    """Keep sequences where length > min_len (mode='L') or <= min_len (mode='S')."""
    op = ">" if mode == "L" else "<="
    records = [r for r in SeqIO.parse(fa_in, "fasta")
               if (len(r.seq) > min_len) == (mode == "L")]
    SeqIO.write(records, fa_out, "fasta")
    log.info("  segment filter (len %s %d): %d sequences kept -> %s",
             op, min_len, len(records), fa_out.name)
    return fa_out


def count_seqs(fa: Path) -> int:
    return sum(1 for _ in SeqIO.parse(fa, "fasta"))


def concatenate(fastas: list[Path], out: Path) -> None:
    with open(out, "w", encoding="utf-8") as fout:
        for fa in fastas:
            fout.write(fa.read_text(encoding="utf-8"))


def dedup_exact(fa_in: Path, fa_out: Path, log: logging.Logger) -> Path:
    """Remove exact-sequence duplicates with seqkit rmdup."""
    run(["seqkit", "rmdup", "-s", "-i", "-o", str(fa_out), str(fa_in)], log)
    before = count_seqs(fa_in)
    after  = count_seqs(fa_out)
    log.info("  exact dedup: %d -> %d sequences", before, after)
    return fa_out


def cluster_mmseqs(fa_in: Path, fa_out: Path, identity: float,
                   coverage: float, threads: int,
                   log: logging.Logger) -> Path:
    """
    Cluster sequences with MMseqs2 easy-cluster and return the representative FASTA.

    Parameters mirror VOTUpipe defaults:
      --min-seq-id  identity   (0.999 for 99.9 % ANI)
      -c            coverage   (0.85, shorter-sequence mode --cov-mode 1)
      --cov-mode    1          coverage of the shorter sequence
    """
    if not shutil.which("mmseqs"):
        raise RuntimeError(
            "mmseqs not found. Install with: brew install mmseqs2  "
            "or: conda install -c bioconda mmseqs2"
        )
    tmp_dir    = fa_out.parent / f"_{fa_out.stem}_mmseqs_tmp"
    out_prefix = fa_out.parent / f"_{fa_out.stem}_mmseqs"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    run([
        "mmseqs", "easy-cluster",
        str(fa_in),
        str(out_prefix),
        str(tmp_dir),
        "--min-seq-id", str(identity),
        "-c",           str(coverage),
        "--cov-mode",   "1",
        "--threads",    str(threads),
    ], log, capture_output=True)

    rep_fasta = Path(f"{out_prefix}_rep_seq.fasta")
    if not rep_fasta.exists():
        raise RuntimeError(f"MMseqs2 representative FASTA not found: {rep_fasta}")
    rep_fasta.rename(fa_out)

    shutil.rmtree(tmp_dir, ignore_errors=True)
    for leftover in fa_out.parent.glob(f"{out_prefix.name}_*.fasta"):
        leftover.unlink(missing_ok=True)
    for leftover in fa_out.parent.glob(f"{out_prefix.name}_*.tsv"):
        leftover.unlink(missing_ok=True)

    before = count_seqs(fa_in)
    after  = count_seqs(fa_out)
    log.info("  mmseqs2 easy-cluster (%.1f%% id, %.0f%% cov): %d -> %d sequences",
             identity * 100, coverage * 100, before, after)
    return fa_out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_target(name: str, cfg: dict, threads: int,
                 log: logging.Logger) -> Path:
    log.info("=== %s ===", name.upper())
    tmp_dir = MFA_DIR / f"_{name}_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # 1. Reference
    ref_fa = fetch_reference(cfg["ref_acc"], REF_DIR, log)

    # 2. Pool from database
    accs = extract_taxid_accessions(cfg, log)
    pool_raw = tmp_dir / "pool_raw.fa"
    seqkit_extract(accs, pool_raw, log)

    # 3. Lassa segment filter
    pool_fa = pool_raw
    if cfg["segment"]:
        pool_filt = tmp_dir / f"pool_{cfg['segment']}.fa"
        filter_by_length(pool_raw, pool_filt, LASSA_L_MIN_LEN,
                         cfg["segment"], log)
        pool_fa = pool_filt

    # 4. Concatenate ref + pool
    raw_out = MFA_DIR / f"{name}.raw.fa"
    concatenate([ref_fa, pool_fa], raw_out)
    log.info("  raw multifasta: %d sequences -> %s",
             count_seqs(raw_out), raw_out.name)

    # 5. Exact dedup
    dedup_out = tmp_dir / "dedup.fa"
    dedup_exact(raw_out, dedup_out, log)

    # 6. Cluster (SARS-CoV-2 only)
    final_out = MFA_DIR / f"{name}.fa"
    if cfg["cluster"]:
        cluster_mmseqs(dedup_out, final_out, identity=0.999,
                       coverage=0.85, threads=threads, log=log)
    else:
        shutil.copy(dedup_out, final_out)
        log.info("  final multifasta: %d sequences -> %s",
                 count_seqs(final_out), final_out.name)

    # Cleanup tmp
    shutil.rmtree(tmp_dir, ignore_errors=True)
    return final_out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=textwrap.dedent(__doc__ or ""),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--targets", nargs="+", choices=list(TARGETS),
                    default=list(TARGETS),
                    help="Which targets to build (default: all).")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    log_dir  = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"sim_01_{time.strftime('%Y%m%d_%H%M%S')}.log"

    fmt = logging.Formatter("[%(asctime)s %(levelname)s] %(message)s",
                            datefmt="%H:%M:%S")
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

    log = logging.getLogger("sim01")
    log.info("log -> %s", log_file)

    REF_DIR.mkdir(parents=True, exist_ok=True)
    MFA_DIR.mkdir(parents=True, exist_ok=True)

    for name in args.targets:
        cfg = TARGETS[name]
        out = build_target(name, cfg, args.threads, log)
        log.info("  -> %s  (%d sequences)", out, count_seqs(out))

    log.info("done — multifastas in %s", MFA_DIR)


if __name__ == "__main__":
    main()
