#!/usr/bin/env python3
"""
sim_bg_download.py — Download background reference genomes for sim_05.

For each entry in jobs/assets/background_refs.yaml:
  1. Query NCBI Assembly via Entrez esummary to resolve the GCF accession to its
     RefSeq FTP path.
  2. Download the genomic FASTA (.fna.gz) to the configured output_dir.
  3. Write the local path back into background_refs.yaml under each entry's
     `fasta` field so that B_sim05_mix_background.py can find the files.

Note: Homo sapiens (GRCh38.p14, ~3 GB) is downloaded via the FTP URL found in
the Assembly record; this may take several minutes on a slow connection.

Usage:
  python3 jobs/framework_B/B_sim03_bg_download.py [-v]

Dependencies:
  python3 + biopython + pyyaml
"""

import argparse
import logging
import sys
import textwrap
import time
import urllib.request
from pathlib import Path

import yaml
from Bio import Entrez

ENTREZ_EMAIL = "your.email@example.org"
BG_YAML = Path("jobs/assets/background_refs.yaml")


# ---------------------------------------------------------------------------
# Entrez helpers
# ---------------------------------------------------------------------------

def resolve_ftp_path(gcf_accession: str, log: logging.Logger) -> str:
    """
    Return the RefSeq FTP directory URL for a GCF assembly accession.
    Raises RuntimeError if not found.
    """
    log.debug("  esearch assembly: %s", gcf_accession)
    handle = Entrez.esearch(db="assembly",
                            term=f"{gcf_accession}[Assembly Accession]")
    record = Entrez.read(handle)
    handle.close()
    ids = record.get("IdList", [])
    if not ids:
        raise RuntimeError(f"No assembly record found for {gcf_accession}")
    uid = ids[0]

    log.debug("  esummary uid: %s", uid)
    handle = Entrez.esummary(db="assembly", id=uid, report="full")
    summary = Entrez.read(handle)
    handle.close()
    doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
    ftp = doc.get("FtpPath_RefSeq", "")
    if not ftp or ftp == "na":
        raise RuntimeError(
            f"No RefSeq FTP path for {gcf_accession} (uid={uid}). "
            "Assembly may not have a RefSeq release."
        )
    return ftp          # e.g. https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/…


def download_fna(ftp_dir: str, dest_dir: Path,
                 log: logging.Logger) -> Path:
    """
    Download <assembly_name>_genomic.fna.gz from ftp_dir into dest_dir.
    Returns the local Path.  Skips if already present.
    """
    basename = ftp_dir.rstrip("/").split("/")[-1]
    filename = f"{basename}_genomic.fna.gz"
    url = f"{ftp_dir}/{filename}"
    dest = dest_dir / filename
    if dest.exists():
        log.info("  already on disk: %s", dest.name)
        return dest
    log.info("  downloading %s …", filename)
    log.debug("    url: %s", url)

    def _progress(block_num, block_size, total_size):
        if total_size > 0 and block_num % 500 == 0:
            pct = min(100, block_num * block_size * 100 // total_size)
            log.debug("    … %d%%", pct)

    urllib.request.urlretrieve(url, dest, reporthook=_progress)
    log.info("  saved: %s (%.1f MB)", dest.name, dest.stat().st_size / 1e6)
    return dest


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description=textwrap.dedent(__doc__ or ""),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    log_dir = Path("logs")
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"sim_bg_download_{time.strftime('%Y%m%d_%H%M%S')}.log"

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

    log = logging.getLogger("bg_download")
    log.info("log -> %s", log_file)

    if not BG_YAML.exists():
        log.error("background_refs.yaml not found: %s", BG_YAML)
        sys.exit(1)

    with open(BG_YAML, encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)

    out_dir = Path(cfg["output_dir"])
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("output dir: %s", out_dir)

    Entrez.email = ENTREZ_EMAIL
    sequences = cfg["sequences"]
    changed = False

    for entry in sequences:
        name = entry["name"]
        acc = entry["accession"]
        log.info("=== %s (%s) ===", name, acc)

        if entry.get("fasta") and Path(entry["fasta"]).exists():
            log.info("  fasta already present: %s", entry["fasta"])
            continue

        try:
            ftp_dir = resolve_ftp_path(acc, log)
            local = download_fna(ftp_dir, out_dir, log)
            entry["fasta"] = str(local)
            changed = True
            time.sleep(0.4)   # be polite to NCBI
        except Exception as exc:
            log.error("  FAILED %s: %s", acc, exc)

    if changed:
        with open(BG_YAML, "w", encoding="utf-8") as fh:
            yaml.dump(cfg, fh, allow_unicode=True, sort_keys=False)
        log.info("updated %s with local fasta paths", BG_YAML)
    else:
        log.info("nothing changed — all fastas already present")

    n_ok = sum(1 for e in sequences if e.get("fasta") and Path(e["fasta"]).exists())
    log.info("done: %d / %d entries have local fasta", n_ok, len(sequences))
    if n_ok < len(sequences):
        log.warning("%d entries still missing fasta — check errors above",
                    len(sequences) - n_ok)
        sys.exit(1)


if __name__ == "__main__":
    main()
