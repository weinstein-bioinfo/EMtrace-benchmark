#!/usr/bin/env python3
"""
D_sim00_taxonomy.py — Framework D, step 0: derive per-source taxonomy labels.

For every simulated-read source used in the combined Nanopore datasets — the four
viral targets in sim_targets.yaml and the nine background organisms in
assets/background_refs.yaml — resolve a ground-truth taxonomy label from its
accession via NCBI Entrez:

    source   accession   taxid   species   division

  * division ∈ {Viruses, Bacteria, Eukaryota} (the superkingdom rank).
  * taxid / species are the species-rank ancestor (clean species-level truth),
    with the organism (strain) scientific name kept for reference.

The result is written to assets/source_taxonomy.tsv and printed for verification.
This file is consumed by D_sim06_combine_labeled.py to label each read. It is kept
separate from the two source YAMLs so their curated comments are not clobbered.

Usage:
  python3 jobs/framework_D/D_sim00_taxonomy.py [-v]

Dependencies: python3 + biopython (Bio.Entrez) + pyyaml. Needs network access.
"""

import argparse
import logging
import sys
import time
from pathlib import Path

import yaml
from Bio import Entrez

ENTREZ_EMAIL = "your.email@example.org"

TARGETS_YAML = Path("jobs/framework_B/sim_targets.yaml")
BG_YAML = Path("jobs/assets/background_refs.yaml")
OUT_TSV = Path("jobs/assets/source_taxonomy.tsv")

#: superkingdom ScientificName → division label used in the ground-truth table.
DIVISION_FROM_SUPERKINGDOM = {
    "Viruses": "Viruses",
    "Bacteria": "Bacteria",
    "Archaea": "Archaea",
    "Eukaryota": "Eukaryota",
}


def nuccore_taxid(accession: str, log: logging.Logger) -> int:
    """Resolve a nuccore accession (e.g. OX637002.1) to its NCBI TaxId."""
    handle = Entrez.esummary(db="nuccore", id=accession)
    rec = Entrez.read(handle)
    handle.close()
    taxid = int(rec[0]["TaxId"])
    log.debug("  nuccore %s → taxid %d", accession, taxid)
    return taxid


def assembly_taxid(gcf_accession: str, log: logging.Logger) -> int:
    """Resolve a GCF assembly accession to its organism TaxId."""
    handle = Entrez.esearch(db="assembly",
                            term=f"{gcf_accession}[Assembly Accession]")
    rec = Entrez.read(handle)
    handle.close()
    ids = rec.get("IdList", [])
    if not ids:
        raise RuntimeError(f"No assembly record for {gcf_accession}")
    handle = Entrez.esummary(db="assembly", id=ids[0], report="full")
    summary = Entrez.read(handle)
    handle.close()
    doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
    taxid = int(doc["Taxid"])
    log.debug("  assembly %s → taxid %d", gcf_accession, taxid)
    return taxid


def resolve_taxonomy(taxid: int, log: logging.Logger) -> dict:
    """efetch the taxonomy record → species (name+taxid), division, organism.

    division is the superkingdom ScientificName; species is the species-rank
    ancestor (or the organism itself when it is already a species).
    """
    handle = Entrez.efetch(db="taxonomy", id=str(taxid))
    rec = Entrez.read(handle)
    handle.close()
    node = rec[0]
    organism = node["ScientificName"]
    lineage = list(node.get("LineageEx", [])) + [
        {"TaxId": node["TaxId"], "ScientificName": organism, "Rank": node.get("Rank", "")}
    ]

    species_name, species_taxid = organism, int(taxid)
    for entry in lineage:
        if entry.get("Rank", "") == "species":
            species_name = entry["ScientificName"]
            species_taxid = int(entry["TaxId"])

    # Division by ScientificName membership in the lineage — robust to NCBI's
    # superkingdom→domain rank rename and to Viruses sitting under a realm/root.
    names = {e["ScientificName"] for e in lineage}
    division = next((DIVISION_FROM_SUPERKINGDOM[d] for d in
                     ("Viruses", "Bacteria", "Archaea", "Eukaryota") if d in names), "")
    if not division:
        log.warning("  taxid %d: no domain/superkingdom in lineage; division blank", taxid)
    return {
        "taxid": species_taxid,
        "species": species_name,
        "division": division,
        "organism": organism,
    }


def collect_sources(log: logging.Logger) -> list[dict]:
    """Build the (source, accession, kind) work list from the two YAMLs."""
    with open(TARGETS_YAML, encoding="utf-8") as fh:
        targets = yaml.safe_load(fh)["targets"]
    with open(BG_YAML, encoding="utf-8") as fh:
        bg = yaml.safe_load(fh)["sequences"]

    rows: list[dict] = []
    for t in targets:
        rows.append({"source": t["name"], "accession": t["ref_acc"], "kind": "nuccore"})
    for s in bg:
        rows.append({"source": s["name"], "accession": s["accession"], "kind": "assembly"})
    log.info("sources: %d viral + %d background = %d",
             len(targets), len(bg), len(rows))
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(levelname)s] %(message)s")
    log = logging.getLogger("D_sim00")

    Entrez.email = ENTREZ_EMAIL
    for yml in (TARGETS_YAML, BG_YAML):
        if not yml.exists():
            log.error("missing config: %s", yml)
            sys.exit(1)

    rows = collect_sources(log)
    out: list[dict] = []
    for r in rows:
        try:
            taxid = (nuccore_taxid(r["accession"], log) if r["kind"] == "nuccore"
                     else assembly_taxid(r["accession"], log))
            tax = resolve_taxonomy(taxid, log)
        except Exception as exc:  # noqa: BLE001 — report and continue
            log.error("  %s (%s): %s", r["source"], r["accession"], exc)
            continue
        out.append({"source": r["source"], "accession": r["accession"], **tax})
        time.sleep(0.34)  # stay under NCBI's 3 req/s without an API key

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    cols = ["source", "accession", "taxid", "species", "division", "organism"]
    with open(OUT_TSV, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in out:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")

    # Verification print
    log.info("=== derived taxonomy (verify before running D_sim06) ===")
    width = max((len(r["source"]) for r in out), default=6)
    for r in out:
        log.info("  %-*s  %-12s  taxid=%-8s  %-9s  %s",
                 width, r["source"], r["accession"], r["taxid"], r["division"], r["species"])
    log.info("wrote %d rows → %s", len(out), OUT_TSV)


if __name__ == "__main__":
    main()
