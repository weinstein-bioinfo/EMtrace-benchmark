#!/bin/bash
# =============================================================================
# Framework A — MEGAN-LR comparator: ONE-TIME LAST protein index build (Roihu)
# =============================================================================
# LAST+MEGAN-LR pipeline (Huson et al. 2018, Fig. 4): align long reads against a
# protein reference DB with LAST, convert MAF -> DAA (maf2daa), meganize (daa-
# meganizer), then export taxonomic counts (daa2info). This script does step 0:
# build the `lastdb` protein index those alignments run against.
#
# PROTEIN DB CHOICE (updated 2026-08): reuses the SAME cross-pipeline-consistent
# protein FASTA that metaMix's DIAMOND DB is now built from — viral + archaea/
# bacteria 1/genus + human, RefSeq release 234, built by
# jobs/framework_A/db_refseq_abvh_protein.{dload,build}.sbatch (Lazypipe-db's
# db_refseq.py --target abvh --genus --seqtype prot).
#
# Consequence: this script REQUIRES
# jobs/framework_A/db_refseq_abvh_protein.build.sbatch to have completed
# first (same $databases/.../refseq.abvh.1g.234/diamond.refseq.abvh.1g.faa
# input — the combined FASTA sits directly in the release subdir.
#
#
# USAGE: sbatch framework_A/megan-lr/builddb.sbatch
# =============================================================================
#SBATCH --job-name=megan_lr_builddb
#SBATCH --account=Project_2003755
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G                           # curated RefSeq subset, not full nr — retune after first run (seff)
#SBATCH --time=12:00:00                     # rough — retime after first run; lastdb protein indexing of the curated FASTA is far smaller than a full-nr build
#SBATCH --output=logs/megan_lr_builddb_%j.out
#SBATCH --error=logs/megan_lr_builddb_%j.err

set -euo pipefail

PROJAPPL="/projappl/project_XXXXXXX"
CONDA_ENV="$PROJAPPL/venv/megan-lr"          # bioconda: last + megan (megan-lr_setup.md item 1)

: "${databases:?set \$databases (Roihu reference-DB root)}"

RELEASE="234"                               # MUST match db_refseq_abvh_protein.{dload,build}.sbatch's RELEASE
PROTEIN_FAA="$databases/2026/NCBI/refseq.abvh.1g.${RELEASE}/diamond.refseq.abvh.1g.faa"
LASTDB_DIR="$databases/megan_lr"
LASTDB_INDEX="$LASTDB_DIR/refseq.abvh.1g.${RELEASE}"

export PATH="$CONDA_ENV/bin:$PATH"

if [[ ! -s "$PROTEIN_FAA" ]]; then
  echo "[error] $PROTEIN_FAA not found — run" \
       "jobs/framework_A/db_refseq_abvh_protein.dload.sbatch then" \
       "db_refseq_abvh_protein.build.sbatch first" >&2
  exit 1
fi

mkdir -p "$LASTDB_DIR"

if [[ -s "${LASTDB_INDEX}.prj" ]]; then
  echo "[skip] ${LASTDB_INDEX}.prj already present"
  exit 0
fi

echo "[1/1] lastdb -p (protein index) <- $PROTEIN_FAA"
lastdb -p -c -P "$SLURM_CPUS_PER_TASK" "$LASTDB_INDEX" "$PROTEIN_FAA"

echo "[done] LAST protein index ready at ${LASTDB_INDEX}.*"