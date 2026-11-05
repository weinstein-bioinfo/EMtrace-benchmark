#!/bin/bash
# ============================================================================
# Build RefSeq abvh.1g (archaea+bacteria 1/genus + viruses + human) — INDEX
# Roihu SLURM job — stage 2 of 2 (mmi + pack)
# ============================================================================
# Run after db_refseq_abvh.dload.sbatch completes successfully. Builds:
#   - Minimap2 .mmi index
#   - tar.gz archive of the Minimap2 files, ready for upload to Allas
#
# No BLASTN index is built here: EMtrace uses nucleotide/minimap2,
# MEGAN-LR/metaMix use protein (full metaMix mode) — none of them need a
# nucleotide BLASTN db.
#
# EPI2ME does NOT use this script's .mmi output — it needs a Kraken2 index
# instead, built from the same dload output by db_refseq_abvh.kraken2.build.sbatch 
#
# EDIT the placeholders marked  <<<  before submitting
#
# Submit:  sbatch jobs/db_refseq_abvh.build.sbatch
# (once you've confirmed dload finished cleanly)
# ----------------------------------------------------------------------------
#SBATCH --job-name=refseq_abvh_build
#SBATCH --account=Project_2003755                # <<< your CSC project
#SBATCH --partition=small
#SBATCH --cpus-per-task=8                         
#SBATCH --mem=48G                                 
#SBATCH --time=04:00:00                           
#SBATCH --output=logs/db_refseq_abvh_build_%j.out
#SBATCH --error=logs/db_refseq_abvh_build_%j.err

set -euo pipefail

PROJAPPL=/projappl/project_XXXXXXX                    # software/repo location
LZPDB_REPO="$PROJAPPL/Lazypipe-db"
LZPDB_VENV="$PROJAPPL/venv/lazypipe-db"
RELEASE=234

# $databases must be set before running:
#   export databases=/scratch/project_XXXXXXX/emtrace-benchmark/db
: "${databases:?set \$databases (Roihu reference-DB root) before submitting}"

mkdir -p logs

# minimap2 must be on PATH:
export PATH="$PROJAPPL/install/minimap2:$PATH"
# (pigz is on PATH by default on Roihu)

source "$LZPDB_VENV/bin/activate"

echo "[build] target=abvh --genus  release=$RELEASE  on $(hostname)  ${SLURM_CPUS_PER_TASK} cpus"

srun python3 "$LZPDB_REPO/python/db_refseq.py" \
  --release "$RELEASE" \
  --target abvh \
  --genus \
  --pipe mmi,pack \
  --pack-mmi \
  --threads "$SLURM_CPUS_PER_TASK" \
  --config "$LZPDB_REPO/python/db_refseq.config.yaml" \
  -vv

echo "[done] build — Minimap2 archive written under $LZPDB_REPO's wrkdir (see db_refseq.config.yaml: wrkdir), ready to rclone to Allas"
