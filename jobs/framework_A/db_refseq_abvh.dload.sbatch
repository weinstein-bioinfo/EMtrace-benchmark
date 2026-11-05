#!/bin/bash
# ============================================================================
# Build RefSeq abvh.1g (archaea+bacteria 1/genus + viruses + human) — DOWNLOAD
# Roihu SLURM job — stage 1 of 2 (dload only)
# ============================================================================
# Calls Lazypipe-db's db_refseq.py (at $LZPDB_REPO) to build a database 
# which is shared across EMtrace / metaMix / MEGAN-LR.
#
# Split from mmi/pack to avoid crashing after download

# --target abvh --genus: archaea+bacteria filtered to Complete Genome +
#   reference/representative, collapsed to 1 genome per genus (viruses
#   unaffected by --genus) + human GRCh38 (GCA_000001405.15).
#
# Downstream consumers of this stage's combined FASTA + .acc2taxid:
#   db_refseq_abvh.build.sbatch          — minimap2 .mmi (EMtrace)
#   db_refseq_abvh.kraken2.build.sbatch  — Kraken2 + Bracken index (EPI2ME)
# Both read the same dload output; this stage only needs to run once.
#
# EDIT the placeholders marked  <<<  before submitting:
#   --account, --output/--error dirs
#
# Submit:  sbatch /jobs/framework_A/db_refseq_abvh.dload.sbatch
# When succesfully completed:    sbatch /jobs/framework_A/db_refseq_abvh.build.sbatch
# ----------------------------------------------------------------------------
#SBATCH --job-name=refseq_abvh_dload
#SBATCH --account=Project_2003755                # <<< your CSC project
#SBATCH --partition=small
#SBATCH --cpus-per-task=8                        # matches download_workers=8 in db_refseq.config.yaml
#SBATCH --mem=64G                                
#SBATCH --time=24:00:00                           
#SBATCH --output=logs/db_refseq_abvh_dload_%j.out
#SBATCH --error=logs/db_refseq_abvh_dload_%j.err

set -euo pipefail

PROJAPPL=/projappl/project_XXXXXXX
LZPDB_REPO="$PROJAPPL/Lazypipe-db"
LZPDB_VENV="$PROJAPPL/venv/lazypipe-db"
RELEASE=234

# $databases must be set before running:
#   export databases=/scratch/project_XXXXXXX/emtrace-benchmark/db
: "${databases:?set \$databases (Roihu reference-DB root) before submitting}"

mkdir -p logs

# wget, pigz, seqkit must be on PATH. wget and pigz are on PATH by default on
# Roihu. seqkit is not (no biokit module on Roihu) — add it:
export PATH="$PROJAPPL/install/seqkit_dir/bin:$PATH"

source "$LZPDB_VENV/bin/activate"

echo "[dload] target=abvh --genus  release=$RELEASE  on $(hostname)  ${SLURM_CPUS_PER_TASK} cpus"

srun python3 "$LZPDB_REPO/python/db_refseq.py" \
  --release "$RELEASE" \
  --target abvh \
  --genus \
  --pipe dload \
  --threads "$SLURM_CPUS_PER_TASK" \
  --config "$LZPDB_REPO/python/db_refseq.config.yaml" \
  -vv

echo "[done] dload — check logs/ for warnings before submitting the build stage"
