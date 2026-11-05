#!/bin/bash -l
# ============================================================================
# Build RefSeq abvh.1g protein DB (for metaMix DIAMOND + MEGAN-LR LAST) — DOWNLOAD
# Roihu SLURM job — stage 1 of 2 (dload only)
# ============================================================================
# Same organism selection as the nucleotide build (db_refseq_abvh.dload.sbatch /
# .build.sbatch: viral + archaea/bacteria 1/genus + human, release 234), but
# --seqtype prot fetches each assembly's *_protein.faa.gz instead of
# *_genomic.fna.gz (db_refseq.py, extended 2026-08 to support --seqtype).
#
# EDIT the placeholders marked  <<<  before submitting.
#
# Submit:  sbatch jobs/framework_A/db_refseq_abvh_protein.dload.sbatch
# Then (after succesfully completed):    
# sbatch jobs/framework_A/db_refseq_abvh_protein.build.sbatch
# ----------------------------------------------------------------------------
#SBATCH --job-name=refseq_abvh_prot_dload
#SBATCH --account=Project_2003755                # <<< your CSC project
#SBATCH --partition=small
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G                                 
#SBATCH --time=24:00:00                           
#SBATCH --output=logs/db_refseq_abvh_prot_dload_%j.out
#SBATCH --error=logs/db_refseq_abvh_prot_dload_%j.err

set -euo pipefail

PROJAPPL=/projappl/project_XXXXXXX                    # software/repo location
LZPDB_REPO="$PROJAPPL/Lazypipe-db"
LZPDB_VENV="$PROJAPPL/venv/lazypipe-db"
RELEASE=234

# $databases must be set before running:
#   export databases=/scratch/project_XXXXXXX/emtrace-benchmark/db
: "${databases:?set \$databases (Roihu reference-DB root) before submitting}"

mkdir -p logs

# wget, pigz, seqkit must be on PATH (see db_refseq_abvh.dload.sbatch). wget and
# pigz are on PATH by default on Roihu; seqkit is not (no biokit module):
export PATH="$PROJAPPL/install/seqkit_dir/bin:$PATH"

source "$LZPDB_VENV/bin/activate"

echo "[dload] target=abvh --genus --seqtype prot  release=$RELEASE  on $(hostname)  ${SLURM_CPUS_PER_TASK} cpus"

srun python3 "$LZPDB_REPO/python/db_refseq.py" \
  --release "$RELEASE" \
  --target abvh \
  --genus \
  --seqtype prot \
  --pipe dload \
  --threads "$SLURM_CPUS_PER_TASK" \
  --config "$LZPDB_REPO/python/db_refseq.config.yaml" \
  -vv

echo "[done] dload — check logs/ for warnings before submitting the build stage"
