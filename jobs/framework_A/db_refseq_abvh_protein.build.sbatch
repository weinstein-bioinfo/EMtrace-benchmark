#!/bin/bash -l
# ============================================================================
# Build RefSeq abvh.1g protein DB (for metaMix DIAMOND + MEGAN-LR LAST) — INDEX
# Roihu SLURM job — stage 2 of 2 (diamonddb + pack)
# ============================================================================
# Run after db_refseq_abvh_protein.dload.sbatch completes. Builds the DIAMOND
# index (diamond makedb --taxonmap prot.accession2taxid.gz ...) and packs the
# combined FASTA + .dmnd into one tar.gz.
#
# Output: $databases/2026/NCBI/refseq.abvh.1g.234/diamond/diamond.refseq.abvh.1g.dmnd
# (the combined protein FASTA lands one level up, at
# $databases/2026/NCBI/refseq.abvh.1g.234/diamond.refseq.abvh.1g.faa — the
# diamond/ subdir holds only the .dmnd index, per Lazypipe-db's build_paths())
#
# Side effect: this run's --taxonmap step also fetches nodes.dmp/names.dmp
# into $databases/2026/NCBI/refseq.abvh.1g.234/taxdump/ — db_refseq_abvh.kraken2.build.sbatch
# (builds the Kraken2 DB EPI2ME uses) reuses that same taxdump if this script
# has already run, instead of re-fetching it.
#
# EDIT the placeholders marked  <<<  before submitting 
#
# Submit:  sbatch jobs/framework_A/db_refseq_abvh_protein.build.sbatch
# ----------------------------------------------------------------------------
#SBATCH --job-name=refseq_abvh_prot_build
#SBATCH --account=Project_2003755                # <<< your CSC project
#SBATCH --partition=small
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G                                 
#SBATCH --time=08:00:00                         
#SBATCH --output=logs/db_refseq_abvh_prot_build_%j.out
#SBATCH --error=logs/db_refseq_abvh_prot_build_%j.err

set -euo pipefail

PROJAPPL=/projappl/project_XXXXXXX                    # software/repo location
LZPDB_REPO="$PROJAPPL/Lazypipe-db"
LZPDB_VENV="$PROJAPPL/venv/lazypipe-db"
RELEASE=234

# $databases must be set before running:
#   export databases=/scratch/project_XXXXXXX/emtrace-benchmark/db
: "${databases:?set \$databases (Roihu reference-DB root) before submitting}"

mkdir -p logs

# Roihu's diamond module is broken (prepends a nonexistent "diamiond" typo'd
# path even after loading successfully) — use the manually-installed binary
# instead.
export PATH="$PROJAPPL/install/diamond_dir:$PATH"
# (pigz is on PATH by default on Roihu)

source "$LZPDB_VENV/bin/activate"

echo "[build] target=abvh --genus --seqtype prot  release=$RELEASE  on $(hostname)  ${SLURM_CPUS_PER_TASK} cpus"

srun python3 "$LZPDB_REPO/python/db_refseq.py" \
  --release "$RELEASE" \
  --target abvh \
  --genus \
  --seqtype prot \
  --pipe diamonddb,pack \
  --threads "$SLURM_CPUS_PER_TASK" \
  --config "$LZPDB_REPO/python/db_refseq.config.yaml" \
  -vv

echo "[done] build — diamond.refseq.abvh.1g.dmnd ready under \$databases/2026/NCBI/refseq.abvh.1g.234/diamond/"
