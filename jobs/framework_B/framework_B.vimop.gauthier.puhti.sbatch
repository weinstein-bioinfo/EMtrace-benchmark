#!/bin/bash
# ============================================================================
# Framework B — ViMOP on Gauthier 2021 SRA, Puhti SLURM array job
# ============================================================================
# One array task per Gauthier accession (21 SARS-CoV-2 SISPA libraries). Each
# task runs ViMOP (Nextflow + apptainer) on one FASTQ with the COVID virus
# filter. Independent samples run in parallel → wall ≈ the slowest single
# sample (heavy SISPA libs are hours), not the sum.
#
# Why Puhti (not local): the heavy libraries (e.g. SRR15356303 = 939 MB gz,
# millions of reads) need native x86_64 + lots of RAM; ViMOP's iterative Canu
# re-assemblies are memory-hungry and would OOM / crawl on the 16 GB laptop.
#
# Container engine: APPTAINER (Docker is unavailable on HPC). ViMOP's
# `apptainer` profile pulls its docker:// images into $APPTAINER_CACHEDIR.
#
# EDIT the placeholders marked  <<<  before submitting.
#
# Submit:    sbatch jobs/framework_B/framework_B.vimop.gauthier.puhti.sbatch
# Test one:  sbatch --array=1 jobs/framework_B/framework_B.vimop.gauthier.puhti.sbatch
# ----------------------------------------------------------------------------
#SBATCH --job-name=gauthier_vimop
#SBATCH --account=project_XXXXXXX          # <<< your CSC compute project
#SBATCH --partition=small
#SBATCH --array=1-21                        # 21 Gauthier accessions
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G                           # ViMOP/Canu re-assembly is memory-hungry
#SBATCH --time=08:00:00                     # heavy SISPA libs can take hours
#SBATCH --output=logs/gauthier_vimop_%A_%a.out
#SBATCH --error=logs/gauthier_vimop_%A_%a.err

set -euo pipefail

# === paths (EDIT) ===========================================================
PROJ=/scratch/project_XXXXXXX/emtrace-benchmark      # working root (inputs+results)
PROJAPPL=/projappl/project_XXXXXXX                            # software location
REPO="$PROJAPPL/EMtrace-benchmark"                               # this repo, cloned on Puhti
VIMOP_DB="$PROJ/ViMOP_DB"                                     # ViMOP DB staged on scratch (~39 GB)
# ===========================================================================

# Nextflow + apptainer from CSC modules. ViMOP's config needs the legacy parser.
module load nextflow apptainer 2>/dev/null || module load nextflow
export NXF_SYNTAX_PARSER=v1
export APPTAINER_CACHEDIR="$PROJ/apptainer_cache"   # image cache (large) on scratch, not $HOME
export NXF_APPTAINER_CACHEDIR="$APPTAINER_CACHEDIR"
mkdir -p "$APPTAINER_CACHEDIR" logs

ACCLIST="$REPO/jobs/assets/gauthier_accs.txt"  # one accession per line (21 lines)
SRA_DIR="$PROJ/sra"                                  # input FASTQs: <acc>.fastq.gz
RESDIR="$PROJ/results/vimop_out/gauthier"
STAGE_ROOT="$PROJ/vimop_in/gauthier"                 # per-task staging dirs
mkdir -p "$RESDIR" "$STAGE_ROOT"

# --- resolve this task's accession -----------------------------------------
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCLIST")
[[ -n "${ACC:-}" ]] || { echo "ERROR: no accession at line $SLURM_ARRAY_TASK_ID of $ACCLIST" >&2; exit 1; }

FASTQ="$SRA_DIR/${ACC}.fastq.gz"
[[ -f "$FASTQ" ]] || { echo "ERROR: input FASTQ not found: $FASTQ" >&2; exit 1; }

# resume: skip if the consensus already exists
if [[ -d "$RESDIR/${ACC}/selected_consensus" ]]; then
  echo "[skip] $ACC (selected_consensus/ exists)"; exit 0
fi

# stage the FASTQ into its own dir so ViMOP names the sample <ACC>
STAGE="$STAGE_ROOT/$ACC"
mkdir -p "$STAGE"
ln -f "$FASTQ" "$STAGE/" 2>/dev/null || cp -f "$FASTQ" "$STAGE/"

echo "[task $SLURM_ARRAY_TASK_ID] $ACC  on $(hostname)  ${SLURM_CPUS_PER_TASK} cpus"

# Execution provenance: per-task trace/report/timeline (citeable runtime + peak-mem
# artifacts). Per-ACC paths so the 21 array tasks never collide.
PROV="$RESDIR/${ACC}/provenance"
mkdir -p "$PROV"

# All Gauthier samples are SARS-CoV-2 → --targets COVID.
# Real Puhti resources clear ViMOP's preflight gates, so no --min_* overrides.
srun nextflow run opr-group-bnitm/vimop \
  -profile apptainer \
  --fastq   "$STAGE" \
  --out_dir "$RESDIR" \
  --base_db "$VIMOP_DB" \
  --targets COVID \
  -with-trace    "$PROV/trace.txt" \
  -with-report   "$PROV/report.html" \
  -with-timeline "$PROV/timeline.html" \
  -resume

echo "[done] $ACC -> $RESDIR/$ACC  (provenance: $PROV)"
