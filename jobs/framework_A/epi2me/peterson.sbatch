#!/bin/bash
# =============================================================================
# Framework A — EPI2ME wf-metagenomics baseline comparator, Peterson 2026 (Roihu)
# =============================================================================
# Runs epi2me-labs/wf-metagenomics (Kraken2 + Bracken, ONT-native) on the RAW
# fetched FASTQ for one accession (this matches Buddle et al. 2024's 
# methodology)
#  
# This script runs epi2me's own bundled --database_set (Standard-8, see below).
#
# Output: wf-metagenomics writes bracken/<ACC>.kraken2_bracken.report (the
# per-sample Bracken TSV — name/taxonomy_id/taxonomy_lvl/kraken_assigned_reads/
# added_reads/new_est_reads/fraction_total_reads columns) under --out_dir.
#
# USAGE: sbatch framework_A/epi2me/peterson.sbatch
# Test one:
#   sbatch --array=1 framework_A/epi2me/peterson.sbatch
# =============================================================================
#SBATCH --job-name=epi2me_peterson
#SBATCH --account=Project_2003755
#SBATCH --partition=small
#SBATCH --array=1-4                         # 4 accs (peterson_2026, no NC)
#SBATCH --cpus-per-task=8                   
#SBATCH --mem=32G                           
#SBATCH --time=02:00:00                     
#SBATCH --output=logs/epi2me_peterson_%A_%a.out
#SBATCH --error=logs/epi2me_peterson_%A_%a.err

set -euo pipefail

PROJ="/scratch/project_XXXXXXX/emtrace-benchmark"
PROJAPPL="/projappl/project_XXXXXXX"
REPO="$PROJAPPL/EMtrace-benchmark"
SRA_DIR="$PROJ/sra"

module load nextflow

# Apptainer directories (preinstalled on Roihu, no module load needed — see docs.csc.fi/apps/nextflow).
# NXF_SINGULARITY_CACHEDIR is shared/persistent (Nextflow reuses converted .sif images across
# runs and array tasks instead of re-pulling each time); APPTAINER_TMPDIR/CACHEDIR are node-local
# scratch, used only while an image is actively being pulled/built.
export NXF_SINGULARITY_CACHEDIR="$PROJ/singularity_cache"
export APPTAINER_TMPDIR="${LOCAL_SCRATCH:-${TMPDIR:-/tmp}}"
export APPTAINER_CACHEDIR="${LOCAL_SCRATCH:-${TMPDIR:-/tmp}}"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# --store_dir: wf-metagenomics' own persistent cache for downloaded/unpacked reference
# data (Kraken2 --database_set index, taxonomy dump), separate from -work-dir and from
# NXF_SINGULARITY_CACHEDIR above. Without this pinned to a SHARED absolute path, it
# defaults to a bare relative "store_dir" — which, under per-ACC LAUNCHDIR (below),
# means every array task re-downloads+unpacks the full Standard-8 DB into its own 
# private directory instead of reusing one copy. 
# Shared across all epi2me datasets since --database_set Standard-8 is the same DB 
# in every framework_A/epi2me/*.sbatch file.
STORE_DIR="$PROJ/epi2me_store"
mkdir -p "$STORE_DIR"

ACCLIST="$REPO/jobs/assets/peterson_accs.txt"
EXCLUDE_HOST=0            # <<< 0 = raw reads (default, matches Buddle et al. 2024 + baseline framing) | 1 = host-deplete via --exclude_host

ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCLIST")
if [[ -z "$ACC" ]]; then
  echo "[error] no accession at line $SLURM_ARRAY_TASK_ID of $ACCLIST" >&2
  exit 1
fi

FASTQ="$SRA_DIR/${ACC}.fastq.gz"
if [[ ! -f "$FASTQ" ]]; then
  echo "[error] $FASTQ not found" >&2
  exit 1
fi

RESDIR="$PROJ/results/framework_a/epi2me/peterson_2026/$ACC"
WORKDIR="$PROJ/work/epi2me/$ACC"
# Nextflow's session/cache lock (.nextflow/cache/<session>/db/LOCK) lives under the
# LAUNCH directory (i.e. $PWD when `nextflow run` is invoked), not under -work-dir.
# Without a per-task launch dir, every array task shares the same $PROJ launchDir and
# collides on that lock file when run concurrently. Isolate it per ACC.
LAUNCHDIR="$PROJ/work/epi2me/$ACC/.launch"
mkdir -p "$RESDIR" "$LAUNCHDIR" logs
cd "$LAUNCHDIR"

BRACKEN_OUT="$RESDIR/${ACC}_epi2me.bracken.tsv"
if [[ -f "$BRACKEN_OUT" ]]; then
  echo "[skip] $BRACKEN_OUT already present"
  exit 0
fi

EXTRA_ARGS=()
if [[ "$EXCLUDE_HOST" == "1" ]]; then
  : "${hostgenomes:?set \$hostgenomes (Roihu host-genome root) when EXCLUDE_HOST=1}"
  EXTRA_ARGS+=(--exclude_host "$hostgenomes/GCF_000001405.40_genomic+rna.fna.gz")
fi

# Serialize until store_dir's DB is first populated.
if [[ -z "$(ls -A "$STORE_DIR/Standard-8" 2>/dev/null)" ]]; then
  echo "[info] Standard-8 DB not yet cached in $STORE_DIR — serializing until the first task finishes populating it"
  exec 9>"$STORE_DIR/.populate.lock"
  flock 9
fi

echo "[run] epi2me-labs/wf-metagenomics: $ACC"
nextflow run epi2me-labs/wf-metagenomics \
  --fastq "$FASTQ" \
  --sample "$ACC" \
  --classifier kraken2 \
  --database_set Standard-8 \
  --store_dir "$STORE_DIR" \
  --out_dir "$RESDIR" \
  --threads "$SLURM_CPUS_PER_TASK" \
  --disable_ping \
  -profile singularity \
  -work-dir "$WORKDIR" \
  -resume \
  "${EXTRA_ARGS[@]}"

REPORT="$RESDIR/bracken/${ACC}.kraken2_bracken.report"
if [[ ! -f "$REPORT" ]]; then
  echo "[error] expected $REPORT after wf-metagenomics run, not found" >&2
  exit 1
fi
cp "$REPORT" "$BRACKEN_OUT"
rm -rf "$WORKDIR" "$LAUNCHDIR"

echo "[done] $ACC -> $BRACKEN_OUT"
