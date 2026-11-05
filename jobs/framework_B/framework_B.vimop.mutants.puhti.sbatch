#!/bin/bash
# ============================================================================
# Framework B — ViMOP on simulated mutant datasets (pure + mixed), Puhti array
# ============================================================================
# One array task per mutant sample: 4 viruses × 7 deviations × {pure, mixed}
# = 56 tasks. Each runs ViMOP (Nextflow + apptainer) with the matching virus
# target filter.
#
# Why Puhti: the MIXED sets embed the target in ~500k background reads each
# (LASV_L_25_mixed = 500,333 reads) → ViMOP's iterative Canu re-assembly is
# heavy and memory-hungry. The --targets filter is ESSENTIAL here: it strips
# the background before assembly (~halves runtime, avoids OOM). Single-target,
# so no co-infection loss.
#
# Inputs (already on disk in the repo, copy to $PROJ on Puhti):
#   pure :  data/petersen_2026/fastq/<V>_<DD>/<V>_<DD>.fastq.gz   (own dir)
#   mixed:  data/petersen_2026/mixed/<V>_<DD>_mixed.fastq.gz      (flat → staged)
# Target filter per virus: COVID→COVID, DENV→DENV, LASV_L→LASV, LASV_S→LASV.
#
# EDIT the placeholders marked  <<<  before submitting.
# Submit:    sbatch jobs/framework_B/framework_B.vimop.mutants.puhti.sbatch
# Test one:  sbatch --array=1 jobs/framework_B/framework_B.vimop.mutants.puhti.sbatch
# ----------------------------------------------------------------------------
#SBATCH --job-name=mutants_vimop
#SBATCH --account=project_XXXXXXX          # <<< your CSC compute project
#SBATCH --partition=small
#SBATCH --array=1-56                        # 4 viruses × 7 devs × {pure,mixed}
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/mutants_vimop_%A_%a.out
#SBATCH --error=logs/mutants_vimop_%A_%a.err

set -euo pipefail

# === paths (EDIT) ===========================================================
PROJ=/scratch/project_XXXXXXX/emtrace-benchmark
PROJAPPL=/projappl/project_XXXXXXX
REPO="$PROJAPPL/EMtrace-benchmark"
VIMOP_DB="$PROJ/ViMOP_DB"
DATA="$PROJ/data/petersen_2026"     # pure: $DATA/fastq/<s>/<s>.fastq.gz ; mixed: $DATA/mixed/<s>_mixed.fastq.gz
# ===========================================================================

module load nextflow apptainer 2>/dev/null || module load nextflow
export NXF_SYNTAX_PARSER=v1
export APPTAINER_CACHEDIR="$PROJ/apptainer_cache"
export NXF_APPTAINER_CACHEDIR="$APPTAINER_CACHEDIR"
mkdir -p "$APPTAINER_CACHEDIR" logs

RESDIR="$PROJ/results/vimop_out"
STAGE_ROOT="$PROJ/vimop_in/mutants"
mkdir -p "$RESDIR" "$STAGE_ROOT"

# --- build the deterministic task list (same order in every array task) -----
VIRUSES=(COVID DENV LASV_L LASV_S)
DEVS=(00 05 10 15 20 25 30)
TASKS=()
for set in pure mixed; do
  for v in "${VIRUSES[@]}"; do
    for d in "${DEVS[@]}"; do
      TASKS+=("$set ${v}_${d}")
    done
  done
done

idx=$(( SLURM_ARRAY_TASK_ID - 1 ))
[[ "$idx" -ge 0 && "$idx" -lt "${#TASKS[@]}" ]] || { echo "ERROR: array idx $SLURM_ARRAY_TASK_ID out of range (1-${#TASKS[@]})" >&2; exit 1; }
read -r SET SAMPLE <<< "${TASKS[$idx]}"

# virus → target filter (strip trailing _DD; LASV_L/LASV_S both use the LASV filter)
virus="${SAMPLE%_[0-9][0-9]}"
case "$virus" in
  COVID)            TARGET=COVID ;;
  DENV)             TARGET=DENV ;;
  LASV_L|LASV_S)    TARGET=LASV ;;
  *) echo "ERROR: cannot map virus '$virus' (sample $SAMPLE) to a ViMOP --targets" >&2; exit 1 ;;
esac

# resolve input + ViMOP --fastq dir
if [[ "$SET" == pure ]]; then
  NAME="$SAMPLE"                                  # e.g. COVID_25
  SRC="$DATA/fastq/$NAME/$NAME.fastq.gz"
  FASTQ_DIR="$DATA/fastq/$NAME"                   # already a per-sample dir
else
  NAME="${SAMPLE}_mixed"                          # e.g. COVID_25_mixed
  SRC="$DATA/mixed/$NAME.fastq.gz"
  FASTQ_DIR="$STAGE_ROOT/$NAME"                   # flat input → stage into own dir
  mkdir -p "$FASTQ_DIR"
  ln -f "$SRC" "$FASTQ_DIR/" 2>/dev/null || cp -f "$SRC" "$FASTQ_DIR/"
fi

[[ -f "$SRC" ]] || { echo "[skip] $SET/$NAME (input not found: $SRC)"; exit 0; }
if [[ -d "$RESDIR/$SET/$NAME/selected_consensus" ]]; then
  echo "[skip] $SET/$NAME (selected_consensus/ exists)"; exit 0
fi

echo "[task $SLURM_ARRAY_TASK_ID] $SET/$NAME  --targets $TARGET  on $(hostname)  ${SLURM_CPUS_PER_TASK} cpus"

# Execution provenance: per-task trace/report/timeline (citeable runtime + peak-mem
# artifacts). Per-sample paths so the 56 array tasks never collide.
PROV="$RESDIR/$SET/$NAME/provenance"
mkdir -p "$PROV"

srun nextflow run opr-group-bnitm/vimop \
  -profile apptainer \
  --fastq   "$FASTQ_DIR" \
  --out_dir "$RESDIR/$SET" \
  --base_db "$VIMOP_DB" \
  --targets "$TARGET" \
  -with-trace    "$PROV/trace.txt" \
  -with-report   "$PROV/report.html" \
  -with-timeline "$PROV/timeline.html" \
  -resume

echo "[done] $SET/$NAME -> $RESDIR/$SET/$NAME  (provenance: $PROV)"
