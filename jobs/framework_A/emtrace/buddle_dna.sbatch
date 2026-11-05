#!/bin/bash
# ============================================================================
# Framework A — EMtrace on Buddle 2024 DNA viruses, Roihu SLURM array job
# ============================================================================
# One array task per accession (5 samples: 4 positives + 1 NC, buddle_dna_a_accs.txt).
# Each task runs the full emtrace pipeline on one FASTQ. Independent samples
# → they run in parallel.
#
# EDIT the placeholders marked  <<<  before submitting:
#   --account, PROJ, PROJAPPL, REPO, VENV, EMTRACE_DB_DIR, --job-name
#
# Submit:   sbatch jobs/framework_A/emtrace/buddle_dna.sbatch
# Test one: sbatch --array=1 jobs/framework_A/emtrace/buddle_dna.sbatch
#
# ----------------------------------------------------------------------------
#SBATCH --job-name=emtrace_buddle_dna
#SBATCH --account=Project_2003755           # <<< your CSC compute project
#SBATCH --partition=small
#SBATCH --array=1-5                         # 5 accessions (4 samples + 1 NC)
#SBATCH --cpus-per-task=24                  # emtrace -t 24
#SBATCH --mem=80G                           # human .mmi (~8.5G) + ABV index + overhead
#SBATCH --time=02:00:00                     # ~30 min/sample; margin for 1st-run indexing
#SBATCH --output=logs/buddle_dna_%A_%a.out    # %A=jobid, %a=array idx
#SBATCH --error=logs/buddle_dna_%A_%a.err

set -euo pipefail

# === paths (EDIT) ===========================================================
PROJ=/scratch/project_XXXXXXX/emtrace-benchmark      # working root (inputs+results)
PROJAPPL=/projappl/project_XXXXXXX                            # software/repo location
REPO="$PROJAPPL/EMtrace-benchmark"                      # this repo, cloned on Roihu
VENV="$PROJAPPL/venv/emtrace"                     # Python venv: emtrace (+ porechop_abi, orfipy, ...)
# Reference DBs live under Roihu's $databases; host genomes under $hostgenomes
# (both set in your Roihu environment). emtrace expands these env vars directly
# in config paths, so we run with a Roihu config (CONFIG below) whose host-genome
# paths point at $hostgenomes and whose annotation/taxonomy paths point at
# $databases. EMTRACE_DB_DIR is set as a fallback for any unoverridden path.
: "${databases:?set \$databases (Roihu reference-DB root) before submitting}"
: "${hostgenomes:?set \$hostgenomes (Roihu host-genome root) before submitting}"
export EMTRACE_DB_DIR="$databases"
# ===========================================================================

# export databases=/scratch/project_XXXXXXX/emtrace-benchmark/db
# export hostgenomes=/scratch/project_XXXXXXX/emtrace-benchmark/db/hostgenomes

ACCLIST="$REPO/jobs/assets/buddle_dna_a_accs.txt"  # one accession per line (5 lines)

# Add dependencies to PATH
export PATH="$PROJAPPL/install/minimap2:$PATH"
export PATH="$PROJAPPL/install/chopper_dir:$PATH"
export PATH="$PROJAPPL/install/samtools-1.23.1:$PATH"
export PATH="$PROJAPPL/install/bcftools-1.23.1:$PATH"
export PATH="$PROJAPPL/install/seqkit_dir/bin:$PATH"


# Activate the emtrace venv (puts emtrace + its Python entry points on PATH).
source "$VENV/bin/activate"
# The venv has the Python tools (emtrace, porechop_abi, orfipy) but NOT the
# external binaries emtrace shells out to — make these available on PATH or the
# run fails (minimap2 at host-filter/annotation, chopper at QC, samtools+bcftools
# at assembly, seqkit). Load them however they're provided on Roihu, e.g.:
# module load biokit            # <<< samtools / bcftools / seqkit
# module load minimap2 chopper  # <<< or your modules / another env for these

SRA_DIR="$PROJ/sra"                                  # input FASTQs: <acc>.fastq.gz
RESDIR="$PROJ/results"
BENCH="$REPO/jobs/bench_time.py"
TIMING="$RESDIR/metrics/buddle_dna_a.timings.csv"

mkdir -p "$RESDIR/metrics" logs

cd "$RESDIR"

# --- Roihu config: ALWAYS generated from the repo config (never forget it) ---
# Rewrite the single $EMTRACE_DB_DIR root into Roihu's two roots:
#   host genomes  -> $hostgenomes   (e.g. /scratch/.../lazydata/hostgen)
#   everything else (annotation, taxonomy) -> $databases
# emtrace expands $hostgenomes / $databases itself at load time. Written to a
# per-task temp file so the 5 array tasks never race on a shared file.
CONFIG="${LOCAL_SCRATCH:-${TMPDIR:-/tmp}}/emtrace.config.roihu.${SLURM_ARRAY_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-0}.yaml"
# shellcheck disable=SC2016  # $hostgenomes/$databases must stay literal — emtrace expands them
sed -e 's#\$EMTRACE_DB_DIR/hostgenomes/#$hostgenomes/#g' \
    -e 's#\$EMTRACE_DB_DIR/#$databases/#g' \
    "$REPO/jobs/emtrace.config.yaml" > "$CONFIG"

# --- resolve this task's accession -----------------------------------------
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCLIST")
if [[ -z "${ACC:-}" ]]; then
  echo "ERROR: no accession at line $SLURM_ARRAY_TASK_ID of $ACCLIST" >&2
  exit 1
fi

FASTQ="$SRA_DIR/${ACC}.fastq.gz"
if [[ ! -f "$FASTQ" ]]; then
  echo "ERROR: input FASTQ not found: $FASTQ" >&2
  exit 1
fi

# resume: skip if results already exists
if [[ -f "$RESDIR/${ACC}/results/reports/abund_table.tsv" ]]; then
  echo "[skip] $ACC (results/reports/abund_table.tsv exists)"
  exit 0
fi

echo "[task $SLURM_ARRAY_TASK_ID] $ACC  on $(hostname)  ${SLURM_CPUS_PER_TASK} cpus"

srun python3 "$BENCH" --csv "$TIMING" --label "$ACC" -- \
  emtrace run \
    -1 "$FASTQ" \
    -s "$ACC" \
    -r "$RESDIR" \
    -t "$SLURM_CPUS_PER_TASK" \
    --bgfilter Homo_sapiens \
    --ann1 minimap.refseq.abv.1g \
    --config "$CONFIG" \
    --pipe main \
    --logs "$RESDIR/${ACC}/pipe.log" \
    -v

echo "[done] $ACC -> $RESDIR/$ACC"
