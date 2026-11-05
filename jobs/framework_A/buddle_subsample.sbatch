#!/bin/bash
# ============================================================================
# Subsample the active Buddle 2024 DNA/RNA accessions to 5 Gb each.
# ============================================================================
# Reads the active sample list straight from config.yaml (buddle_2024_dna +
# buddle_2024_rna only — datasets/reps commented out there are skipped).
#
# For each accession:
#   1. streams gunzip | subsample_buddle.py directly
#   2. moves the full-size original aside to sra/full_originals/<acc>.fastq.gz
#   3. moves the new 5Gb.fastq.gz into sra/<acc>.fastq.gz — the exact path every
#      framework_A/<pipeline>/*.sbatch script already expects
#
# Submit: sbatch jobs/framework_A/buddle_subsample.sbatch
# ----------------------------------------------------------------------------
#SBATCH --job-name=buddle_subsample
#SBATCH --account=Project_2003755
#SBATCH --partition=small
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/buddle_subsample_%j.out
#SBATCH --error=logs/buddle_subsample_%j.err

set -euo pipefail

REPO="/projappl/project_XXXXXXX/EMtrace-benchmark"
CONFIG="$REPO/jobs/config.yaml"
INPUT_DIR="/scratch/project_XXXXXXX/emtrace-benchmark/sra"
ARCHIVE_DIR="${INPUT_DIR}/full_originals"   # full-size originals preserved here, not deleted
SUBSAMPLE_SCRIPT="$REPO/jobs/framework_A/subsample_buddle.py"
TARGET_BASES=5000000000

module load python-data

mkdir -p "$ARCHIVE_DIR" logs

# ---------------------------------------------------------------------------
# Active Buddle accessions, straight from config.yaml
# Only buddle_2024_dna / buddle_2024_rna
# ---------------------------------------------------------------------------
mapfile -t SAMPLES < <(python3 - "$CONFIG" <<'PY'
import sys, yaml
with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)
seen = set()
for name in ("buddle_2024_dna", "buddle_2024_rna"):
    ds = (cfg.get("datasets") or {}).get(name) or {}
    for sample in (ds.get("samples") or []):
        for key in ("sra", "nc_sra"):
            acc = sample.get(key)
            if not acc or acc == "PLACEHOLDER" or acc in seen:
                continue
            seen.add(acc)
            print(acc)
PY
)

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "ERROR: no active accessions found in $CONFIG for buddle_2024_dna/rna" >&2
    exit 1
fi

echo "Started $(date) — ${#SAMPLES[@]} accession(s): ${SAMPLES[*]}"

for ACC in "${SAMPLES[@]}"; do
    SRC="${INPUT_DIR}/${ACC}.fastq.gz"
    ARCHIVED="${ARCHIVE_DIR}/${ACC}.fastq.gz"
    STAGED="${INPUT_DIR}/.${ACC}.5gb.fastq"

    if [[ -f "$ARCHIVED" ]]; then
        echo "[skip] $ACC (already subsampled — original archived at $ARCHIVED)"
        continue
    fi

    if [[ ! -f "$SRC" ]]; then
        echo "WARNING: $ACC not found at $SRC — skipping" >&2
        continue
    fi

    echo "[$ACC] subsampling"
    set +o pipefail
    if gunzip -c "$SRC" | python3 "$SUBSAMPLE_SCRIPT" /dev/stdin "$STAGED" "$TARGET_BASES"; then
        set -o pipefail
    else
        set -o pipefail
        echo "ERROR: subsample failed for $ACC — original left untouched" >&2
        rm -f "$STAGED"
        continue
    fi

    gzip -f "$STAGED"   # -> ${STAGED}.gz

    if [[ ! -s "${STAGED}.gz" ]]; then
        echo "ERROR: subsampled output empty for $ACC — original left untouched" >&2
        rm -f "${STAGED}.gz"
        continue
    fi

    mv "$SRC" "$ARCHIVED"          # preserve the full-size original
    mv "${STAGED}.gz" "$SRC"       # swap the 5Gb version into the path emtrace reads

    echo "[done] $ACC -> $SRC (5Gb); original archived -> $ARCHIVED"
done

echo
echo "Finished $(date)"
echo
echo "Full-size originals preserved under: $ARCHIVE_DIR"
