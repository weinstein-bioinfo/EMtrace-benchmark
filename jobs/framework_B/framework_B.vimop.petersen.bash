#!/usr/bin/env bash
# ============================================================================
# Framework B — ViMOP consensus assembly on the Petersen 2026 SRA set (LOCAL)
# ============================================================================
#
# Runs ViMOP (Nextflow + Docker) on the 4 Petersen samples on THIS machine.
# Only Petersen runs locally; Gauthier + mutants run on Puhti (see the
# framework_B.vimop.*.puhti.sbatch scripts).
#
# Per sample we apply ViMOP's virus target filter (--targets) — these are
# single-target libraries, so the filter ~halves runtime and cuts memory with
# no loss (co-infection detection is not needed). Filter names come from
# ~/ViMOP_DB/virus/<NAME>.fasta.
#
#   SRR31174237  Zika        --targets ZIKA   (96,903 reads)
#   ERR14941561  West Nile    --targets WNV    (403 reads — already done)
#   SRR15356294  SARS-CoV-2   --targets COVID  (85,593 reads)
#   SRR32419728  Dengue-2     --targets DENV   (297,971 reads)
#
# ViMOP params (see the project notes):
#   --fastq   : a DIRECTORY of FASTQs (not a file); the dir basename becomes the
#               sample/output name, so we stage each SRA into its own dir.
#   --out_dir : parent dir; ViMOP writes <out_dir>/<acc>/ per sample.
#   NXF_SYNTAX_PARSER=v1 : ViMOP's config needs Nextflow's legacy parser.
#   --min_ram_gb 16 --min_cpus 12 : lower the preflight gates to this Mac.
#
# Resume: a sample is skipped if its selected_consensus/ dir already exists.
#
# RUN ME under caffeinate + nohup so a sleeping laptop doesn't kill the loop:
#   caffeinate -ism nohup bash jobs/framework_B/framework_B.vimop.petersen.bash \
#     > logs/vimop/petersen.$(date +%Y%m%d_%H%M%S).out 2>&1 &
#
# Score afterwards with: bash jobs/framework_B/framework_B.vimop.metrics.bash
# ============================================================================
set -uo pipefail

export DOCKER_DEFAULT_PLATFORM=linux/amd64
export NXF_SYNTAX_PARSER=v1

REPO_ROOT="$(pwd)"
cd "$REPO_ROOT" || { echo "cannot cd to $REPO_ROOT"; exit 1; }

SRA_DIR="data/sra"
STAGE_DIR="data/vimop_in"
OUT_DIR="results/vimop_out/petersen"
LOG_DIR="logs/vimop/petersen"
mkdir -p "$STAGE_DIR" "$OUT_DIR" "$LOG_DIR"

# sample -> ViMOP virus filter target
declare -A TARGETS=(
  [SRR31174237]=ZIKA     # Zika
  [ERR14941561]=WNV      # West Nile
  [SRR15356294]=COVID    # SARS-CoV-2
  [SRR32419728]=DENV     # Dengue-2
)

# resolve FASTQ path (.fastq then .fastq.gz)
_fastq() { local b="$SRA_DIR/$1"; [[ -f "${b}.fastq" ]] && echo "${b}.fastq" || echo "${b}.fastq.gz"; }

for SRA in "${!TARGETS[@]}"; do
  TARGET="${TARGETS[$SRA]}"
  FASTQ=$(_fastq "$SRA")

  if [[ ! -f "$FASTQ" ]]; then
    echo "[skip] $SRA (FASTQ not found: $FASTQ — fetch with fasterq-dump $SRA -O $SRA_DIR/)"
    continue
  fi
  if [[ -d "$OUT_DIR/$SRA/selected_consensus" ]]; then
    echo "[skip] $SRA (selected_consensus/ exists)"
    continue
  fi

  # stage the single FASTQ into its own dir so ViMOP names the sample <SRA>
  mkdir -p "$STAGE_DIR/$SRA"
  ln -f "$FASTQ" "$STAGE_DIR/$SRA/" 2>/dev/null || cp -f "$FASTQ" "$STAGE_DIR/$SRA/"

  echo ">>> ViMOP $SRA  (--targets $TARGET)  $(date '+%H:%M:%S')"
  # -ansi-log false + </dev/null: never touch the controlling TTY, so a
  # backgrounded run isn't SIGTTOU-stopped (plain line logs go to the tee file).
  nextflow run opr-group-bnitm/vimop \
    --fastq   "$STAGE_DIR/$SRA" \
    --out_dir "$OUT_DIR" \
    --targets "$TARGET" \
    --min_ram_gb 16 --min_cpus 12 \
    -ansi-log false \
    -resume \
    </dev/null 2>&1 | tee "$LOG_DIR/${SRA}.log"

  if [[ -d "$OUT_DIR/$SRA/selected_consensus" ]]; then
    echo "[done] $SRA -> $OUT_DIR/$SRA"
  else
    echo "[WARN] $SRA finished but no selected_consensus/ — check $LOG_DIR/${SRA}.log"
  fi
done

echo "All Petersen ViMOP runs attempted. Score with: bash jobs/framework_B/framework_B.vimop.metrics.bash"
