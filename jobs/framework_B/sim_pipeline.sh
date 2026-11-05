#!/usr/bin/env bash
#
# sim_pipeline.sh — Framework B genome simulation orchestrator.
#
# Chains the five simulation steps in sequence:
#   sim_01  B_sim01_build_multifastas.py     fetch refs + build multifastas (local mode only)
#   sim_02  B_sim02_align_and_weights.py     alignment + positional variability weights
#   sim_03  B_sim03_generate_mutants.py      generate artificial genomes (28 FASTAs)
#   sim_04  B_sim04_simulate_reads.py        PBSIM3 read simulation (28 × FASTQ.gz)
#   sim_05  B_sim05_mix_background.py        build background-only + mixed datasets
#
# Before the background step (sim_04 --mode background, sim_05), download the
# background genomes once with B_sim03_bg_download.py.
#
# Usage:
#   bash jobs/framework_B/sim_pipeline.sh
#        [--targets COVID DENV LASV_L LASV_S]
#        [--source vimop|local]     # alignment source for sim_02 (default: vimop)
#        [--mode pure|background|both]  # read output mode for sim_05 (default: both)
#        [--skip 1,4,5]             # comma-separated step numbers to skip
#        [--seed 42] [--threads 8] [--depth 30]
#
# When --source vimop (default), step 1 is automatically skipped since ViMOP's
# curated alignments replace the GenBank fetch + MAFFT pipeline.
#
# Outputs:
#   data/petersen_2026/weights/     — positional variability weights
#   data/petersen_2026/mutants/     — 28 artificial genome FASTAs
#   data/petersen_2026/fastq/       — 28 simulated read sets (sim_04)
#   data/petersen_2026/mixed/       — background-embedded datasets (sim_05)
#   logs/sim_pipeline_<date>.log
#
# Dependencies:
#   python3 + biopython + numpy + pyyaml
#   seqkit          (brew install seqkit)          — sim_01
#   mafft           (brew install mafft)           — sim_02 local mode
#   pbsim3          (conda install -c bioconda pbsim3) — sim_04
#   seqkit          (brew install seqkit)           — sim_01, sim_05

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

LOG_DIR="logs"
mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/sim_pipeline_$(date +%Y%m%d_%H%M%S).log"

# ---- defaults ---------------------------------------------------------------
SOURCE="vimop"
MODE="both"
SEED=42
THREADS=8
DEPTH=30
TARGETS_ARG=()
SKIP_STEPS=()

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
die() { log "ERROR: $*"; exit 1; }

should_run() {   # returns 0 (run) or 1 (skip)
    local step="$1"
    for s in "${SKIP_STEPS[@]:-}"; do [[ "$s" == "$step" ]] && return 1; done
    return 0
}

# ---- parse args -------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --targets) shift; while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                       TARGETS_ARG+=("$1"); shift; done ;;
        --source)  SOURCE="$2";  shift 2 ;;
        --mode)    MODE="$2";    shift 2 ;;
        --seed)    SEED="$2";    shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --depth)   DEPTH="$2";   shift 2 ;;
        --skip)    IFS=',' read -ra SKIP_STEPS <<< "$2"; shift 2 ;;
        *) die "Unknown argument: $1" ;;
    esac
done

# Auto-skip step 1 when using vimop source
if [[ "$SOURCE" == "vimop" ]] && ! printf '%s\n' "${SKIP_STEPS[@]:-}" | grep -qx "1"; then
    SKIP_STEPS+=("1")
    log "source=vimop: step 1 (sim_01) automatically skipped"
fi

TARGETS_FLAG=()
[[ ${#TARGETS_ARG[@]} -gt 0 ]] && TARGETS_FLAG=(--targets "${TARGETS_ARG[@]}")

log "=== sim_pipeline ==="
log "source:  $SOURCE"
log "mode:    $MODE"
log "seed:    $SEED"
log "threads: $THREADS"
log "depth:   $DEPTH"
log "skip:    ${SKIP_STEPS[*]:-none}"
log "log:     $LOG"

# ---- dependency checks ------------------------------------------------------
log "checking core dependencies..."
command -v python3 &>/dev/null || die "python3 not found"
python3 -c "import Bio, numpy, yaml" 2>/dev/null \
    || die "missing python packages (pip install biopython numpy pyyaml)"

# ---- step 1: build multifastas (local mode only) ----------------------------
if should_run 1; then
    log "--- step 1: build multifastas (local) ---"
    python3 "$SCRIPT_DIR/B_sim01_build_multifastas.py" \
        "${TARGETS_FLAG[@]}" --threads "$THREADS" 2>&1 | tee -a "$LOG"
    log "step 1 complete"
else
    log "--- step 1: skipped ---"
fi

# ---- step 2: alignment + weights --------------------------------------------
if should_run 2; then
    log "--- step 2: alignment + weights (source=$SOURCE) ---"
    python3 "$SCRIPT_DIR/B_sim02_align_and_weights.py" \
        "${TARGETS_FLAG[@]}" --source "$SOURCE" --threads "$THREADS" \
        2>&1 | tee -a "$LOG"
    log "step 2 complete"
else
    log "--- step 2: skipped ---"
fi

# ---- step 3: generate mutant genomes ----------------------------------------
if should_run 3; then
    log "--- step 3: generate mutant genomes ---"
    python3 "$SCRIPT_DIR/B_sim03_generate_mutants.py" \
        "${TARGETS_FLAG[@]}" --seed "$SEED" 2>&1 | tee -a "$LOG"
    log "step 3 complete"
else
    log "--- step 3: skipped ---"
fi

# ---- step 4: simulate reads -------------------------------------------------
if should_run 4; then
    log "--- step 4: simulate reads (PBSIM3) ---"
    python3 "$SCRIPT_DIR/B_sim04_simulate_reads.py" \
        "${TARGETS_FLAG[@]}" --seed "$SEED" --depth "$DEPTH" \
        --threads "$THREADS" 2>&1 | tee -a "$LOG"
    log "step 4 complete"
else
    log "--- step 4: skipped ---"
fi

# ---- step 5: mix background -------------------------------------------------
if should_run 5; then
    log "--- step 5: mix background ($MODE) ---"
    python3 "$SCRIPT_DIR/B_sim05_mix_background.py" \
        "${TARGETS_FLAG[@]}" --mode "$MODE" --seed "$SEED" \
        2>&1 | tee -a "$LOG"
    log "step 5 complete"
else
    log "--- step 5: skipped ---"
fi

# ---- summary ----------------------------------------------------------------
DATA="data/petersen_2026"
log "=== simulation complete ==="
log "outputs:"
log "  weights:    $(find $DATA/weights  -name '*.npy' 2>/dev/null | wc -l | tr -d ' ') npy files"
log "  mutants:    $(find $DATA/mutants  -name '*.fasta' 2>/dev/null | wc -l | tr -d ' ') fasta files"
log "  fastq:      $(find $DATA/fastq   -name '*.fastq.gz' 2>/dev/null | wc -l | tr -d ' ') fastq.gz files"
log "  mixed:      $(find $DATA/mixed   -name '*.fastq.gz' 2>/dev/null | wc -l | tr -d ' ') fastq.gz files"
log "  summary:    $DATA/mutants/summary.tsv"
log "  log:        $LOG"
