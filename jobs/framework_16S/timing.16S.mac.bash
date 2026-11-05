#!/usr/bin/env bash
# ============================================================================
# 16S — matched-hardware timing, EMtrace on the Mac   [local Mac, 8 threads]
# ============================================================================
# The 16S cells timed EMtrace on Roihu (32 cores) and Emu /
# Spaghetti on the Mac (8 threads), so their costs are not comparable. This
# script times EMtrace on the Mac, same threads, same minimap2, same
# reference (Emu's own sequences), one run at a time. Frozen Panel B.2 cell
# (--anns 16S.emu, -N 25, em.temp 10); the Roihu cells are untouched:
#
#   PART emtrace_b79    barcode79 (100,128 reads, real qualities)
#                    -> barcode79/emtrace/emu.N25.T10.mac
#   PART emtrace_sg     synthgut, SRR14326315 (521,661 reads). The Mac copy is
#                    SRA Lite with synthetic qualities, so chopper's quality
#                    filter and quality trimming are off (suffix noQ = "no
#                    quality filtering", as emu.N25.T10.noQ); the
#                    1,200-1,800 bp length window is kept.
#                    -> synthgut/emtrace/emu.N25.T10.noQ.mac
#   PART emtrace_gutmock ZymoBIOMICS Gut Microbiome Standard, R10.4.1,
#                    gutmock_sup_r1..r4 (real qualities)
#                    -> gutmock_sup_rN/emtrace/emu.N25.T10.mac
#
# Timings: metrics/framework_16S/timings/framework_16S.mac.timings.csv
# Memory guard as in polish.match-read-rerun.bash (16 GB Mac).
#
# Usage (from the repo root). caffeinate -is keeps the Mac awake (no idle or
# system sleep) for as long as the script runs:
#   caffeinate -is bash jobs/framework_16S/timing.16S.mac.bash                      # parts whose RUN_* flag is true
#   caffeinate -is bash jobs/framework_16S/timing.16S.mac.bash --part emtrace_gutmock  # that part only
# Run flags (edit below, or override from the environment, e.g. RUN_EMTRACE_B79=true):
#   RUN_EMTRACE_B79=false (done), RUN_EMTRACE_SG=true, RUN_EMTRACE_GUTMOCK=true
# A part whose abund_table.tsv exists is skipped anyway.
# Detached and awake (survives the terminal / agent session, no sleep), at
# normal priority — launched through bash because zsh renices background jobs
# (+5, option BG_NICE), and an unprivileged process cannot undo that:
#   bash -c 'nohup caffeinate -is bash jobs/framework_16S/timing.16S.mac.bash  > logs/framework_16S/timing.16S.mac.driver.log 2>&1 < /dev/null &'
#   bash -c 'nohup caffeinate -is bash jobs/framework_16S/timing.16S.mac.bash --part emtrace_gutmock  > logs/framework_16S/timing.16S.mac.driver.log 2>&1 < /dev/null &'
# ============================================================================

set -uo pipefail
[[ -d jobs ]] || { echo "error: run from the repository root" >&2; exit 1; }
# Timings must be taken at normal priority (see Usage)
NICE=$(ps -o nice= -p $$ | tr -d ' ')
[[ "$NICE" == 0 ]] || { echo "error: running at nice $NICE; launch as in Usage (bash -c '...')" >&2; exit 1; }
# Matched-hardware timing: refuse to run beside another copy
LOCK="logs/framework_16S/.timing.16S.mac.lock"
mkdir -p "$(dirname "$LOCK")"
if ! mkdir "$LOCK" 2>/dev/null; then
  if kill -0 "$(cat "$LOCK/pid" 2>/dev/null)" 2>/dev/null; then
    echo "error: another timing.16S.mac.bash is running (pid $(cat "$LOCK/pid"))" >&2; exit 1
  fi
  echo "removing stale lock $LOCK" >&2
fi
echo $$ > "$LOCK/pid"
trap 'rm -rf "$LOCK"' EXIT

RUN_EMTRACE_B79="${RUN_EMTRACE_B79:-false}"
RUN_EMTRACE_SG="${RUN_EMTRACE_SG:-true}"
RUN_EMTRACE_GUTMOCK="${RUN_EMTRACE_GUTMOCK:-true}"

PART=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --part) PART="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 1 ;;
  esac
done
case "$PART" in
  "") ;;
  emtrace_b79)     RUN_EMTRACE_B79=true;  RUN_EMTRACE_SG=false; RUN_EMTRACE_GUTMOCK=false ;;
  emtrace_sg)      RUN_EMTRACE_B79=false; RUN_EMTRACE_SG=true;  RUN_EMTRACE_GUTMOCK=false ;;
  emtrace_gutmock) RUN_EMTRACE_B79=false; RUN_EMTRACE_SG=false; RUN_EMTRACE_GUTMOCK=true  ;;
  *) echo "unknown part: $PART" >&2; exit 1 ;;
esac

THREADS="${THREADS:-8}"
MIN_FREE_PCT="${MIN_FREE_PCT:-10}"
LOG_DIR="logs/framework_16S"
TIMINGS="metrics/framework_16S/timings/framework_16S.mac.timings.csv"
export EMTRACE_DB_DIR="${EMTRACE_DB_DIR:-$HOME/databases/emtrace}"
mkdir -p "$LOG_DIR" "$(dirname "$TIMINGS")"

declare -A READS=(
  ["barcode79"]="data/mock_community.16S/barcode79.fq.gz"
  ["synthgut"]="data/sra/16S/SRR14326315.fastq"
  ["gutmock_sup_r1"]="data/sra/16S/gutmock_sup_r1.fastq.gz"
  ["gutmock_sup_r2"]="data/sra/16S/gutmock_sup_r2.fastq.gz"
  ["gutmock_sup_r3"]="data/sra/16S/gutmock_sup_r3.fastq.gz"
  ["gutmock_sup_r4"]="data/sra/16S/gutmock_sup_r4.fastq.gz"
)

# guarded <label> <cmd...>: run cmd in its own process group; stop it if free
# memory falls below MIN_FREE_PCT. Returns the command's exit status (143 if stopped).
guarded() {
  local label="$1"; shift
  set -m
  "$@" &
  local pid=$! pgid; pgid=$(ps -o pgid= -p "$pid" | tr -d ' ')
  set +m
  while kill -0 "$pid" 2>/dev/null; do
    local free
    free=$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub("%","",$2); print $2}')
    if [[ -n "$free" && "$free" -lt "$MIN_FREE_PCT" ]]; then
      echo "[GUARD] $label: free memory ${free}% < ${MIN_FREE_PCT}% — stopping it" >&2
      kill -TERM -- "-$pgid" 2>/dev/null
      sleep 5; kill -KILL -- "-$pgid" 2>/dev/null
      wait "$pid" 2>/dev/null
      return 143
    fi
    sleep 5
  done
  wait "$pid"
}

# run_emtrace <library> [noQ]: the frozen Panel B.2 cell; "noQ" switches chopper's
# quality filter and quality trimming off (synthetic qualities)
run_emtrace() {
  local lib="$1" noq="${2:-}" reads="${READS[$1]}" cell="emu.N25.T10${2:+.$2}.mac"
  local outdir="results/framework_16S/$lib/emtrace/$cell"
  local cfg="$LOG_DIR/configs/emtrace.config.16S.$cell.yaml"
  local abund="$outdir/$lib/results/reports/abund_table.tsv"
  [[ -f "$EMTRACE_DB_DIR/minimap.emu.fasta.mmi" ]] || {
    echo "Emu-reference index not found in $EMTRACE_DB_DIR" >&2; return 1; }
  [[ -s "$reads" ]] || { echo "reads not found: $reads" >&2; return 1; }
  if [[ -f "$abund" ]]; then
    echo "[skip] EMtrace $lib $cell (rm -rf $outdir to redo)"; return 0
  fi
  mkdir -p "$outdir" "$(dirname "$cfg")"
  python3 - jobs/framework_16S/emtrace.config.16S.yaml "$cfg" "$noq" <<'CFGPY' || return 1
import sys, yaml
src, dst, noq = sys.argv[1:]
cfg = yaml.safe_load(open(src))
cfg["annotation"]["minimap"]["max_target_seqs"] = 25   # -N 25 (frozen)
cfg["em"]["temp"] = 10                                 # em.temp 10 (frozen)
if noq:
    ch = cfg["preprocess"]["chopper"]
    ch["quality"] = 0          # no mean-quality filter
    ch["trim_approach"] = ""   # no quality end-trimming (flag omitted)
    ch["cutoff"] = 0           # (flag omitted)
    assert ch["minlength"] == 1200 and ch["maxlength"] == 1800, "length window changed"
yaml.safe_dump(cfg, open(dst, "w"), sort_keys=False, default_flow_style=False)
CFGPY
  echo "=== EMtrace $lib $cell  ($(date +%H:%M)) ==="
  guarded "emtrace.$lib" python3 jobs/bench_time.py --csv "$TIMINGS" --label "emtrace.$lib.$cell" -- \
    emtrace run -1 "$reads" -s "$lib" -r "$outdir" -t "$THREADS" \
      --config "$cfg" --anns 16S.emu --pipe pore,chop,ann,em,rep,plot \
      --logs "$LOG_DIR/$cell.$lib" -v \
    > "$LOG_DIR/$cell.$lib.log" 2>&1
  local rc=$?
  [[ $rc -eq 0 && -f "$abund" ]] || {
    echo "[FAIL] EMtrace $lib exit $rc (143 = memory guard) — see $LOG_DIR/$cell.$lib.log" >&2; return 1; }
}

echo "minimap2 $(minimap2 --version); threads $THREADS; free memory now: $(memory_pressure | awk -F': ' '/free percentage/{print $2}')"
echo "parts: emtrace_b79=$RUN_EMTRACE_B79 emtrace_sg=$RUN_EMTRACE_SG emtrace_gutmock=$RUN_EMTRACE_GUTMOCK"
[[ "$RUN_EMTRACE_B79" == true ]] && { run_emtrace barcode79 || exit 1; }
[[ "$RUN_EMTRACE_SG"  == true ]] && { run_emtrace synthgut noQ || exit 1; }
if [[ "$RUN_EMTRACE_GUTMOCK" == true ]]; then
  for r in 1 2 3 4; do run_emtrace "gutmock_sup_r$r" || exit 1; done
fi
echo "done $(date +%H:%M); timings in $TIMINGS"
