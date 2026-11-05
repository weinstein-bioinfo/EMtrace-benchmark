#!/usr/bin/env bash
# ============================================================================
# 16S — matched-read re-run (synthetic gut, Emu reference)   [local Mac]
# ============================================================================
# Plan: the project notes. On the synthetic gut
# (SRR14326315, 521,661 reads) EMtrace's Panel B.2 cell `emu.N25` classified
# the 324,400 reads its preprocessing kept (porechop_abi + chopper Q10,
# quality end-trimming, 1,200-1,800 bp), while Emu classified all 521,661. The
# species-L1 gap (Emu 0.501, EMtrace 0.527) therefore mixes the method with
# the read input. Two runs complete a 2 x 2 (tool x read set):
#
#                          all reads (Emu's input)   EMtrace's filtered reads
#   EMtrace (emu ref)   PART emtrace  (new)          emu.N25        (existing)
#   Emu (own DB)           emu/       (existing)     PART emu       (new)
#
# PART emu — Emu at its defaults (as framework_16S.comparators.bash), own DB,
#   on EMtrace's filtered reads of emu.N25 (reads/reads.chopped.fq.gz,
#   324,400 reads). --keep-read-assignments is added for the planned bootstrap
#   (the project notes); it does not change the abundances. Emu peaked at
#   5.5 GB on all 521,661 reads.
#   Output: results/framework_16S/synthgut/emu.emtracereads/
#
# PART emtrace — the frozen Panel B.2 cell (--anns 16S.emu, -N 25, em.temp 10)
#   with chopper's QUALITY filtering switched off: --quality 0 and no quality
#   end-trimming (trim_approach '', cutoff 0). The 1,200-1,800 bp length window
#   and porechop_abi adapter trimming are kept.
#   Reads: the Mac copy data/sra/16S/SRR14326315.fastq. It is an SRA Lite
#   copy with synthetic quality strings; with the quality filter off nothing downstream reads quality
#   (porechop_abi, chopper's length filter, minimap2 and the EM do not), and
#   its sequences are byte-identical to the Roihu copy.
#   Output: results/framework_16S/synthgut/emtrace/emu.N25.T10.noQ/
#
# MEMORY (16 GB Mac). emu.N25 peaked at 15.7 GB on Roihu with 324,400 reads;
# without the quality filter EMtrace sees more reads, so it may not fit.
# The parts therefore run one after the other, never together, and a guard
# stops a part (SIGTERM to its process group) when free memory falls below
# MIN_FREE_PCT — a stopped run is better than a frozen machine (the Mac had
# two kernel watchdog panics from memory exhaustion on 29 Sep 2026). Close
# other large applications first. If the guard trips on PART emtrace, run it on
# Roihu instead (framework_16S.emtrace.bash conventions).
#
# Usage (from the repo root):
#   bash jobs/framework_16S/polish.match-read-rerun.bash              # emu, then emtrace
#   bash jobs/framework_16S/polish.match-read-rerun.bash --part emu
#   bash jobs/framework_16S/polish.match-read-rerun.bash --part emtrace
# Then add the printed manifest rows to framework_16S.manifest.tsv and score:
#   python3 jobs/framework_16S/framework_16S.score.py --manifest jobs/framework_16S/framework_16S.manifest.tsv
# ============================================================================

set -uo pipefail
[[ -d jobs ]] || { echo "error: run from the repository root" >&2; exit 1; }

PART=all
while [[ $# -gt 0 ]]; do
  case "$1" in
    --part) PART="$2"; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ "$PART" == all || "$PART" == emu || "$PART" == emtrace ]] || { echo "usage: $0 [--part emu|emtrace]" >&2; exit 1; }

SAMPLE=synthgut
CELL=emu.N25.T10.noQ
THREADS="${THREADS:-8}"
MIN_FREE_PCT="${MIN_FREE_PCT:-10}"   # memory guard threshold (macOS memory_pressure)
LOG_DIR="logs/framework_16S"
TIMINGS="metrics/framework_16S/timings/framework_16S.comparators.timings.csv"
export EMTRACE_DB_DIR="${EMTRACE_DB_DIR:-$HOME/databases/emtrace}"
mkdir -p "$LOG_DIR" "$(dirname "$TIMINGS")"

# guarded <label> <cmd...>: run cmd in its own process group; stop it if free
# memory falls below MIN_FREE_PCT. Returns the command's exit status (143 if stopped).
guarded() {
  local label="$1"; shift
  set -m                                   # own process group for the job
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

# ============================================================================
# PART emu — Emu on EMtrace's filtered reads
# ============================================================================
run_emu() {
  local EMU_PY="${EMU_PY:-data/comparators/emu-venv/bin/python}"
  local EMU_SCRIPT="${EMU_SCRIPT:-data/comparators/emu-src/emu}"
  local EMU_DB="${EMU_DB:-data/comparators/emu_db}"
  local READS="results/framework_16S/$SAMPLE/emtrace/emu.N25/reads/reads.chopped.fq.gz"
  local OUT="results/framework_16S/$SAMPLE/emu.emtracereads"

  [[ -x "$EMU_PY" && -f "$EMU_SCRIPT" ]] || { echo "Emu not found ($EMU_PY, $EMU_SCRIPT)" >&2; return 1; }
  [[ -f "$EMU_DB/species_taxid.fasta" && -f "$EMU_DB/taxonomy.tsv" ]] || { echo "Emu DB incomplete: $EMU_DB" >&2; return 1; }
  [[ -s "$READS" ]] || { echo "EMtrace filtered reads not found: $READS" >&2; return 1; }
  command -v minimap2 > /dev/null || { echo "minimap2 not on PATH (Emu shells out to it)" >&2; return 1; }
  mkdir -p "$OUT"

  if [[ -f "$OUT/${SAMPLE}_rel-abundance.tsv" ]]; then
    echo "[skip] Emu output exists ($OUT; remove it to redo)"
  else
    echo "=== PART emu: Emu on EMtrace's filtered reads ($READS) ==="
    guarded "emu" python3 jobs/bench_time.py --csv "$TIMINGS" --label "emu.emtracereads.$SAMPLE" -- \
      "$EMU_PY" "$EMU_SCRIPT" abundance \
        --db "$EMU_DB" \
        --type map-ont \
        --output-basename "$SAMPLE" \
        --threads "$THREADS" \
        --keep-counts \
        --keep-read-assignments \
        --output-dir "$OUT" \
        "$READS" \
      > "$LOG_DIR/emu.emtracereads.$SAMPLE.log" 2>&1
    local rc=$?
    [[ $rc -eq 0 && -f "$OUT/${SAMPLE}_rel-abundance.tsv" ]] || {
      echo "[FAIL] Emu exit $rc (143 = stopped by the memory guard) — see $LOG_DIR/emu.emtracereads.$SAMPLE.log" >&2; return 1; }
  fi

  {
    echo "date:        $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "emu:         $("$EMU_PY" "$EMU_SCRIPT" --version 2>&1 | head -1)"
    echo "database:    $EMU_DB"
    echo "reads:       $READS  (EMtrace emu.N25 preprocessing: porechop_abi + chopper Q10, 1200-1800 bp)"
    echo "n_reads:     $(gzip -dc "$READS" | awk 'END{print NR/4}')"
    echo "reads_md5:   $(md5 -q "$READS")"
    echo "threads:     $THREADS"
    echo "thresholds:  Emu defaults; --keep-read-assignments for the bootstrap"
  } > "$OUT/emu.provenance.txt"

  echo "Manifest row:"
  printf 'Emu.emtracereads\t%s\tB.2\temu\tncbi\t%s/%s_rel-abundance.tsv\tMATCHED READS: Emu defaults, own DB, on EMtrace.emu filtered reads (324,400)\n' \
    "$SAMPLE" "$OUT" "$SAMPLE"
}

# ============================================================================
# PART emtrace — EMtrace on all reads, quality filtering off
# ============================================================================
run_emtrace() {
  local READS="data/sra/16S/SRR14326315.fastq"
  local OUTDIR="results/framework_16S/$SAMPLE/emtrace/$CELL"
  local CFG_DIR="$LOG_DIR/configs"
  local BASE_CFG="jobs/framework_16S/emtrace.config.16S.yaml"
  local CFG="$CFG_DIR/emtrace.config.16S.$CELL.yaml"
  local ABUND="$OUTDIR/$SAMPLE/results/reports/abund_table.tsv"

  local tool
  for tool in emtrace minimap2 porechop_abi chopper; do
    command -v "$tool" > /dev/null || { echo "Missing tool: $tool" >&2; return 1; }
  done
  [[ -f "$READS" ]] || { echo "reads not found: $READS" >&2; return 1; }
  [[ -f "$EMTRACE_DB_DIR/minimap.emu.fasta.mmi" ]] || {
    echo "Emu-reference index not found: $EMTRACE_DB_DIR/minimap.emu.fasta.mmi" >&2; return 1; }
  mkdir -p "$OUTDIR" "$CFG_DIR"

  if [[ -f "$ABUND" ]]; then
    echo "[skip] $CELL (abund_table.tsv exists; rm -rf $OUTDIR to redo)"
  else
    # frozen Panel B.2 configuration, quality filtering off, length window kept
    python3 - "$BASE_CFG" "$CFG" <<'CFGPY' || { echo "[FATAL] config derivation failed" >&2; return 1; }
import sys, yaml
src, dst = sys.argv[1:]
cfg = yaml.safe_load(open(src))
cfg["annotation"]["minimap"]["max_target_seqs"] = 25   # -N 25 (frozen)
cfg["em"]["temp"] = 10                                 # em.temp 10 (frozen)
ch = cfg["preprocess"]["chopper"]
ch["quality"] = 0          # no mean-quality filter
ch["trim_approach"] = ""   # no quality end-trimming (flag omitted)
ch["cutoff"] = 0           # (flag omitted)
assert ch["minlength"] == 1200 and ch["maxlength"] == 1800, "length window changed"
yaml.safe_dump(cfg, open(dst, "w"), sort_keys=False, default_flow_style=False)
CFGPY
    echo "=== PART emtrace: $CELL  --anns 16S.emu  config=$(basename "$CFG") ==="
    guarded "emtrace" python3 jobs/bench_time.py --csv "$TIMINGS" --label "$SAMPLE.$CELL" -- \
      emtrace run \
        -1 "$READS" \
        -s "$SAMPLE" \
        -r "$OUTDIR" \
        -t "$THREADS" \
        --config "$CFG" \
        --anns 16S.emu \
        --pipe pore,chop,ann,em,rep,plot \
        --logs "$LOG_DIR/$CELL" \
        -v \
      > "$LOG_DIR/$CELL.log" 2>&1
    local rc=$?
    [[ $rc -eq 0 && -f "$ABUND" ]] || {
      echo "[FAIL] EMtrace exit $rc (143 = stopped by the memory guard) — see $LOG_DIR/$CELL.log" >&2; return 1; }
  fi

  echo "reads after chopper: $(grep -m1 "chopper: kept" "$LOG_DIR/$CELL.log" | sed 's/.*chopper: //')"
  echo "Manifest row:"
  printf 'EMtrace.emu.noQ\t%s\tB.2\temtrace\tncbi\t%s\tMATCHED READS: as EMtrace.emu, chopper quality filter and quality trimming off, length 1200-1800 kept\n' \
    "$SAMPLE" "$ABUND"
}

echo "memory guard: stop a part if free memory < ${MIN_FREE_PCT}%   (now: $(memory_pressure | awk -F': ' '/free percentage/{print $2}'))"
[[ "$PART" == all || "$PART" == emu  ]] && { run_emu  || exit 1; }
[[ "$PART" == all || "$PART" == emtrace ]] && { run_emtrace || exit 1; }
