#!/usr/bin/env bash
# ============================================================================
# 16S — two-pass EM test (Emu reference, both libraries)   [local Mac]
# ============================================================================
# Plan: the project notes. Re-runs only the EM and
# the report (--pipe em,rep) of the two Panel B.2 cells behind manuscript
# Table 1, on their existing annotation (symlinked, not recomputed), with
# the frozen configuration (-N 25, em.temp 10) plus a two-pass EM:
#   em.prune_min_reads 10 — after the first EM, taxa with < 10 expected reads
#   are dropped and the EM is re-run so their reads are re-assigned (Emu's
#   rule; emtrace analysis/em.py::em_loop).
#
#   synthgut   emtrace/emu.N25          -> emtrace/emu.N25.T10.2pass
#   barcode79  emtrace/emu.N25.T10      -> emtrace/emu.N25.T10.2pass
#
# Memory: the EM holds the annotation (0.4 GB synthgut, 1.2 GB barcode79 TSV)
# in memory; cells run one at a time and a guard stops a cell if free memory
# falls below MIN_FREE_PCT (the Mac had two memory panics on 29 Sep 2026).
#
# Then add the manifest rows it prints and score:
#   python3 jobs/framework_16S/framework_16S.score.py --manifest jobs/framework_16S/framework_16S.manifest.tsv
# Run from the repository root:
#   bash jobs/framework_16S/polish.two-pass-em.bash
# ============================================================================

set -uo pipefail
[[ -d jobs ]] || { echo "error: run from the repository root" >&2; exit 1; }

R=results/framework_16S
CELL=emu.N25.T10.2pass
PRUNE="${PRUNE:-10}"
THREADS="${THREADS:-8}"
MIN_FREE_PCT="${MIN_FREE_PCT:-10}"
LOG_DIR="logs/framework_16S"
CFG="$LOG_DIR/configs/emtrace.config.16S.$CELL.yaml"
export EMTRACE_DB_DIR="${EMTRACE_DB_DIR:-$HOME/databases/emtrace}"
mkdir -p "$LOG_DIR/configs"

python3 - jobs/framework_16S/emtrace.config.16S.yaml "$CFG" "$PRUNE" <<'CFGPY' || exit 1
import sys, yaml
src, dst, prune = sys.argv[1], sys.argv[2], float(sys.argv[3])
cfg = yaml.safe_load(open(src))
cfg["annotation"]["minimap"]["max_target_seqs"] = 25   # frozen
cfg["em"]["temp"] = 10                                 # frozen
cfg["em"]["prune_min_reads"] = prune                   # two-pass EM
yaml.safe_dump(cfg, open(dst, "w"), sort_keys=False, default_flow_style=False)
CFGPY

# guarded <label> <cmd...>: stop cmd if free memory < MIN_FREE_PCT
guarded() {
  local label="$1"; shift
  set -m; "$@" & local pid=$! pgid; pgid=$(ps -o pgid= -p "$pid" | tr -d ' '); set +m
  while kill -0 "$pid" 2>/dev/null; do
    local free; free=$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub("%","",$2); print $2}')
    if [[ -n "$free" && "$free" -lt "$MIN_FREE_PCT" ]]; then
      echo "[GUARD] $label: free memory ${free}% < ${MIN_FREE_PCT}% — stopping" >&2
      kill -TERM -- "-$pgid" 2>/dev/null; sleep 5; kill -KILL -- "-$pgid" 2>/dev/null
      wait "$pid" 2>/dev/null; return 143
    fi
    sleep 5
  done
  wait "$pid"
}

# run_cell <library> <reads> <source annotation dir>
run_cell() {
  local lib="$1" reads="$2" src="$3"
  local out="$R/$lib/emtrace/$CELL"
  local abund="$out/$lib/results/reports/abund_table.tsv"
  if [[ -f "$abund" ]]; then echo "[skip] $lib/$CELL exists"; return 0; fi
  [[ -f "$src/annot1.tsv" ]] || { echo "[FAIL] no annotation in $src" >&2; return 1; }
  mkdir -p "$out/$lib"
  ln -sfn "$PWD/$src" "$out/$lib/annotation"
  echo "[run ] $lib/$CELL (em.prune_min_reads $PRUNE)"
  guarded "$lib" emtrace run -1 "$reads" -s "$lib" -r "$out" -t "$THREADS" \
      --config "$CFG" --pipe em,rep --logs "$LOG_DIR/$CELL.$lib" -v \
      > "$LOG_DIR/$CELL.$lib.log" 2>&1
  local rc=$?
  [[ $rc -eq 0 && -f "$abund" ]] || { echo "[FAIL] $lib exit $rc (143 = memory guard) — $LOG_DIR/$CELL.$lib.log" >&2; return 1; }
  grep -h "EM second pass" "$LOG_DIR/$CELL.$lib.log" | sed 's/.*| //'
  printf 'Manifest row:\nEMtrace.emu.2pass\t%s\tdiag\temtrace\tncbi\t%s\ttwo-pass EM (em.prune_min_reads %s), otherwise as EMtrace.emu (the project notes)\n' \
    "$lib" "$abund" "$PRUNE"
}

run_cell synthgut  data/sra/16S/SRR14326315.fastq          "$R/synthgut/emtrace/emu.N25/annotation"
run_cell barcode79 data/mock_community.16S/barcode79.fq.gz "$R/barcode79/emtrace/emu.N25.T10/barcode79/annotation"
