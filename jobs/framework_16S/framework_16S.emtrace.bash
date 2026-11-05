#!/bin/bash -l
#SBATCH --job-name=framework_16S
#SBATCH --account=project_XXXXXXX
#SBATCH --partition=interactive
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --chdir=/projappl/project_XXXXXXX/EMtrace-benchmark
#SBATCH --output=/scratch/project_XXXXXXX/emtrace-benchmark/logs/framework_16S.%j.out
# ============================================================================
# Framework 16S — EMtrace cells: every library x every configuration (Roihu)
# ============================================================================
#
# The one EMtrace runner for the 16S benchmark. Writes one result tree per
# cell for jobs/framework_16S/framework_16S.score.py. Design: the project notes
#
# LIBRARIES (switches below)
#   barcode79  ZymoBIOMICS mock, 100 k reads   -- the Panel A tuning library
#   synthgut   synthetic gut mock, 522 k reads -- the validation library
#
# CELL GROUPS. Each runs for every selected library. A cell whose
# abund_table.tsv exists is skipped, so leaving groups on only fills gaps.
#   RUN_FROZEN    <reference>.N25.T10 -- the configuration Panel A settled on,
#                 run unchanged on each selected reference. The `emu` cell is
#                 Panel B.2: aligning against Emu's own sequences is what makes
#                 a EMtrace-vs-Emu difference a method difference rather
#                 than a database one (roadmap section 7.2).
#   RUN_MATRIX    Panel A axis 1 -- one cell per reference at the shipped
#                 defaults. Axes 3-4 are recomputed by the scorer, free.
#   RUN_TOPHITS   axis 9 -- --tophits variants; re-annotates, because the
#                 filter acts in the annotate step.
#   RUN_SWEEPS    axis 6 on one base cell (rrnDB.dedup): -N 0/25/50, 3 annotations
#   RUN_TUNED     -N 25 x em.temp 50/100, em,rep over the -N 25 annotation.
#
# Cell directories are named for their configuration, not for the run that
# produced them: <reference>[.N<n>.T<temp>][.modifier].
#
# --pipe pore,chop,ann,em,rep,plot -- `flt` omitted: no host in a mock.
# --config jobs/framework_16S/emtrace.config.16S.yaml -- chopper 1200-1800,
#   minimap s=200, -N 5, report.min_qcov 0.40, em.temp 10.0.
#
# Skip logic: a cell is skipped when its abund_table.tsv exists, so an
# interrupted run resumes. To force one:  rm -rf $RES_ROOT/<cell>
#
# Layout -- EVERYTHING on /scratch, nothing written into the repo, so one rsync
# of $SCRATCH_BASE pulls a whole run to the Mac and on to Allas:
#   $SCRATCH_BASE = /scratch/project_XXXXXXX/emtrace-benchmark
#     results/framework_16S/<library>/emtrace/<cell>/<library>/results/
#     logs/framework_16S.<jobid>.out          Slurm stdout
#     logs/framework_16S/<library>/           bash + pipeline logs
#     logs/framework_16S/configs/             derived configs, provenance
#     metrics/framework_16S/timings/          one row per cell
#
# Create $SCRATCH_BASE/logs BEFORE the first submission: Slurm opens --output
# before the script runs and will not create the directory itself.
#
# Submit:
#   sbatch jobs/framework_16S/framework_16S.emtrace.bash            # everything
#   RUN_SYNTHGUT=false RUN_SWEEPS=false sbatch ...             # one library
#
# Then score and plot (one manifest covers both libraries):
#   python3 jobs/framework_16S/framework_16S.score.py --manifest jobs/framework_16S/framework_16S.manifest.tsv
#   Rscript rplots/mock16S_composition.R
#   LIBRARY=synthgut Rscript rplots/mock16S_composition.R
#
# Resources: 32 cores, 60 GB, 36 h. 60 GB is set by the SILVA cell, which
# indexes 14 GB of FASTA in memory; every other cell peaks under 16 GB, and
# memory does not scale with read count (15.7 GB on a 5.2x larger library).
# Measured on Roihu, 32 cores, 100 k reads: rrnDB 1 h 46 m | rrnDB.dedup 27 m |
# refseq 8 m | chain1 1 h 41 m | silva 1 h 05 m | em,rep re-runs ~1.5 min.
# A full Panel A pass is ~8 h; the frozen cells alone are well under an hour
# per library.
# ============================================================================

# =============================================================================
# Environment on Roihu
# =============================================================================
# Load minimap2, samtools and seqkit. Do not run `module purge` on Roihu,
# it removes the module tree these tools come from.
module load bio-apps minimap2 samtools seqkit

# chopper is not available as a module, so add its folder to PATH.
export PATH="/projappl/project_XXXXXXX/bin:/projappl/project_XXXXXXX/install/chopper_dir:$PATH"

# Activate EMtrace (also provides porechop_abi).
source /projappl/project_XXXXXXX/venv/emtrace/bin/activate

# The modules above set PYTHONPATH to a different numpy, which breaks EMtrace.
# Clear it, and do it after the activate line.
unset PYTHONPATH

# Where EMtrace finds its databases.
export EMTRACE_DB_DIR=/scratch/project_XXXXXXX/databases/emtrace

set -uo pipefail

# Stop early if a tool is missing, rather than hours into the run.
for tool in emtrace minimap2 porechop_abi chopper; do
  if ! command -v "$tool" > /dev/null; then
    echo "Missing tool: $tool" >&2
    exit 1
  fi
done

# Stop early if EMtrace predates the changes this run depends on.
if ! emtrace run --help | grep -q -- "--tophits"; then
  echo "EMtrace is too old for run v3 (no --tophits). Update the venv first." >&2
  exit 1
fi

echo "Job ${SLURM_JOB_ID:-local} on $(hostname), started $(date)"
echo

BASE_CFG="jobs/framework_16S/emtrace.config.16S.yaml"
SRA_DIR="/scratch/project_XXXXXXX/emtrace-benchmark/sra"

# Reads per library. Adding a library is one line here plus one switch below.
# Accession or file name per library, when it is not the library name itself.
declare -A LIB_FILE=(
  ["synthgut"]="$SRA_DIR/SRR14326315.fastq"
)

# Resolve whichever extension is actually on disk: barcode79 is .fq.gz,
# synthgut is a plain .fastq. Hard-coding one suffix silently skips a library.
reads_for() {
  local lib="$1" c
  for c in "$SRA_DIR/$lib.fq.gz" "$SRA_DIR/$lib.fastq.gz" \
           "$SRA_DIR/$lib.fq" "$SRA_DIR/$lib.fastq" \
           "${LIB_FILE[$lib]:-}"; do
    [[ -n "$c" && -f "$c" ]] && { echo "$c"; return 0; }
  done
  echo ""
}
# Results live on /scratch, not under the repo on /projappl (small quota; the
# previous run nearly filled it). Logs, configs and timings stay in the repo.
# Everything this job writes goes to /scratch, logs and metrics included, so a
# run is pulled to the Mac (and on to Allas) with one rsync of one directory.
# NOTE: Slurm opens --output before the script runs and does NOT create its
# directory, so SCRATCH_BASE/logs must already exist at submit time.
SCRATCH_BASE="/scratch/project_XXXXXXX/emtrace-benchmark"
RESULTS_BASE="$SCRATCH_BASE/results"
LOG_BASE="$SCRATCH_BASE/logs/framework_16S"
CFG_DIR="$LOG_BASE/configs"         # derived config variants, kept for provenance
METRICS_DIR="$SCRATCH_BASE/metrics/framework_16S/timings"
TIMING_CSV="$METRICS_DIR/framework_16S.timings.csv"
# RES_ROOT and LOG_DIR are per-library and set by the loop at the end.
THREADS=${SLURM_CPUS_PER_TASK:-8}   # follows --cpus-per-task

mkdir -p "$RES_ROOT" "$LOG_DIR" "$CFG_DIR" "$METRICS_DIR"

# Refuse to write results onto /projappl, however the path was reached.
case "$(cd "$RES_ROOT" && pwd -P)" in
  /projappl/*) echo "[FATAL] RES_ROOT resolves under /projappl: $RES_ROOT" >&2; exit 1 ;;
esac

# =============================================================================
# Control switches — set to true/false before running
# =============================================================================
# CURRENT STATE (22 Sep 2026): set to the ONE cell still missing.
#
#   barcode79 x emu  ->  emu.N25.T10     ~25 min, the Panel B.2 cell for the
#                                        tuning library
#
# Everything else in this script has already run and is on disk:
#   barcode79  24 Panel A cells (matrix, tophits, -N sweep, tuned)
#   synthgut   rrnDB.dedup.N25, emu.N25, refseq.N25, silva.N25
#
# The skip-on-abund_table.tsv check would make a full run safe anyway, but two
# things argue for narrowing it here. The barcode79 Panel A cells carry their
# OLD names (rrnDB.dedup, sweep.secondary.N25), so RUN_FROZEN + RUN_RRNDB_DEDUP
# would not recognise them and would re-run ~30 min of annotation to produce
# rrnDB.dedup.N25.T10 -- the same configuration under a new name. And a full
# pass re-derives every sweep config, which is noise in the log when one cell
# is wanted.
#
# TO RESTORE A FULL PASS: set every switch in this block back to true. Each one
# also honours an environment override, so a one-off needs no edit:
#   RUN_SYNTHGUT=true RUN_MATRIX=true sbatch <this script>

# WHICH LIBRARIES
RUN_BARCODE79="${RUN_BARCODE79:-true}"  # ZymoBIOMICS mock, 100 k reads (Panel A tuning library)
RUN_SYNTHGUT="${RUN_SYNTHGUT:-false}"   # synthetic gut mock, 522 k reads -- all cells present

# WHICH CELL GROUPS. Each runs for every selected library; a cell whose
# abund_table.tsv exists is skipped, so leaving them on only fills gaps.
RUN_FROZEN="${RUN_FROZEN:-true}"    # <ref>.N25.T10 per reference -- the only group with work left
RUN_MATRIX="${RUN_MATRIX:-false}"   # Panel A axis 1: done for barcode79
RUN_TOPHITS="${RUN_TOPHITS:-false}" # axis 9, --tophits: done
RUN_SWEEPS="${RUN_SWEEPS:-false}"   # axis 6: done
RUN_TUNED="${RUN_TUNED:-false}"     # -N 25 x em.temp 50/100: done

# Per-reference switches. Only emu has a missing cell.
RUN_EMU="${RUN_EMU:-true}"              # Emu's own reference, repackaged unchanged.
                       # Panel B.2: aligning against the exact sequences Emu
                       # does is what makes a EMtrace-vs-Emu difference a
                       # method difference (roadmap section 7.2).
RUN_RRNDB="${RUN_RRNDB:-false}"          # rrnDB 5.10, all 265 k sequences
RUN_RRNDB_DEDUP="${RUN_RRNDB_DEDUP:-false}"  # rrnDB 5.10 de-duplicated, 79 k seqs
RUN_REFSEQ16S="${RUN_REFSEQ16S:-false}"  # RefSeq 16S alone
RUN_CHAIN="${RUN_CHAIN:-false}"          # rrnDB -> RefSeq chain. On this mock the RefSeq
                       # pass got 4 unassigned reads and added nothing, so the
                       # cell equals rrnDB; ~1 h 45 m if switched on.
RUN_SILVA="${RUN_SILVA:-false}"          # SILVA 144 (reports through SILVA taxonomy)

# Cells the axis-9 (--tophits) variants are derived from; one annotation each.
# rrnDB.dedup is cheap. Adding "rrnDB" (~1 h 45 m) reproduces the earlier Roihu
# top-hit baseline under the new tie rule.
TOPHITS_CELLS=("rrnDB.dedup")

# Which sweeps to run in Section B
SWEEP_SECONDARY=true   # axis 6 — minimap2 -N: 0 / 25 / 50 (5 = base cell)

# em.temp values crossed with -N 25 in the tuned cells
TUNED_TEMPS=(50.0 100.0)

# Base of the Section B sweeps. rrnDB.dedup: ~3.6x cheaper to annotate than
# rrnDB, and de-duplication is what lets -N reach genuinely different
# alternatives instead of identical copies (axis 6 is only meaningful on it).
# Axes 6-8 are NOT crossed with axis 1.
# The frozen configuration (roadmap section 6.2), used by RUN_FROZEN.
FROZEN_N=25
FROZEN_TEMP=10

SWEEP_BASE_ANNS="16S.rrnDB.dedup"
SWEEP_BASE_CELL="rrnDB.dedup"   # Section A cell produced by SWEEP_BASE_ANNS
# =============================================================================

# anns value of each Section A cell, for the axis-9 re-annotations
declare -A CELL_ANNS=(
  ["rrnDB"]="16S.rrnDB"
  ["rrnDB.dedup"]="16S.rrnDB.dedup"
  ["refseq"]="16S.refseq"
  ["chain1"]="16S.chain1"
  ["silva"]="16S.silva"
  ["emu"]="16S.emu"
)


# --- Derive a config variant -------------------------------------------------
# Takes dotted key=value overrides and writes a new YAML under $CFG_DIR.
# Overrides go through a YAML parse rather than sed because keys such as
# `max_target_seqs` and `s` are not unique as plain text.
mkcfg() {
  local name="$1"; shift
  local out="$CFG_DIR/emtrace.config.16S.${name}.yaml"
  python3 - "$BASE_CFG" "$out" "$@" <<'CFGPY'
import sys, yaml
src, dst, *overrides = sys.argv[1:]
cfg = yaml.safe_load(open(src))
for ov in overrides:
    key, _, raw = ov.partition("=")
    try:
        val = yaml.safe_load(raw)
    except Exception:
        val = raw
    node = cfg
    parts = key.split(".")
    for p in parts[:-1]:
        if p not in node:
            raise SystemExit(f"config path not found: {key} (missing section {p!r})")
        node = node[p]
    if parts[-1] not in node:
        raise SystemExit(f"config key not found: {key}")
    node[parts[-1]] = val
yaml.safe_dump(cfg, open(dst, "w"), sort_keys=False, default_flow_style=False)
CFGPY
  local status=$?
  # A typo in an override key must stop the run rather than silently falling
  # through to a stale or absent config.
  if [[ $status -ne 0 ]]; then
    echo "[FATAL] mkcfg failed for '$name' ($*)" >&2
    exit 1
  fi
  echo "$out"
}


# --- Run one cell ------------------------------------------------------------
# run_cell <cell-name> <--anns value> [config path] [extra flags...]
run_cell() {
  local cell="$1" anns="$2" cfg="${3:-$BASE_CFG}"
  shift $(( $# < 3 ? $# : 3 ))
  local extra=("$@")
  local outdir="$RES_ROOT/$cell"
  local abund="$outdir/$SAMPLE/results/reports/abund_table.tsv"

  if [[ -f "$abund" ]]; then
    echo "[skip] $cell (abund_table.tsv exists)"
    return 0
  fi
  echo "[run ] $cell  --anns $anns ${extra[*]:-} config=$(basename "$cfg")"
  mkdir -p "$outdir" "$LOG_DIR/$cell"

  python3 jobs/bench_time.py --csv "$TIMING_CSV" --label "$cell" -- \
    emtrace run \
      -1 "$READS" \
      -s "$SAMPLE" \
      -r "$outdir" \
      -t "$THREADS" \
      --config "$cfg" \
      --anns "$anns" \
      "${extra[@]}" \
      --pipe pore,chop,ann,em,rep,plot \
      --logs "$LOG_DIR/$cell" \
      -v \
    2>&1 | tee "$LOG_DIR/${cell}.log"

  [[ -f "$abund" ]] || echo "[WARN] $cell produced no abund_table.tsv — see $LOG_DIR/${cell}.log"
}


# --- Re-score an existing cell without re-annotating --------------------------
# rerun_cell <cell-name> <source-cell> <config path> [extra flags...]
#
# For settings that act downstream of annotation: --nospnorep (em + rep),
# report.min_qcov (rep) and em.temp (em). The source cell's annotation/ is
# symlinked, not copied; `--pipe em,rep` only reads it, and report.sort_reads is
# false in the 16S config, so nothing writes back into the source tree.
#
# NOT usable for --tophits, -N or s: those change what annotation writes.
rerun_cell() {
  local cell="$1" src="$2" cfg="$3"; shift 3
  local extra=("$@")
  local outdir="$RES_ROOT/$cell"
  local abund="$outdir/$SAMPLE/results/reports/abund_table.tsv"
  local src_annot="$RES_ROOT/$src/$SAMPLE/annotation"

  if [[ -f "$abund" ]]; then
    echo "[skip] $cell (abund_table.tsv exists)"
    return 0
  fi
  if [[ ! -d "$src_annot" ]]; then
    echo "[WARN] $cell: source annotation missing ($src_annot) — run cell '$src' first" >&2
    return 0
  fi

  echo "[rep ] $cell  reuse=$src  config=$(basename "$cfg") ${extra[*]:-}"
  mkdir -p "$outdir/$SAMPLE" "$LOG_DIR/$cell"
  # Absolute target: the symlink is read through a relative -r path.
  ln -sfn "$(cd "$src_annot" && pwd)" "$outdir/$SAMPLE/annotation"

  python3 jobs/bench_time.py --csv "$TIMING_CSV" --label "$cell" -- \
    emtrace run \
      -1 "$READS" \
      -s "$SAMPLE" \
      -r "$outdir" \
      -t "$THREADS" \
      --config "$cfg" \
      "${extra[@]}" \
      --pipe em,rep \
      --logs "$LOG_DIR/$cell" \
      -v \
    2>&1 | tee "$LOG_DIR/${cell}.log"

  [[ -f "$abund" ]] || echo "[WARN] $cell produced no abund_table.tsv — see $LOG_DIR/${cell}.log"
}


# =============================================================================
# Every selected cell group, for every selected library
# =============================================================================
LIBRARIES=()
[[ "$RUN_BARCODE79" == true ]] && LIBRARIES+=("barcode79")
[[ "$RUN_SYNTHGUT"  == true ]] && LIBRARIES+=("synthgut")
[[ -n "${SAMPLE:-}" ]] && LIBRARIES=("$SAMPLE")   # explicit override
[[ ${#LIBRARIES[@]} -gt 0 ]] || { echo "[FATAL] no library selected." >&2; exit 1; }

echo "Libraries: ${LIBRARIES[*]}"
echo

for SAMPLE in "${LIBRARIES[@]}"; do

READS="${READS_OVERRIDE:-$(reads_for "$SAMPLE")}"
if [[ ! -f "$READS" ]]; then
  echo "[WARN] $SAMPLE: reads not found (${READS:-<unknown library>}) — skipping."
  echo "       stage them, e.g.  rsync -avP <local>.fq.gz roihu:$SRA_DIR/"
  continue
fi

# Per-library roots. One tree keyed library -> tool -> configuration, so a path
# states which cell it is without consulting a log (roadmap section 9).
RES_ROOT="$RESULTS_BASE/framework_16S/$SAMPLE/emtrace"
LOG_DIR="$LOG_BASE/$SAMPLE"
mkdir -p "$RES_ROOT" "$LOG_DIR" "$CFG_DIR" "$METRICS_DIR"
case "$(cd "$RES_ROOT" && pwd -P)" in
  /projappl/*) echo "[FATAL] RES_ROOT resolves under /projappl: $RES_ROOT" >&2; exit 1 ;;
esac

echo "=== $SAMPLE  ($(basename "$READS")) ==="


# --- Frozen cells: <reference>.N25.T10 ---------------------------------------
# The configuration Panel A settled on, run unchanged on every selected
# reference. emu is the Panel B.2 cell; rrnDB.dedup is the Panel A frozen cell.
if [[ "$RUN_FROZEN" == true ]]; then
  frozen_cfg=$(mkcfg "frozen.N${FROZEN_N}.T${FROZEN_TEMP}" \
                 "annotation.minimap.max_target_seqs=${FROZEN_N}" \
                 "em.temp=${FROZEN_TEMP}")
  for ref in rrnDB.dedup emu rrnDB refseq chain1 silva; do
    case "$ref" in
      rrnDB.dedup) [[ "$RUN_RRNDB_DEDUP" == true ]] || continue ;;
      emu)         [[ "$RUN_EMU"         == true ]] || continue ;;
      rrnDB)       [[ "$RUN_RRNDB"       == true ]] || continue ;;
      refseq)      [[ "$RUN_REFSEQ16S"   == true ]] || continue ;;
      chain1)      [[ "$RUN_CHAIN"       == true ]] || continue ;;
      silva)       [[ "$RUN_SILVA"       == true ]] || continue ;;
    esac
    run_cell "${ref}.N${FROZEN_N}.T${FROZEN_TEMP}" "${CELL_ANNS[$ref]}" "$frozen_cfg"
  done
fi


# =============================================================================
# Panel A, axis 1 — one cell per reference DB at the shipped defaults
# =============================================================================
if [[ "$RUN_MATRIX" == true ]]; then

# De-duplicated rrnDB first: it is the Section B base, and the cheapest rrnDB cell.
[[ "$RUN_RRNDB_DEDUP" == true ]] && run_cell "rrnDB.dedup" "16S.rrnDB.dedup"
[[ "$RUN_RRNDB"       == true ]] && run_cell "rrnDB"       "16S.rrnDB"
[[ "$RUN_REFSEQ16S"   == true ]] && run_cell "refseq"      "16S.refseq"

# Chain order is not cosmetic: rrnDB first. The legacy RefSeq-first order pushed
# 47.8 % of reads into Other on the same two databases (the project notes).
[[ "$RUN_CHAIN" == true ]] && run_cell "chain1" "16S.chain1"

[[ "$RUN_SILVA" == true ]] && run_cell "silva" "16S.silva"

fi  # RUN_MATRIX


# =============================================================================
# Axis 9 — hit mode: all hits (Section A default) vs --tophits
#
# Top-hit filtering happens in the annotate step, so each variant re-annotates.
# annot_table.tsv and abund_table readn should match the all-hits cell except
# for tied reads; EM prob is what differs (the project notes, hit-mode change).
# =============================================================================
if [[ "$RUN_TOPHITS" == true ]]; then
for src in "${TOPHITS_CELLS[@]}"; do
  run_cell "${src}.tophits" "${CELL_ANNS[$src]}" "$BASE_CFG" "--tophits"
done
fi


# =============================================================================
# Section B — secondary alignments (axis 6) on $SWEEP_BASE_CELL
# =============================================================================
if [[ "$RUN_SWEEPS" == true ]]; then

# --- Axis 6: secondary alignments retained per read --------------------------
# 0, 25 and 50 (5 = base cell). 0 and 50 are the operating points of Spaghetti
# (--secondary=no) and Emu. With all hits kept, extra alignments reach EM as
# ranked alternatives; on the de-duplicated reference they are genuinely
# different sequences rather than identical copies. Re-annotates.
if [[ "$SWEEP_SECONDARY" == true ]]; then
for N in 0 25 50; do
  CFG=$(mkcfg "N${N}" "annotation.minimap.max_target_seqs=${N}")
  run_cell "sweep.secondary.N${N}" "$SWEEP_BASE_ANNS" "$CFG"
done
fi

fi  # RUN_SWEEPS


# =============================================================================
# Tuned cells — best axis-6 and axis-8 settings combined
#
# Run v3 measured -N 25 and em.temp 50 separately (each the best of its axis on
# the EM posterior; the em.temp trend had not flattened at 50, hence 100 too).
# em.temp acts only in EM, so each cell re-runs em,rep over the
# sweep.secondary.N25 annotation. The derived config also sets
# max_target_seqs=25 so it records the full setting, although em,rep never
# reads it.
# =============================================================================
if [[ "$RUN_TUNED" == true ]]; then
for T in "${TUNED_TEMPS[@]}"; do
  CFG=$(mkcfg "N25_temp${T}" "annotation.minimap.max_target_seqs=25" "em.temp=${T}")
  rerun_cell "tuned.N25.temp${T}" "sweep.secondary.N25" "$CFG"
done
fi  # RUN_TUNED


done   # SAMPLE


# =============================================================================
echo
echo "Results:         $RESULTS_BASE/framework_16S/<library>/emtrace/"
echo "Derived configs: $CFG_DIR/"
[[ -f "$TIMING_CSV" ]] && echo "Timings:         $TIMING_CSV"
echo
echo "Next:"
echo "  point jobs/framework_16S/framework_16S.manifest.tsv at the new cells, then"
echo "  python3 jobs/framework_16S/framework_16S.score.py --manifest jobs/framework_16S/framework_16S.manifest.tsv"
echo "  Rscript rplots/mock16S_composition.R"
