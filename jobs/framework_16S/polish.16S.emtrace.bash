#!/bin/bash -l
#SBATCH --job-name=polish_16S
#SBATCH --account=project_XXXXXXX
#SBATCH --partition=interactive
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=24G
#SBATCH --chdir=/projappl/project_XXXXXXX/EMtrace-benchmark
#SBATCH --output=/scratch/project_XXXXXXX/emtrace-benchmark/logs/polish_16S.%j.out
# ============================================================================
# 16S polishing — R10.4.1 benchmark: EMtrace cells   (Roihu)
# ============================================================================
# Plan: the project notes. ONT open data zymo_16s_2025.09
# (R10.4.1, Dorado 1.1.1, SQK-MAB114.24), one library per replicate:
#   gutmock_<model>_r1..r4   ZymoBIOMICS Gut Microbiome Standard (D6331)
#   msplus_<model>_r1..r6    Zymo DNA standard + 4 ATCC species
#
# FROZEN CONFIGURATION ONLY, no retuning (the project notes): the Panel A
# settings, unchanged -- porechop_abi, chopper Q10 and 1,200-1,800 bp, -N 25,
# all hits, em.temp 10. Two cells per library, as in Table 1 and
# Supplementary Table 4:
#   rrnDB.dedup.N25.T10   --anns 16S.rrnDB.dedup   Panel B.1 (own database)
#   emu.N25.T10           --anns 16S.emu           Panel B.2 (Emu's reference)
#
# READS: $SRA_DIR/<library>.fastq.gz, written on the Mac by
#   bash jobs/framework_16S/polish.16S.comparators.bash --part fetch
# and staged with
#   rsync -avP data/sra/16S/gutmock_sup_r*.fastq.gz roihu:$SRA_DIR/
# The FASTQ comments keep Dorado's st and qs tags; chopper filters on the
# quality strings, not on qs.
#
# Switches (environment overrides, no edit needed):
#   COMMUNITY=gutmock|msplus   (default gutmock)
#   MODEL=sup|hac              (default sup)
#   REPS="1 2 3 4"             (default: every replicate of the community)
#   RUN_RRNDB_DEDUP=true|false, RUN_EMU=true|false
#
# Submit:
#   sbatch jobs/framework_16S/polish.16S.emtrace.bash
#   COMMUNITY=msplus sbatch jobs/framework_16S/polish.16S.emtrace.bash
# (sbatch passes the submitting shell's environment by default.)
#
# Skip logic: a cell whose abund_table.tsv exists is skipped, so a resubmit
# only fills gaps. To force one: rm -rf $RESULTS_BASE/framework_16S/<library>/emtrace/<cell>
#
# Layout -- everything on /scratch, as framework_16S.emtrace.bash:
#   $SCRATCH_BASE/results/framework_16S/<library>/emtrace/<cell>/<library>/results/
#   $SCRATCH_BASE/logs/polish_16S/<library>/          bash + pipeline logs
#   $SCRATCH_BASE/logs/polish_16S/configs/            derived config, provenance
#   $SCRATCH_BASE/metrics/framework_16S/timings/polish.16S.timings.csv
# Pull to the Mac with
#   rsync -avP roihu:$SCRATCH_BASE/results/framework_16S/{gutmock,msplus}_* results/framework_16S/
#
# Resources: ~73-110 k reads per replicate. On barcode79 (100 k reads, 32
# cores) rrnDB.dedup took 27 min at 14 GB; the Emu reference is smaller. 8
# cells for gutmock is ~4 h; msplus (12 cells) ~6 h. 24 GB covers the peak.
#
# Real-time replay of one replicate is not run here: it follows the Framework C
# conventions (framework_C.replay.bash, laptop), after batching the FASTQ by
# its st tag.
# ============================================================================

# =============================================================================
# Environment on Roihu (as framework_16S.emtrace.bash)
# =============================================================================
# Do not run `module purge` on Roihu: it removes the module tree these come from.
module load bio-apps minimap2 samtools seqkit
export PATH="/projappl/project_XXXXXXX/bin:/projappl/project_XXXXXXX/install/chopper_dir:$PATH"
source /projappl/project_XXXXXXX/venv/emtrace/bin/activate
# The modules set PYTHONPATH to a different numpy, which breaks EMtrace.
unset PYTHONPATH
export EMTRACE_DB_DIR=/scratch/project_XXXXXXX/databases/emtrace

set -uo pipefail

for tool in emtrace minimap2 porechop_abi chopper; do
  command -v "$tool" > /dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done
for idx in minimap.rrnDB_16S_rRNA.dedup.fasta.mmi minimap.emu.fasta.mmi; do
  [[ -f "$EMTRACE_DB_DIR/$idx" ]] || echo "[WARN] index not found: $EMTRACE_DB_DIR/$idx"
done

echo "Job ${SLURM_JOB_ID:-local} on $(hostname), started $(date)"

COMMUNITY="${COMMUNITY:-gutmock}"
MODEL="${MODEL:-sup}"
case "$COMMUNITY" in
  gutmock) NREPS=4 ;;
  msplus)  NREPS=6 ;;
  *) echo "[FATAL] COMMUNITY must be gutmock|msplus" >&2; exit 1 ;;
esac
case "$MODEL" in sup|hac) ;; *) echo "[FATAL] MODEL must be sup|hac" >&2; exit 1 ;; esac
REPS="${REPS:-$(seq 1 "$NREPS" | tr '\n' ' ')}"
RUN_RRNDB_DEDUP="${RUN_RRNDB_DEDUP:-true}"
RUN_EMU="${RUN_EMU:-true}"

# The frozen configuration (the project notes).
FROZEN_N=25
FROZEN_TEMP=10
BASE_CFG="jobs/framework_16S/emtrace.config.16S.yaml"

SCRATCH_BASE="/scratch/project_XXXXXXX/emtrace-benchmark"
SRA_DIR="$SCRATCH_BASE/sra"
RESULTS_BASE="$SCRATCH_BASE/results"
LOG_BASE="$SCRATCH_BASE/logs/polish_16S"
CFG_DIR="$LOG_BASE/configs"
METRICS_DIR="$SCRATCH_BASE/metrics/framework_16S/timings"
TIMING_CSV="$METRICS_DIR/polish.16S.timings.csv"
THREADS=${SLURM_CPUS_PER_TASK:-8}
mkdir -p "$CFG_DIR" "$METRICS_DIR"

# --- Derive the frozen config ------------------------------------------------
# Same derivation as framework_16S.emtrace.bash (mkcfg): YAML parse, not sed,
# and a missing key stops the run.
FROZEN_CFG="$CFG_DIR/emtrace.config.16S.frozen.N${FROZEN_N}.T${FROZEN_TEMP}.yaml"
python3 - "$BASE_CFG" "$FROZEN_CFG" "$FROZEN_N" "$FROZEN_TEMP" <<'CFGPY' || { echo "[FATAL] config derivation failed" >&2; exit 1; }
import sys, yaml
src, dst, n, temp = sys.argv[1:]
cfg = yaml.safe_load(open(src))
cfg["annotation"]["minimap"]["max_target_seqs"] = int(n)
cfg["em"]["temp"] = float(temp)
ch = cfg["preprocess"]["chopper"]
assert (ch["quality"], ch["minlength"], ch["maxlength"]) == (10, 1200, 1800), \
    f"base config no longer frozen: chopper {ch}"
yaml.safe_dump(cfg, open(dst, "w"), sort_keys=False, default_flow_style=False)
CFGPY
echo "config: $FROZEN_CFG"

# --- Run one cell ------------------------------------------------------------
# run_cell <cell-name> <--anns value>
run_cell() {
  local cell="$1" anns="$2"
  local outdir="$RES_ROOT/$cell"
  local abund="$outdir/$SAMPLE/results/reports/abund_table.tsv"
  if [[ -f "$abund" ]]; then
    echo "[skip] $SAMPLE/$cell (abund_table.tsv exists)"; return 0
  fi
  echo "[run ] $SAMPLE/$cell  --anns $anns"
  mkdir -p "$outdir" "$LOG_DIR/$cell"
  python3 jobs/bench_time.py --csv "$TIMING_CSV" --label "$SAMPLE.$cell" -- \
    emtrace run \
      -1 "$READS" \
      -s "$SAMPLE" \
      -r "$outdir" \
      -t "$THREADS" \
      --config "$FROZEN_CFG" \
      --anns "$anns" \
      --pipe pore,chop,ann,em,rep,plot \
      --logs "$LOG_DIR/$cell" \
      -v \
    2>&1 | tee "$LOG_DIR/${cell}.log"
  [[ -f "$abund" ]] || echo "[WARN] $SAMPLE/$cell produced no abund_table.tsv — see $LOG_DIR/${cell}.log"
}

# =============================================================================
# Every replicate
# =============================================================================
MANIFEST_ROWS=()
for r in $REPS; do
  SAMPLE="${COMMUNITY}_${MODEL}_r${r}"
  READS="$SRA_DIR/$SAMPLE.fastq.gz"
  if [[ ! -s "$READS" ]]; then
    echo "[WARN] $SAMPLE: $READS not found — skipping."
    echo "       stage it: rsync -avP data/sra/16S/$SAMPLE.fastq.gz roihu:$SRA_DIR/"
    continue
  fi
  RES_ROOT="$RESULTS_BASE/framework_16S/$SAMPLE/emtrace"
  LOG_DIR="$LOG_BASE/$SAMPLE"
  mkdir -p "$RES_ROOT" "$LOG_DIR"
  case "$(cd "$RES_ROOT" && pwd -P)" in
    /projappl/*) echo "[FATAL] RES_ROOT resolves under /projappl: $RES_ROOT" >&2; exit 1 ;;
  esac
  echo; echo "=== $SAMPLE  ($(basename "$READS"), $(( $(gzip -dc "$READS" | wc -l) / 4 )) reads) ==="

  CELL_B1="rrnDB.dedup.N${FROZEN_N}.T${FROZEN_TEMP}"
  CELL_B2="emu.N${FROZEN_N}.T${FROZEN_TEMP}"
  if [[ "$RUN_RRNDB_DEDUP" == true ]]; then
    run_cell "$CELL_B1" "16S.rrnDB.dedup"
    MANIFEST_ROWS+=("$(printf 'EMtrace.rrnDB.dedup\t%s\tdiag\temtrace\tncbi\tresults/framework_16S/%s/emtrace/%s/%s/results/reports/abund_table.tsv\tR10.4.1 (the project notes): frozen, rrnDB.dedup, -N 25, em.temp 10' "$SAMPLE" "$SAMPLE" "$CELL_B1" "$SAMPLE")")
  fi
  if [[ "$RUN_EMU" == true ]]; then
    run_cell "$CELL_B2" "16S.emu"
    MANIFEST_ROWS+=("$(printf 'EMtrace.emu\t%s\tdiag\temtrace\tncbi\tresults/framework_16S/%s/emtrace/%s/%s/results/reports/abund_table.tsv\tR10.4.1 (the project notes): frozen, Emu reference, -N 25, em.temp 10' "$SAMPLE" "$SAMPLE" "$CELL_B2" "$SAMPLE")")
  fi
  # Record what the frozen chopper settings kept.
  grep -h -m1 "chopper: kept" "$LOG_DIR"/*.log 2>/dev/null | sed "s/^/  $SAMPLE /" || true
done

# =============================================================================
echo
echo "Results: $RESULTS_BASE/framework_16S/${COMMUNITY}_${MODEL}_r*/emtrace/"
[[ -f "$TIMING_CSV" ]] && echo "Timings: $TIMING_CSV"
echo
echo "Manifest rows (panel diag; paths as on the Mac after rsync):"
for row in "${MANIFEST_ROWS[@]+"${MANIFEST_ROWS[@]}"}"; do echo "$row"; done
echo
echo "Next, on the Mac:"
echo "  rsync -avP roihu:$RESULTS_BASE/framework_16S/${COMMUNITY}_${MODEL}_r* results/framework_16S/"
echo "  python3 jobs/framework_16S/framework_16S.score.py --manifest jobs/framework_16S/framework_16S.manifest.tsv"
