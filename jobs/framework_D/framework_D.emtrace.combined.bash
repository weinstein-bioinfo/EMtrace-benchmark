#!/usr/bin/env bash
# ============================================================================
# Framework D — host-filter -s optimization on combined Nanopore datasets (final)
# ============================================================================
#
# Sweeps the minimap2 DP alignment-score threshold (-s, with -c on so true AS
# scores gate filtering) to find the -s that best separates host from virus on
# realistic Nanopore data. Runs ONLY the background-filter (flt) step of
# emtrace on the labeled multi-virus datasets built by
# D_sim06_combine_labeled.py: each `combined_<rate>.fastq.gz` holds COVID + DENV
# + LASV_L + LASV_S (all at one mutation rate) + human+bacteria background, with
# a per-read division truth table (`combined_<rate>.truth.tsv`).
#
# Sweep:  -s ∈ {0,30,50,75,100,125}, -c on, -m held at 0 (so -s is the sole gate),
#         in two seeding configs:
#           default   --preset map-ont                (k/w = minimap default)
#           k11w5     --preset map-ont -k 11 -w 5     (sensitive seeding)
#         minimap2 -c/-s/-k/-w are NOT on the emtrace CLI; they live only in
#         the config (preprocess.bgfilter.minimap.*), so each setting is run by
#         patching the base config into a generated YAML (gen_config) → --config.
#         The k/w-aware .mmi cache (<db>.mmi vs <db>.k11w5.mmi) is built once per
#         config; the -s sweep reuses it (only -s changes at map time).
#
# --pipe flt ONLY: no Porechop/Chopper before, no annotation/EM/assembly after.
# With preprocess.subsample.active=false the flt step reads the raw input FASTQ
# directly (resolve_preprocess_input → input_reads), so all ~502k reads are
# filtered untrimmed.
#
# Scoring: per rate, all (config × -s) removed-read FASTQs are scored against
# that rate's division truth with `ngs-benchmark bgfilter-divbatch`, giving one
# table of fraction-removed per division (Viruses / Bacteria / Eukaryota) — pick
# the -s that maximises Bacteria+Eukaryota removal at minimal Viruses loss.
#
# Skip logic: a setting is skipped (run) if its removed FASTQ already exists, but
# it is still added to the scoring manifest. Force a re-run by deleting
# results/framework_d/combined/<sample>.
#
# Outputs:
#   metrics/framework_D/framework_d_combined.<rate>.table.tsv    optimization table
#   metrics/framework_D/framework_d_combined.<rate>.manifest.tsv
#   results/framework_d/combined/<sample>/reads/...          filtered reads
#   logs under logs/{,emtrace/}framework_d/

set -euo pipefail

# =============================================================================
# Control switches
# =============================================================================
RATES=(00 10)                      # mutation rate(s) to optimize on (subset of 00 10 20 30)
RUN_DEFAULT=true                # config: map-ont (k/w = minimap default)
RUN_K11W5=false                 # config: map-ont -k 11 -w 5
THRESH=(0 30 50 75 100 125)     # minimap2 -s DP-score thresholds to sweep
M_FIXED=0                       # minimap2 -m held fixed so -s is the sole gate
THREADS=8
# =============================================================================

DATA_DIR="data/petersen_2026/combined"
RESULTS_DIR="results/framework_d/combined"
BASE_CONFIG="jobs/emtrace.config.yaml"
BGFILTER="Homo_sapiens"

METRICS_DIR="metrics/framework_D"
LOG_DIR="logs/framework_d"
PIPE_LOG_DIR="logs/emtrace/framework_d"
WORK_DIR="$RESULTS_DIR/_configs"     # generated per-setting config YAMLs

mkdir -p "$RESULTS_DIR" "$METRICS_DIR" "$LOG_DIR" "$PIPE_LOG_DIR" "$WORK_DIR"

# -----------------------------------------------------------------------------
# gen_config CFG_LABEL S K W  ->  prints generated config path
#   Patches preprocess.bgfilter.minimap (c/s/m/k/w) onto the full base config
#   (host.databases etc. preserved — emtrace does NOT merge partials).
# -----------------------------------------------------------------------------
gen_config() {
  local cfg="$1" s="$2" k="$3" w="$4"
  local out="$WORK_DIR/config.${cfg}.AS${s}.yaml"
  BASE="$BASE_CONFIG" OUT="$out" S="$s" K="$k" W="$w" M="$M_FIXED" python3 - <<'PY'
import os, yaml
with open(os.environ["BASE"]) as fh:
    cfg = yaml.safe_load(fh)
mm = cfg.setdefault("preprocess", {}).setdefault("bgfilter", {})
mm["aligner"] = "minimap2"
mm["bgfilter"] = "Homo_sapiens"
mp = mm.setdefault("minimap", {})
mp["c"] = True                       # full DP -> AS:i: scores produced
mp["s"] = int(os.environ["S"])       # swept DP alignment-score threshold
mp["m"] = int(os.environ["M"])       # held fixed so -s is the sole gate
mp["k"] = int(os.environ["K"])       # 0 = minimap default (15)
mp["w"] = int(os.environ["W"])       # 0 = minimap default (10)
mp["save_mmi_index"] = True          # build/reuse a (k/w-aware) .mmi index
with open(os.environ["OUT"], "w") as fh:
    yaml.safe_dump(cfg, fh, sort_keys=False)
PY
  echo "$out"
}

# -----------------------------------------------------------------------------
# run_setting RATE CFG_LABEL K W S MANIFEST
#   Runs one flt setting (skips if already done) and appends it to MANIFEST.
# -----------------------------------------------------------------------------
run_setting() {
  local rate="$1" cfg="$2" k="$3" w="$4" s="$5" manifest="$6"
  local sample="combined_${rate}_${cfg}_AS${s}"
  local reads="$DATA_DIR/combined_${rate}.fastq.gz"
  local removed="$RESULTS_DIR/$sample/reads/hflt.${BGFILTER}.reads.fq.gz"

  if [[ -f "$removed" ]]; then
    echo "  [skip] $sample (removed reads already present)"
  else
    local conf; conf=$(gen_config "$cfg" "$s" "$k" "$w")
    echo "  [run]  $sample  (c=on -s $s, k=$k w=$w)"
    emtrace run \
      -1 "$reads" \
      -s "$sample" \
      -r "$RESULTS_DIR" \
      -t "$THREADS" \
      --bgfilter "$BGFILTER" \
      --config "$conf" \
      --pipe flt \
      --logs "$PIPE_LOG_DIR/${sample}.log" \
      -v \
      &> "$LOG_DIR/${sample}.log" || { echo "  [FAIL] $sample (see log)"; return 0; }
  fi
  printf '%s\t%s\n' "${cfg}_AS${s}" "$removed" >> "$manifest"
}

# =============================================================================
# Run + score, per rate
# =============================================================================
for rate in "${RATES[@]}"; do
  reads="$DATA_DIR/combined_${rate}.fastq.gz"
  truth="$DATA_DIR/combined_${rate}.truth.tsv"
  if [[ ! -f "$reads" || ! -f "$truth" ]]; then
    echo "[skip] rate $rate: dataset/truth missing — run D_sim06_combine_labeled.py" >&2
    continue
  fi

  echo "=== rate $rate: sweeping -s over ${THRESH[*]} ==="
  manifest="$METRICS_DIR/framework_d_combined.${rate}.manifest.tsv"
  printf 'filter\tremoved_fastq\n' > "$manifest"

  for s in "${THRESH[@]}"; do
    [[ "$RUN_DEFAULT" == true ]] && run_setting "$rate" "default" 0  0 "$s" "$manifest"
    [[ "$RUN_K11W5"   == true ]] && run_setting "$rate" "k11w5"   11 5 "$s" "$manifest"
  done

  echo "  [score] rate $rate -> ngs-benchmark bgfilter-divbatch"
  ngs-benchmark bgfilter-divbatch \
    --manifest "$manifest" \
    --truth-table "$truth" \
    --outdir "$WORK_DIR/score.$rate"
  mv "$WORK_DIR/score.$rate/table.tsv" "$METRICS_DIR/framework_d_combined.${rate}.table.tsv"
  echo "  [done] rate $rate -> $METRICS_DIR/framework_d_combined.${rate}.table.tsv"
done

echo "optimization tables in: $METRICS_DIR/framework_d_combined.*.table.tsv"
