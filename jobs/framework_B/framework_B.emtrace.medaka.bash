#!/usr/bin/env bash
# ============================================================================
# Framework B — medaka-polished EMtrace assemblies (v20260924 batches)
# ============================================================================
# Re-runs ONLY the assembly step (`--pipe ass --medaka`) on finished
# v20260924 samples, so the polished consensus can be scored against the
# unpolished one on identical inputs:
#   reference  the sample's own assembly/combined_ref.fa (--refseqs), i.e. the
#              reference ref_auto chose in the original run; ref_auto is off
#   reads      reads/reads.trim.hflt.fq.gz (assembly.input_reads = bgfiltered,
#              as in the original run)
#   settings   min_coverage MIN_COVERAGE (10), everything else from CONFIG
#   polishing  (steps/assembly.py, _polish_medaka)
#                draft   = dominant-allele VCF on the UNMASKED reference
#                reads   -> draft (minimap2 -k 11 -w 5, assembly.minimap)
#                medaka inference on draft regions with depth >= min_coverage
#                medaka sequence --fill_char N (rest of the draft -> N)
#                reads   -> polished; depth < min_coverage re-masked with N
#
# Output: {batch}/{sample}/assembly.medaka/ — the existing assembly/ is never
# touched.  The run happens in a staging copy (STAGE_ROOT/{batch}/{sample}:
# job.yaml + a symlink to the sample's reads/), whose assembly/ is then moved
# to assembly.medaka/.  Key files there:
#   consensus.fa             polished, masked consensus (score this)
#   consensus.bcftools.fa    unpolished; should equal assembly/consensus.fa
#                            (checked below and recorded in the stamp as
#                            repro=identical|DIFFERS)
#   consensus.medaka.unmasked.fa, reads.polished.md.bam, *.lowcovmask.bed
# Samples without assembly/combined_ref.fa (ref_auto selected no reference)
# are skipped: there is nothing to polish.
#
# Medaka model:
#   simulated mutants  auto — the reads carry a basecaller tag
#                      (r1041_e82_400bps_hac_g632)
#   SRA libraries      headers stripped by SRA, so set per accession in
#                      model_for() below.  No chemistry metadata is available;
#                      2020 MinION runs (Gauthier, SRR15356294) are R9.4.1 ->
#                      r941_min_hac_g507; the 2024-25 Petersen libraries (median
#                      read Q 24-36) are taken as R10.4.1 sup -> medaka default.
#
# Medaka lives in its own conda env (MEDAKA_BIN), appended to PATH so it
# never shadows the Homebrew samtools/bcftools/python.
#
# Skip logic: a sample is skipped once <sample>/assembly.medaka/.done.medaka
# exists.  To redo: rm -rf <batch>/<sample>/assembly.medaka
#
# Timing:  metrics/framework_B/<batch>.medaka.timings.csv (bench_time.py)
# Logs:    logs/<batch>.medaka/<sample>.log            (tee)
#          logs/emtrace/<batch>.medaka/<sample>/    (--logs)
#
# Run from the repository root:
#   bash jobs/framework_B/framework_B.emtrace.medaka.bash
# ============================================================================

set -o pipefail

[[ -d jobs ]] || { echo "error: run this script from the repository root"; exit 1; }

VERSION="v20260924"
RESULTS_ROOT="results/framework_b"
STAGE_ROOT="$RESULTS_ROOT/.medaka_stage"
LOG_DIR="logs"
PIPE_LOG_DIR="logs/emtrace"
METRICS_DIR="metrics/framework_B"
CONFIG="jobs/emtrace.config.yaml"
MEDAKA_BIN="$HOME/miniconda3/envs/medaka/bin"
THREADS=8
MIN_COVERAGE=10          # must match the original runs (assembly/ holds d10)

# =============================================================================
# Control switches — set to true/false before running
# =============================================================================
# simulated mutants, per annotation cell (see framework_B.emtrace.v20260924.bash)
RUN_PURE_REFAUTO=false     # pure_v20260924_refauto   (28 samples)
RUN_PURE_CHAIN2=false      # pure_v20260924_chain2    (28)
RUN_PURE_K11W5=false      # pure_v20260924_k11w5     (14, Lassa only)
RUN_MIXED_REFAUTO=false    # mixed_v20260924_refauto  (28)
RUN_MIXED_CHAIN2=false     # mixed_v20260924_chain2   (28)
RUN_MIXED_K11W5=false      # mixed_v20260924_k11w5    (14, Lassa only)
# real libraries
RUN_PETERSEN_REFAUTO=false # sra_petersen_2026_v20260924_refauto (4)
RUN_PETERSEN_CHAIN2=false  # sra_petersen_2026_v20260924_chain2  (4)
RUN_GAUTHIER_REFAUTO=true # sra_gauthier_2021_v20260924_refauto (21, R9.4.1)

DRY_RUN="${DRY_RUN:-false}"   # list what would run; touch nothing (env DRY_RUN=true)
# =============================================================================

# model_for <sample>: medaka model for --medaka_model
model_for() {
  case "$1" in
    SRR153562*|SRR153563*) echo "r941_min_hac_g507" ;;   # Gauthier 2021 / PRJNA752146, R9.4.1
    SRR*|ERR*)             echo "r1041_e82_400bps_sup_v5.2.0" ;;
    *)                     echo "auto" ;;
  esac
}

[[ -x "$MEDAKA_BIN/medaka" ]] || { echo "error: medaka not found in $MEDAKA_BIN"; exit 1; }
export PATH="$PATH:$MEDAKA_BIN"

RAN=0; SKIPPED=0; FAILED=0; NOREF=0; DIFFERS=0
FAILED_SAMPLES=(); DIFFERS_SAMPLES=()

# ---------------------------------------------------------------------------
# run_sample <batch> <sample>
# ---------------------------------------------------------------------------
run_sample() {
  local batch="$1" sample="$2"
  local rdir="$RESULTS_ROOT/$batch/$sample"
  local out="$rdir/assembly.medaka" stamp="$rdir/assembly.medaka/.done.medaka"
  local ref="$rdir/assembly/combined_ref.fa"

  if [[ -f "$stamp" ]]; then
    echo "[skip] $batch/$sample — done"; SKIPPED=$((SKIPPED + 1)); return 0
  fi
  if [[ ! -s "$ref" ]]; then
    echo "[skip] $batch/$sample — no reference selected in the original run"
    NOREF=$((NOREF + 1)); return 0
  fi
  local fastq
  fastq=$(awk -F': ' '$1=="input_reads"{print $2}' "$rdir/job.yaml")
  if [[ ! -s "$rdir/reads/reads.trim.hflt.fq.gz" ]]; then
    echo "[skip] $batch/$sample — reads/reads.trim.hflt.fq.gz missing"; SKIPPED=$((SKIPPED + 1)); return 0
  fi
  local model; model=$(model_for "$sample")
  if [[ "$DRY_RUN" == true ]]; then
    echo "[dry-run] $batch/$sample  --refseqs $ref --medaka_model $model"
    return 0
  fi

  # staging copy: job.yaml + reads symlink; nothing under $rdir is written
  local stage="$STAGE_ROOT/$batch/$sample"
  rm -rf "$stage"; mkdir -p "$stage"
  cp "$rdir/job.yaml" "$stage/"
  ln -s "$PWD/$rdir/reads" "$stage/reads"

  mkdir -p "$LOG_DIR/$batch.medaka" "$PIPE_LOG_DIR/$batch.medaka"
  echo "[run]  $batch/$sample  (model $model)"
  local rc=0
  python3 jobs/bench_time.py --csv "$METRICS_DIR/${batch}.medaka.timings.csv" --label "$sample" -- \
    emtrace run \
    -1 "$fastq" \
    -s "$sample" \
    -r "$STAGE_ROOT/$batch" \
    -t "$THREADS" \
    --refseqs "$PWD/$ref" \
    --no-ref_auto \
    --min_coverage "$MIN_COVERAGE" \
    --medaka --medaka_model "$model" \
    --config "$CONFIG" \
    --pipe ass \
    --logs "$PIPE_LOG_DIR/$batch.medaka/$sample" \
    -v \
    2>&1 | tee "$LOG_DIR/$batch.medaka/$sample.log" > /dev/null || rc=$?

  if [[ $rc -ne 0 || ! -f "$stage/assembly/consensus.bcftools.fa" ]]; then
    [[ $rc -eq 0 ]] && echo "[FAIL] $batch/$sample — no polished consensus (medaka failed or skipped; see log)"
    [[ $rc -ne 0 ]] && echo "[FAIL] $batch/$sample — exit $rc (see $LOG_DIR/$batch.medaka/$sample.log)"
    FAILED=$((FAILED + 1)); FAILED_SAMPLES+=("$batch/$sample"); return 0
  fi

  rm -rf "$out"
  mv "$stage/assembly" "$out"
  rm -rf "$stage"

  local repro="identical"
  if ! cmp -s "$out/consensus.bcftools.fa" "$rdir/assembly/consensus.fa"; then
    repro="DIFFERS"; DIFFERS=$((DIFFERS + 1)); DIFFERS_SAMPLES+=("$batch/$sample")
    echo "[note] $batch/$sample — unpolished consensus differs from assembly/consensus.fa"
  fi
  local used
  used=$(grep -o "polishing (model [^)]*)" "$LOG_DIR/$batch.medaka/$sample.log" | head -1)
  printf 'date=%s\nmodel=%s\nrepro=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "${used:-$model}" "$repro" > "$stamp"
  RAN=$((RAN + 1))
}

# run_batch <batch>: every sample directory of the batch
run_batch() {
  local batch="$1" d
  if [[ ! -d "$RESULTS_ROOT/$batch" ]]; then
    echo "=== $batch — not found, skipped ==="; return 0
  fi
  echo "=== $batch ==="
  for d in "$RESULTS_ROOT/$batch"/*/; do
    [[ -f "$d/job.yaml" ]] || continue
    run_sample "$batch" "$(basename "$d")"
  done
}

mkdir -p "$METRICS_DIR"
echo "Framework B medaka polishing of $VERSION assemblies, min_coverage $MIN_COVERAGE"
echo "EMtrace: $(emtrace --version 2>/dev/null || echo '?');  $("$MEDAKA_BIN/medaka" --version 2>/dev/null)"

[[ "$RUN_PURE_REFAUTO"     == true ]] && run_batch "pure_${VERSION}_refauto"
[[ "$RUN_PURE_CHAIN2"      == true ]] && run_batch "pure_${VERSION}_chain2"
[[ "$RUN_PURE_K11W5"       == true ]] && run_batch "pure_${VERSION}_k11w5"
[[ "$RUN_MIXED_REFAUTO"    == true ]] && run_batch "mixed_${VERSION}_refauto"
[[ "$RUN_MIXED_CHAIN2"     == true ]] && run_batch "mixed_${VERSION}_chain2"
[[ "$RUN_MIXED_K11W5"      == true ]] && run_batch "mixed_${VERSION}_k11w5"
[[ "$RUN_PETERSEN_REFAUTO" == true ]] && run_batch "sra_petersen_2026_${VERSION}_refauto"
[[ "$RUN_PETERSEN_CHAIN2"  == true ]] && run_batch "sra_petersen_2026_${VERSION}_chain2"
[[ "$RUN_GAUTHIER_REFAUTO" == true ]] && run_batch "sra_gauthier_2021_${VERSION}_refauto"
rmdir "$STAGE_ROOT"/* "$STAGE_ROOT" 2>/dev/null

echo
echo "=== summary ==="
echo "ran: $RAN   skipped: $SKIPPED   no reference: $NOREF   failed: $FAILED   unpolished differs: $DIFFERS"
[[ $FAILED  -gt 0 ]] && printf '  failed: %s\n' "${FAILED_SAMPLES[@]}"
[[ $DIFFERS -gt 0 ]] && printf '  differs: %s\n' "${DIFFERS_SAMPLES[@]}"
echo "timings: $METRICS_DIR/<batch>.medaka.timings.csv"
