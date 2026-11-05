#!/usr/bin/env bash
# ============================================================================
# Framework C — genome coverage over a sequencing run (viral replays)
# ============================================================================
# Plan: the project notes. For each replayed viral library,
# the host-filtered reads of batches 1..k (the replay's own per-batch output,
# pipeline/<lib>/batches/batch_NNNN/reads/reads.trim.hflt.fq.gz) are pooled
# and assembled at log-spaced batches k, so genome coverage can be plotted
# against sequencing time. Watch mode itself is not changed.
#
#   reference   per library (segments 1–10 for the orbiviruses, one for TBEV):
#               Orbivirus_kevo  Kevo orbivirus fi663 PQ165789.1–PQ165798.1 (published case;
#                               not in the replay's annotation db, so assembly is the only
#                               place the true genome enters)
#               TBEV_FBC92293   TBEV Neudoerfl NC_001672.1 (European subtype, 11,141 nt;
#                               published case). The reads sit at ~95.5 % identity to it
#                               (R10.4.1 hac), so the strains differ from Neudoerfl by a few
#                               per cent: identity measures that divergence, not assembly
#                               error. The library pools four isolates (IL209-11-15,
#                               IL232-5, IL239-1, IL256-1), so the consensus is theirs jointly
#               ticks_35-33     CCHFV (Aigai virus) isolate CAP12, L/M/S complete cds
#                               MN811031.1, MN811037.1, MN811034.1: nearest complete relative
#                               of the pools' sequences in Ergunay 2023 Fig. 2
#               ticks_58-37     the same CCHFV plus JMTV strain Tulcea1, segments 1–4 complete
#                               cds MW561147.1–MW561150.1 (Ergunay 2023 Fig. 1), one assembly
#                               (the project notes). The tick runs use the RefSeq-only
#                               replay's reads: host-filtered reads are the same with DIAMOND
#   assembly    emtrace run --pipe ass, --refseqs <reference> (no ref_auto),
#               min_coverage 10 for every library (3 was tried for the tick pools:
#               coverage_over_time.mincov3/, not used),
#               Framework B config (minimap2 -k 11 -w 5), no polishing at the time points
#   end point   the last batch is assembled once more WITH medaka polishing
#               (b<NNNN>.medaka/), to show the accuracy of the final genome; not for
#               the tick pools, whose basecaller model is not known (SRA deposit)
#   scoring     framework_C.coverage_over_time.score.py, per segment, with the
#               Framework B scorer (MAFFT G-INS-i, one alignment at a time)
#
# Output:
#   results/framework_C/<lib>/coverage_over_time/b<NNNN>[.medaka]/
#       consensus.fa, reads.refgen.lowcovmask.bed, .done
#   metrics/framework_C/coverage_over_time.tsv   (one row per lib x batch x segment;
#       rows of libraries switched off here are kept, the rest are re-scored)
# Figure:  Rscript rplots/frameworkC_coverage_over_time.R
#
# Skip logic: a point is skipped once its .done exists. To redo:
#   rm -rf results/framework_C/<lib>/coverage_over_time/b<NNNN>
#
# Run from the repository root:
#   bash jobs/framework_C/framework_C.coverage_over_time.bash
# ============================================================================

set -o pipefail
[[ -d jobs ]] || { echo "error: run this script from the repository root"; exit 1; }

FC="results/framework_C"
ORBI_REF="data/orbivirus.kevo/refs/kevoov_PQ165789-PQ165798.fasta"
# The KEVOOV reads carry the Mk1C on-device model id 2021-05-17_dna_r9.4.1_minion_96_29d8704b,
# which medaka's `auto` cannot map; use the R9.4.1 model the manuscript uses for R9.4.1 data.
ORBI_MEDAKA_MODEL="r941_min_hac_g507"
TBEV_REF="data/tbev.inhouse/refs/tbev_NC_001672.fasta"
# The TBEV replay's synthetic batches rewrite headers to `@id start_time=...`, dropping the
# Dorado RG tag medaka's `auto` reads; the reads are dna_r10.4.1_e8.2_400bps_hac@v5.2.0.
TBEV_MEDAKA_MODEL="r1041_e82_400bps_hac_v5.2.0"
TICK_CCHFV_REF="data/ticks.ergunay/refs/cchfv_aigai_CAP12_MN811031-MN811037-MN811034.fasta"
TICK_58_REF="data/ticks.ergunay/refs/ticks_58-37_cchfv_CAP12+jmtv_tulcea1.fasta"
TICK_MIN_COVERAGE=10
CONFIG="jobs/emtrace.config.yaml"
METRICS="metrics/framework_C/coverage_over_time.tsv"
LOG_DIR="logs/framework_C.coverage_over_time"
MEDAKA_BIN="$HOME/miniconda3/envs/medaka/bin"
THREADS=8
MIN_COVERAGE=10

# =============================================================================
# Control switches
# =============================================================================
RUN_ORBI=false            # Orbivirus_kevo: 116 batches, 19.0 h (published viral case); uses the
                          # RefSeq replay's reads (host-filtered reads are the same in both db cases)
RUN_TBEV=false            # TBEV_FBC92293: 79 batches, 4.98 h (published viral case, 4 isolates pooled)
RUN_TICKS=true            # ticks_35-33 (55 batches), ticks_58-37 (12): field pools, no medaka
N_POINTS=20               # log-spaced batches per library (all batches if fewer)
POLISH_END=true           # also assemble the last batch with medaka
SCORE=true                # score all points into $METRICS
# =============================================================================

[[ -x "$MEDAKA_BIN/medaka" ]] && export PATH="$PATH:$MEDAKA_BIN"
# The scorer imports ngs_benchmark: run it with the ngs-benchmark launcher's own
# interpreter, whatever python is first on PATH.
NB_PY=$(head -1 "$(command -v ngs-benchmark)" | sed -E 's/^#! *//; s/^\/usr\/bin\/env +//')

# points <n_batches>: log-spaced batch numbers 1..n, always including 1 and n
points() {
  python3 -c "
import math, sys
n, k = int(sys.argv[1]), int(sys.argv[2])
if n <= k: print(*range(1, n + 1)); sys.exit()
print(*sorted({round(math.exp(i * math.log(n) / (k - 1))) for i in range(k)}))" "$1" "$N_POINTS"
}

# assemble_point <lib> <ref> <batch> <polish: true|false> [medaka model, default auto]
#                [min coverage, default $MIN_COVERAGE]
assemble_point() {
  local lib="$1" ref="$2" k="$3" polish="$4" model="${5:-auto}" mincov="${6:-$MIN_COVERAGE}"
  local tag; tag=$(printf 'b%04d' "$k"); [[ "$polish" == true ]] && tag="$tag.medaka"
  local out="$FC/$lib/coverage_over_time/$tag"
  [[ -f "$out/.done" ]] && { echo "[skip] $lib $tag"; return 0; }

  local stage="$FC/$lib/coverage_over_time/.stage"
  local sample="$lib.$tag"
  rm -rf "$stage/$sample"; mkdir -p "$stage/$sample/reads" "$LOG_DIR"
  local i
  for ((i = 1; i <= k; i++)); do
    cat "$FC/$lib/pipeline/$lib/batches/$(printf 'batch_%04d' "$i")/reads/reads.trim.hflt.fq.gz"
  done > "$stage/$sample/reads/reads.trim.hflt.fq.gz"

  local polish_flag=(--no-medaka)
  [[ "$polish" == true ]] && polish_flag=(--medaka --medaka_model "$model")
  echo "[run]  $lib $tag"
  local rc=0
  emtrace run \
    -1 "$stage/$sample/reads/reads.trim.hflt.fq.gz" \
    -s "$sample" \
    -r "$stage" \
    -t "$THREADS" \
    --refseqs "$PWD/$ref" \
    --no-ref_auto \
    --min_coverage "$mincov" \
    "${polish_flag[@]}" \
    --config "$CONFIG" \
    --pipe ass \
    --logs "$LOG_DIR/$sample" \
    -v > "$LOG_DIR/$sample.log" 2>&1 || rc=$?

  local asm="$stage/$sample/assembly"
  if [[ $rc -ne 0 || ! -f "$asm/consensus.fa" ]]; then
    echo "[FAIL] $lib $tag — exit $rc (see $LOG_DIR/$sample.log)"; return 0
  fi
  if [[ "$polish" == true && ! -f "$asm/consensus.bcftools.fa" ]]; then
    echo "[FAIL] $lib $tag — not polished (see $LOG_DIR/$sample.log)"; return 0
  fi
  rm -rf "$out"; mkdir -p "$out"
  cp "$asm/consensus.fa" "$asm/reads.refgen.lowcovmask.bed" "$out/"
  date -u +%Y-%m-%dT%H:%M:%SZ > "$out/.done"
  rm -rf "$stage/$sample"
}

# run_lib <lib> <ref> [medaka model, default auto; none = no polished end point]
#         [min coverage, default $MIN_COVERAGE]
run_lib() {
  local lib="$1" ref="$2" model="${3:-auto}" mincov="${4:-$MIN_COVERAGE}" n
  [[ -s "$ref" ]] || { echo "error: reference not found: $ref"; exit 1; }
  # A replay still running (or failed) has only some of its batches: wait for Module C's tables.
  [[ -f "$FC/$lib/timeseries.tsv" ]] || { echo "error: $lib replay not finished (no $FC/$lib/timeseries.tsv)"; exit 1; }
  n=$(ls -d "$FC/$lib/pipeline/$lib/batches"/batch_* | wc -l | tr -d ' ')
  local pts; pts=$(points "$n")
  echo "=== $lib: $n batches, points: $pts ==="
  local k
  for k in $pts; do assemble_point "$lib" "$ref" "$k" false auto "$mincov"; done
  [[ "$POLISH_END" == true && "$model" != none ]] && assemble_point "$lib" "$ref" "$n" true "$model" "$mincov"
  rmdir "$FC/$lib/coverage_over_time/.stage" 2>/dev/null
  LIBS+=("$lib"); REFS+=("$ref")
}

LIBS=(); REFS=()
[[ "$RUN_ORBI"      == true ]] && run_lib Orbivirus_kevo.refseq_abv1g_k11w5 "$ORBI_REF" "$ORBI_MEDAKA_MODEL"
[[ "$RUN_TBEV"      == true ]] && run_lib TBEV_FBC92293.refseq_abv1g "$TBEV_REF" "$TBEV_MEDAKA_MODEL"
[[ "$RUN_TICKS"     == true ]] && run_lib ticks_35-33.refseq_abv1g_k11w5 "$TICK_CCHFV_REF" none "$TICK_MIN_COVERAGE"
[[ "$RUN_TICKS"     == true ]] && run_lib ticks_58-37.refseq_abv1g_k11w5 "$TICK_58_REF" none "$TICK_MIN_COVERAGE"

# The scorer takes one reference per call: score each library on its own, then
# rebuild $METRICS from the new rows plus the old rows of libraries not run here.
if [[ "$SCORE" == true && ${#LIBS[@]} -gt 0 ]]; then
  mkdir -p "$(dirname "$METRICS")"
  : > "$LOG_DIR/score.log"
  rc=0
  for i in "${!LIBS[@]}"; do
    "$NB_PY" jobs/framework_C/framework_C.coverage_over_time.score.py "${REFS[$i]}" "${LIBS[$i]}" \
      > "$METRICS.${LIBS[$i]}.tmp" 2>> "$LOG_DIR/score.log" || { rc=1; echo "[FAIL] scoring ${LIBS[$i]} (see $LOG_DIR/score.log)"; }
  done
  if [[ $rc -eq 0 ]]; then
    {
      head -1 "$METRICS.${LIBS[0]}.tmp"
      [[ -f "$METRICS" ]] && awk -F'\t' -v libs="${LIBS[*]}" \
        'BEGIN { n = split(libs, a, " "); for (i = 1; i <= n; i++) skip[a[i]] = 1 }
         NR > 1 && !($1 in skip)' "$METRICS"
      for lib in "${LIBS[@]}"; do tail -n +2 "$METRICS.$lib.tmp"; done
    } > "$METRICS.tmp" && mv "$METRICS.tmp" "$METRICS" && echo "wrote $METRICS"
  fi
  rm -f "$METRICS".*.tmp
fi
