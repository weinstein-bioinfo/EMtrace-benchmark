#!/usr/bin/env bash
# ============================================================================
# Framework B — Assembly quality metrics for INSaFLU/TELEVIR results
# ============================================================================
#
# Peer of framework_B.metrics.bash (EMtrace) and framework_B.vimop.metrics.bash
# (ViMOP). TELEVIR produces mapping output but NO consensus FASTA, so this script
# has TWO parts:
#
#   PART 1 — consensus creation
#     For every sample pulled from the INSaFLU web app under
#       results/insaflu_out/<set>/<sample>/
#     rebuild the consensus from the mapped-reads BAM with `samtools consensus`,
#     writing <sample>/televir.consensus.fa. Samples with no BAM/BAI are skipped
#     — these are the divergent-LASV runs where TELEVIR produced no target
#     mapping (empty report, or target absent among background hits), i.e. a
#     genuine "no consensus" outcome.
#
#   PART 2 — scoring
#     Score each consensus against the SAME benchmark reference EMtrace/ViMOP
#     used, so the three pipelines are directly comparable (joinable summaries):
#       sra    : fixed reference data/sra/refs/<set>/<acc>.fasta
#                (fallback: <emtrace_batch>/<acc>/assembly/refaccs.fetched.fa)
#       mutant : simulated mutated template
#                data/petersen_2026/mutants/<virus>/<base>.fasta
#                (base = sample minus _mixed; virus = base minus trailing _NN)
#     via `ngs-benchmark assembly-batch` (MAFFT global align; identity excludes
#     consensus N; positions_called = non-N / ref_len).
#
# Outputs (per label) to metrics/framework_B/:
#   televir_<label>.manifest.tsv   sample_id, consensus, reference
#   televir_<label>.summary.tsv    one row per sample   (consensus vs benchmark ref)
#   televir_<label>.coverage.tsv   one row per segment per sample
# Labels match ViMOP (peterson_2026, gauthier_2021, pure, mixed) so the tables
# join with vimop_<label>.* and sra_<label>.* / B*_*.summary.tsv on sample_id.
#
# NOTE — this is the rigorous, apples-to-apples leg (consensus vs benchmark
# reference). The companion extractor (same dir) extract_televir_metrics.py emits the
# platform's SELF-reported statistics (Cov % / 100-Mismatch rate vs TELEVIR's
# own auto-selected reference) to televir_<set>.summary.tsv (pure_mutants,
# mixed_mutants, petersen2026, gauthier2021). Different filenames — no collision.
#
# Samples with no consensus (no BAM) or no resolvable reference are skipped
# with a warning.
# ============================================================================
set -uo pipefail

REPO_ROOT="$(pwd)"
cd "$REPO_ROOT" || { echo "cannot cd to $REPO_ROOT"; exit 1; }

METRICS_DIR="metrics/framework_B"
INSAFLU_ROOT="results/insaflu_out"
MUTANTS_DIR="data/petersen_2026/mutants"
mkdir -p "$METRICS_DIR"

# TELEVIR delivers only the mapped-reads BAM and its reference, with no
# consensus and no depth threshold, so the consensus rules are ours. They match
# EMtrace as benchmarked (--min_coverage 10, dominant allele > 50 %):
#   samtools consensus -m simple -d MIN_DEPTH -c CALL_FRACT
# Simple mode is required: the default Bayesian mode ignores --min-depth, which
# made the pre-2026-09-28 rebuild call bases at 1-9x (SARS-CoV-2 SRR15356294:
# 97.8 % coverage from a BAM with 0 % of positions at >= 20x).
# ViMOP is scored on its own consensus (built-in 20x).
MIN_DEPTH=10
CALL_FRACT=0.5

# =============================================================================
# Control switches — set to true/false before running
# =============================================================================
RUN_PETERSEN=false   # Petersen 2026 SRA (4) — rescored 2026-09-28 on isolate-matched refs
RUN_GAUTHIER=true    # Gauthier 2021 SRA (16)
RUN_PURE=false       # pure simulated mutants (28: 4 viruses x 7 rates)
RUN_MIXED=false      # mixed simulated mutants (28: 4 viruses x 7 rates)
# =============================================================================

# =============================================================================
# PART 1 — consensus creation
# =============================================================================
# make_consensus SET_DIR   (SET_DIR under $INSAFLU_ROOT, e.g. petersen2026)
make_consensus() {
  local set_dir="$1"
  local root="$INSAFLU_ROOT/$set_dir"
  if [[ ! -d "$root" ]]; then
    echo "[SKIP consensus] $set_dir: dir not found: $root" >&2
    return 0
  fi
  local n=0 n_skip=0
  for d in "$root"/*/; do
    [[ -d "$d" ]] || continue
    local s; s=$(basename "$d")
    local bam bai
    bam=$(find "$d" -maxdepth 1 -name '*.bam' | head -1)
    bai=$(find "$d" -maxdepth 1 -name '*.bai' | head -1)
    if [[ -z "$bam" || -z "$bai" ]]; then
      echo "  [skip] $set_dir/$s: no bam/bai — no target consensus" >&2
      (( n_skip++ )); continue
    fi
    local out="$d/televir.consensus.fa"
    if samtools consensus -a -f fasta -m simple -d "$MIN_DEPTH" -c "$CALL_FRACT" "$bam" -o "$out"; then
      (( n++ ))
    else
      echo "  [err]  $set_dir/$s: samtools consensus failed" >&2
      rm -f "$out"
    fi
  done
  echo "[consensus $set_dir] $n written, $n_skip skipped (simple, min-depth $MIN_DEPTH, call-fract $CALL_FRACT)"
}

# =============================================================================
# PART 2 — scoring
# =============================================================================
# resolve_reference SET_TYPE SAMPLE ACCESSION EMTRACE_DIR  -> echoes ref path or empty
resolve_reference() {
  local set_type="$1" sample="$2" acc="$3" emtrace_dir="${4:-}"
  if [[ "$set_type" == mutant ]]; then
    local base="${sample%_mixed}"
    local virus="${base%_[0-9][0-9]}"
    echo "$MUTANTS_DIR/${virus}/${base}.fasta"
  else  # sra: <dir>/<acc>.fasta, else parallel EMtrace fetched reference (keyed by accession)
    if [[ -f "$emtrace_dir/$acc.fasta" ]]; then echo "$emtrace_dir/$acc.fasta"
    else echo "$emtrace_dir/$acc/assembly/refaccs.fetched.fa"; fi
  fi
}

# run_batch LABEL SET_TYPE SET_DIR [EMTRACE_DIR]
run_batch() {
  local label="$1" set_type="$2" set_dir="$3" emtrace_dir="${4:-}"
  local root="$INSAFLU_ROOT/$set_dir"

  if [[ ! -d "$root" ]]; then
    echo "[SKIP] $label: dir not found: $root" >&2
    return 1
  fi

  local manifest="$METRICS_DIR/televir_${label}.manifest.tsv"
  printf 'sample_id\tconsensus\treference\n' > "$manifest"

  local n=0 n_skip=0
  for d in "$root"/*/; do
    [[ -d "$d" ]] || continue
    local s; s=$(basename "$d")

    local consensus="$d/televir.consensus.fa"
    if [[ ! -s "$consensus" ]]; then
      echo "  [skip] $s: no consensus (no bam / not built)" >&2
      (( n_skip++ )); continue
    fi

    # sample_id joins with vimop/emtrace: accession for sra, dir name for mutants
    local sample_id acc=""
    if [[ "$set_type" == mutant ]]; then
      sample_id="$s"
    else
      acc="${s#*_}"        # strip leading virus/patient label -> accession
      sample_id="$acc"
    fi

    local reference; reference=$(resolve_reference "$set_type" "$s" "$acc" "$emtrace_dir")
    if [[ -z "$reference" || ! -f "$reference" ]]; then
      echo "  [skip] $s: reference not found (looked for ${reference:-<none>})" >&2
      (( n_skip++ )); continue
    fi

    printf '%s\t%s\t%s\n' "$sample_id" "$consensus" "$reference" >> "$manifest"
    (( n++ ))
  done

  if [[ "$n" -eq 0 ]]; then
    echo "[$label] no scorable samples in $root ($n_skip skipped)" >&2
    return 1
  fi

  echo "[$label] $n samples → ngs-benchmark assembly-batch ($n_skip skipped)"
  ngs-benchmark assembly-batch --manifest "$manifest" --outdir "$METRICS_DIR"
  mv "$METRICS_DIR/summary.tsv"  "$METRICS_DIR/televir_${label}.summary.tsv"
  mv "$METRICS_DIR/coverage.tsv" "$METRICS_DIR/televir_${label}.coverage.tsv"
  echo "[$label] done → $METRICS_DIR/televir_${label}.{manifest,summary,coverage}.tsv"
}

# =============================================================================
# Run
# =============================================================================
echo "=== PART 1: consensus creation (samtools consensus -m simple, min-depth $MIN_DEPTH, call-fract $CALL_FRACT) ==="
[[ "$RUN_PETERSEN" == true ]] && make_consensus petersen2026
[[ "$RUN_GAUTHIER" == true ]] && make_consensus gauthier2021
[[ "$RUN_PURE"     == true ]] && make_consensus pure_mutants
[[ "$RUN_MIXED"    == true ]] && make_consensus mixed_mutants

echo ""
echo "=== PART 2: scoring (ngs-benchmark assembly-batch vs benchmark reference) ==="
# SRA sets are scored against fixed, pipeline-independent references,
# data/sra/refs/<set>/<acc>.fasta, passed per batch (SRR15356294 is in
# both sets, with different references):
#   petersen_2026  isolate-matched complete genomes (Supplementary Table 5)
#   gauthier_2021  OX637002.1 for all 21 libraries, the reference of the June
#                  2026 EMtrace run, staged here so the scoring no longer
#                  depends on that run's folder. To be revisited (NC_045512.2)
#                  with the scheduled Gauthier re-run (the project notes), for all three pipelines together.
[[ "$RUN_PETERSEN" == true ]] && run_batch peterson_2026 sra    petersen2026  "data/sra/refs/petersen_2026"
[[ "$RUN_GAUTHIER" == true ]] && run_batch gauthier_2021 sra    gauthier2021  "data/sra/refs/gauthier_2021"
[[ "$RUN_PURE"     == true ]] && run_batch pure          mutant pure_mutants
[[ "$RUN_MIXED"    == true ]] && run_batch mixed         mutant mixed_mutants

# combined tier-2 summary across all scored labels — mirrors the tier-1
# televir_reported.all.summary.tsv, with a leading `label` column.
combined="$METRICS_DIR/televir_all.summary.tsv"
first=true
for label in peterson_2026 gauthier_2021 pure mixed; do
  f="$METRICS_DIR/televir_${label}.summary.tsv"
  [[ -s "$f" ]] || continue
  if $first; then { printf 'label\t'; head -1 "$f"; } > "$combined"; first=false; fi
  awk -F'\t' -v L="$label" 'NR>1{print L"\t"$0}' "$f" >> "$combined"
done
[[ -s "$combined" ]] && echo "combined tier-2 summary → $combined"

echo ""
echo "metrics written to: $METRICS_DIR"
