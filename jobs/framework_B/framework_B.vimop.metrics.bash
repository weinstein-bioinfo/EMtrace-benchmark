#!/usr/bin/env bash
# ============================================================================
# Framework B — Assembly quality metrics for ViMOP consensus results
# ============================================================================
#
# Mirrors framework_B.sra.metrics.bash / framework_B.metrics.bash but scores
# ViMOP's per-organism deliverable, selected_consensus/*.fasta, against the
# SAME truth reference EMtrace used — so the two pipelines are directly
# comparable (joinable summary tables).
#
# Reference resolution (per set type):
#   sra    : the parallel EMtrace fetched reference
#            results/framework_b/<emtrace_batch>/<sample>/assembly/refaccs.fetched.fa
#            (fallback: data/truth/vimop_ref_map.tsv  sample_id<TAB>reference_path)
#   mutant : the simulated mutated template
#            data/petersen_2026/mutants/<virus>/<base>.fasta
#            (base = sample minus _mixed; virus = base minus trailing _NN)
#
# Consensus: selected_consensus/*.fasta (one record/file per organism; with the
# single-target --targets filter there is normally exactly one).
#
# Scores with the ngs-benchmark Module B CLI (MAFFT global align; identity
# excludes consensus N; positions_called = non-N / ref_len):
#   ngs-benchmark assembly-batch --manifest M --outdir O
#
# Outputs (per label) to metrics/framework_B/:
#   vimop_<label>.manifest.tsv   sample_id, consensus, reference
#   vimop_<label>.summary.tsv    one row per sample
#   vimop_<label>.coverage.tsv   one row per segment per sample
# Columns match the EMtrace sra_*/B*_*.summary.tsv schema (viral_reads is
# absent — ViMOP has no stats_div.csv; join on sample for the comparison).
#
# Samples with no selected_consensus/ (run incomplete) or no resolvable
# reference are skipped with a warning.
# ============================================================================
set -uo pipefail

REPO_ROOT="$(pwd)"
cd "$REPO_ROOT" || { echo "cannot cd to $REPO_ROOT"; exit 1; }

METRICS_DIR="metrics/framework_B"
VIMOP_ROOT="results/vimop_out"
MUTANTS_DIR="data/petersen_2026/mutants"
REF_MAP="data/truth/vimop_ref_map.tsv"
mkdir -p "$METRICS_DIR"

# =============================================================================
# Control switches — set to true/false before running
# =============================================================================
RUN_PETERSEN=true    # ViMOP Petersen (local runs; rescored 2026-09-28 on isolate-matched refs)
RUN_GAUTHIER=false   # ViMOP Gauthier (Puhti results synced to vimop_out/gauthier)
RUN_PURE=false       # ViMOP pure mutants (Puhti results synced to vimop_out/pure)
RUN_MIXED=false      # ViMOP mixed mutants (Puhti results synced to vimop_out/mixed)
# =============================================================================

# --- ref_map fallback (sample_id -> reference path) --------------------------
declare -A REF_OVERRIDE=()
if [[ -f "$REF_MAP" ]]; then
  while IFS=$'\t' read -r s r _; do
    [[ -z "$s" || "$s" == \#* || "$s" == sample_id ]] && continue
    REF_OVERRIDE["$s"]="$r"
  done < "$REF_MAP"
fi

# pick_consensus SAMPLE_DIR -> echoes the consensus FASTA to score (empty if none).
#   1) ViMOP's selected per-organism deliverable: selected_consensus/*.fasta
#   2) fallback (when ViMOP selected none, e.g. genome recovery below its threshold —
#      the SARS-CoV-2 SRR15356294 near-non-recovery case): the per-reference
#      consensus/*.consensus.fasta with the most called (non-N) bases.
pick_consensus() {
  local d="$1" best="" bestn=-1 f n
  best=$(find "$d/selected_consensus" -maxdepth 1 -name '*.fasta' 2>/dev/null | head -1)
  if [[ -n "$best" ]]; then echo "$best"; return; fi
  while IFS= read -r f; do
    [[ -n "$f" ]] || continue
    n=$(awk '!/^>/{s=s$0} END{gsub(/[Nn]/,"",s); print length(s)}' "$f")
    if (( n > bestn )); then bestn=$n; best=$f; fi
  done < <(find "$d/consensus" -maxdepth 1 -name '*.consensus.fasta' 2>/dev/null)
  echo "$best"
}

# resolve_reference SET_TYPE SAMPLE EMTRACE_DIR  -> echoes ref path or empty
resolve_reference() {
  local set_type="$1" sample="$2" emtrace_dir="${3:-}"
  if [[ "$set_type" == mutant ]]; then
    local base="${sample%_mixed}"
    local virus="${base%_[0-9][0-9]}"
    echo "$MUTANTS_DIR/${virus}/${base}.fasta"
  else  # sra: <dir>/<sample>.fasta, else parallel EMtrace fetched reference, then map override
    local auto="$emtrace_dir/$sample/assembly/refaccs.fetched.fa"
    if [[ -f "$emtrace_dir/$sample.fasta" ]]; then echo "$emtrace_dir/$sample.fasta"
    elif [[ -f "$auto" ]]; then echo "$auto"
    else echo "${REF_OVERRIDE[$sample]:-}"; fi
  fi
}

# run_batch LABEL SET_TYPE VIMOP_DIR [EMTRACE_DIR]
run_batch() {
  local label="$1" set_type="$2" vimop_dir="$3" emtrace_dir="${4:-}"

  if [[ ! -d "$vimop_dir" ]]; then
    echo "[SKIP] $label: ViMOP results dir not found: $vimop_dir" >&2
    return 1
  fi

  local manifest="$METRICS_DIR/vimop_${label}.manifest.tsv"
  printf 'sample_id\tconsensus\treference\n' > "$manifest"

  local n=0 n_skip=0
  for sample_dir in "$vimop_dir"/*/; do
    [[ -d "$sample_dir" ]] || continue
    local sample; sample=$(basename "$sample_dir")

    # ViMOP's per-organism deliverable, with per-reference fallback
    local consensus; consensus=$(pick_consensus "$sample_dir")
    if [[ -z "$consensus" ]]; then
      echo "  [skip] $sample: no selected_consensus/ and no consensus/*.consensus.fasta (run incomplete?)" >&2
      (( n_skip++ )); continue
    fi
    if [[ "$consensus" != *"/selected_consensus/"* ]]; then
      echo "  [note] $sample: ViMOP selected no consensus — fallback to best per-ref $(basename "$consensus")" >&2
    fi

    local reference; reference=$(resolve_reference "$set_type" "$sample" "$emtrace_dir")
    if [[ -z "$reference" || ! -f "$reference" ]]; then
      echo "  [skip] $sample: no truth reference (looked for ${reference:-<none>})" >&2
      (( n_skip++ )); continue
    fi

    printf '%s\t%s\t%s\n' "$sample" "$consensus" "$reference" >> "$manifest"
    (( n++ ))
  done

  if [[ "$n" -eq 0 ]]; then
    echo "[$label] no scorable samples in $vimop_dir ($n_skip skipped)" >&2
    return 1
  fi

  echo "[$label] $n samples → ngs-benchmark assembly-batch ($n_skip skipped)"
  ngs-benchmark assembly-batch --manifest "$manifest" --outdir "$METRICS_DIR"
  mv "$METRICS_DIR/summary.tsv"  "$METRICS_DIR/vimop_${label}.summary.tsv"
  mv "$METRICS_DIR/coverage.tsv" "$METRICS_DIR/vimop_${label}.coverage.tsv"
  echo "[$label] done → $METRICS_DIR/vimop_${label}.{manifest,summary,coverage}.tsv"
}

# =============================================================================
# Run
# =============================================================================
# Petersen 2026 is scored against the isolate-matched complete genomes in
# data/sra/refs/petersen_2026/<acc>.fasta (rationale: Supplementary Table 5); passed as the reference dir of that batch only, since
# SRR15356294 is also a Gauthier sample.
[[ "$RUN_PETERSEN" == true ]] && run_batch "peterson_2026" sra    "$VIMOP_ROOT/petersen" "data/sra/refs/petersen_2026"
[[ "$RUN_GAUTHIER" == true ]] && run_batch "gauthier_2021" sra    "$VIMOP_ROOT/gauthier" "results/framework_b/sra_gauthier_2021"
[[ "$RUN_PURE"     == true ]] && run_batch "pure"          mutant "$VIMOP_ROOT/pure"
[[ "$RUN_MIXED"    == true ]] && run_batch "mixed"         mutant "$VIMOP_ROOT/mixed"

echo "metrics written to: $METRICS_DIR"
