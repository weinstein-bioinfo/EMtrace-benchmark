#!/usr/bin/env bash
# ============================================================================
# Framework B — Assembly quality metrics for EMtrace consensus results
# ============================================================================
#
# For each Framework B results subdir, builds a manifest of
#   sample_id, consensus, reference, stats_div
# and scores it with the ngs-benchmark Module B CLI:
#   ngs-benchmark assembly-batch --manifest M --outdir O
#
# Metrics (ROADMAP Module B): genome coverage (%) and identity (%) from
# consensus.fa aligned to the simulated mutant genome (global pairwise MAFFT;
# identity excludes consensus N positions).
# Reference: data/petersen_2026/mutants/<virus>/<sample_id>.fasta
#   virus subdir = sample_id with _mixed and trailing _NN deviation stripped.
#   coverage_pct = covered_bp / ref_len    × 100   (breadth over the reference)
#   identity_pct = matching_bp / aligned_bp × 100   (accuracy of called bases)
# viral_reads is the read_count for division=Viruses in
# results/reports/stats_div.csv ("no hits" when the file is absent).
#
# Pipeline map — edit PIPELINE_DIRS to add/rename batches:
#   B6_pure         pure viral reads (4 viruses × 7 mutation rates, default config)
#   B7_mixed        background-embedded reads (same matrix, default config)
#   B6_pure_k11w5   pure viral reads, k11/w5 minimap2 config
#   B7_mixed_k11w5  background-embedded reads, k11/w5 minimap2 config
#   B6_pure_v20260924, B7_mixed_v20260924
#                   full re-run v20260924 (framework_B.emtrace.v20260924.bash)
#   B6_pure_v20260924_k11w5, B7_mixed_v20260924_k11w5
#                   same, CELL=k11w5: sensitive annotation index
#                   minimap.refseq.abv.1g.k11w5; Lassa segments only
#   B6_pure_v20260924_chain2, B7_mixed_v20260924_chain2
#                   same, CELL=chain2: nucleotide+DIAMOND chain annotation,
#                   reference from assembly.database only
#   B6_pure_v20260924_refprim, B7_mixed_v20260924_refprim
#                   same, CELL=refprim: ref_auto also searches
#                   assembly.database at k11/w5 (use_assembly_db true)
#
# Results layout (three flat TSVs per batch <label>, no per-sample dirs):
#   (<label> gets a .d<N> suffix when MIN_COVERAGE is set — see below)
#   metrics/framework_B/<label>.manifest.tsv  input manifest (kept for debugging)
#   metrics/framework_B/<label>.summary.tsv   one row per sample
#   metrics/framework_B/<label>.coverage.tsv  one row per segment per sample
#
# <label>.summary.tsv columns:
#   sample, n_segments, total_ref_len, mean_coverage_pct, mean_identity_pct,
#   weighted_coverage_pct, weighted_identity_pct, viral_reads
# <label>.coverage.tsv columns:
#   sample, segment, ref_len, covered_bp, aligned_bp, matching_bp,
#   coverage_pct, identity_pct
#
# Samples whose mutant FASTA is absent are skipped (name doesn't match any known
# mutant). Missing consensus.fa → NA metrics; missing stats_div.csv → "no hits".

# =============================================================================
# Pipeline map — label → results subdir
# =============================================================================
declare -A PIPELINE_DIRS=(
  [B6_pure]="results/framework_b/pure"
  [B7_mixed]="results/framework_b/mixed"
  [B6_pure_k11w5]="results/framework_b/pure_k11w5s100"
  [B7_mixed_k11w5]="results/framework_b/mixed_k11w5s100"
  [B6_pure_v20260924]="results/framework_b/pure_v20260924_refauto"    # framework_B.emtrace.v20260924.bash
  [B7_mixed_v20260924]="results/framework_b/mixed_v20260924_refauto"
  [B6_pure_v20260924_refprim]="results/framework_b/pure_v20260924_refprim"    # REF_VARIANT=refprim
  [B7_mixed_v20260924_refprim]="results/framework_b/mixed_v20260924_refprim"
  [B6_pure_v20260924_k11w5]="results/framework_b/pure_v20260924_k11w5"      # CELL=k11w5
  [B7_mixed_v20260924_k11w5]="results/framework_b/mixed_v20260924_k11w5"
  [B6_pure_v20260924_chain2]="results/framework_b/pure_v20260924_chain2"      # CELL=chain2
  [B7_mixed_v20260924_chain2]="results/framework_b/mixed_v20260924_chain2"
)

# =============================================================================
# Control switches — set to true/false before running
# =============================================================================
RUN_B6=false          # B.6  pure viral reads          (28 samples: 4 viruses × 7 deviations)
RUN_B7=false          # B.7  background-mixed reads    (same matrix)
RUN_B6_K11W5=false   # B.6  pure, k11/w5 minimap2 config
RUN_B7_K11W5=false   # B.7  mixed, k11/w5 minimap2 config
RUN_B6_V20260924=false  # pure,  full re-run v20260924 (ref_auto, d10, indel calling)
RUN_B7_V20260924=false  # mixed, full re-run v20260924
RUN_B6_V20260924_REFPRIM=false  # pure,  v20260924 + primary ref search (k11/w5)
RUN_B7_V20260924_REFPRIM=false  # mixed, v20260924 + primary ref search (k11/w5)
RUN_B6_V20260924_K11W5=false     # pure,  v20260924 + sensitive annotation index (Lassa only)
RUN_B7_V20260924_K11W5=false     # mixed, v20260924 + sensitive annotation index (Lassa only)
RUN_B6_V20260924_CHAIN2=true    # pure,  v20260924 + nucleotide/DIAMOND chain
RUN_B7_V20260924_CHAIN2=true    # mixed, v20260924 + nucleotide/DIAMOND chain
# =============================================================================

# Per-threshold scoring: set to a value the re-assembly script was run with
# (framework_B.emtrace.assemblies.bash) to score <SAMPLE>/assembly.d<N>/
# instead of <SAMPLE>/assembly/, writing <label>.d<N>.{manifest,summary,coverage}.tsv.
# Empty = score assembly/ into the unsuffixed <label>.* files (as before).
# METRICS_TAG (env, optional) is appended to the label, e.g. METRICS_TAG=ontsup →
# <label>.d10.ontsup.*, to score a new code version without overwriting older scores.
MIN_COVERAGE="${MIN_COVERAGE:-}"   # or: MIN_COVERAGE=10 bash <this script>
ASSEMBLY_DIR="assembly${MIN_COVERAGE:+.d$MIN_COVERAGE}"
LABEL_SUFFIX="${MIN_COVERAGE:+.d$MIN_COVERAGE}${METRICS_TAG:+.$METRICS_TAG}"

METRICS_DIR="metrics/framework_B"
mkdir -p "$METRICS_DIR"


# -----------------------------------------------------------------------------
# run_batch LABEL RESULTS_DIR
#   Builds a manifest from the sample subdirs, runs ngs-benchmark
#   assembly-batch over it, writing three files to $METRICS_DIR/<label>.*:
#     <label>.manifest.tsv   kept for debugging
#     <label>.summary.tsv    one row per sample
#     <label>.coverage.tsv   one row per segment per sample
# -----------------------------------------------------------------------------
run_batch() {
  local label="$1$LABEL_SUFFIX"
  local results_dir="$2"

  if [[ ! -d "$results_dir" ]]; then
    echo "[SKIP] $label: results dir not found: $results_dir" >&2
    return 1
  fi

  local manifest="$METRICS_DIR/${label}.manifest.tsv"
  printf 'sample_id\tconsensus\treference\tstats_div\n' > "$manifest"

  local mutants_dir="data/petersen_2026/mutants"

  local n=0 n_skip=0
  for sample_dir in "$results_dir"/*/; do
    [[ -d "$sample_dir" ]] || continue
    local sample
    sample=$(basename "$sample_dir")

    # Derive virus subdir: strip _mixed suffix, then trailing _NN deviation.
    local base="${sample%_mixed}"
    local virus="${base%_[0-9][0-9]}"
    local reference="$mutants_dir/${virus}/${base}.fasta"

    local consensus="$sample_dir/$ASSEMBLY_DIR/consensus.fa"
    local stats_div="$sample_dir/results/reports/stats_div.csv"

    if [[ ! -f "$reference" ]]; then
      echo "  [skip] $sample: mutant FASTA not found: $reference" >&2
      (( n_skip++ ))
      continue
    fi

    printf '%s\t%s\t%s\t%s\n' "$sample" "$consensus" "$reference" "$stats_div" >> "$manifest"
    (( n++ ))
  done

  if [[ "$n" -eq 0 ]]; then
    echo "[${label}] no scorable samples in $results_dir ($n_skip skipped)" >&2
    return 1
  fi

  echo "[${label}] $n samples → ngs-benchmark assembly-batch ($n_skip skipped)"
  ngs-benchmark assembly-batch \
    --manifest "$manifest" \
    --outdir   "$METRICS_DIR"

  # ngs-benchmark writes summary.tsv / coverage.tsv; rename to labelled files
  mv "$METRICS_DIR/summary.tsv"  "$METRICS_DIR/${label}.summary.tsv"
  mv "$METRICS_DIR/coverage.tsv" "$METRICS_DIR/${label}.coverage.tsv"

  echo "[${label}] done → $METRICS_DIR/${label}.{manifest,summary,coverage}.tsv"
}


# =============================================================================
# Run
# =============================================================================
[[ "$RUN_B6"       == true ]] && run_batch "B6_pure"         "${PIPELINE_DIRS[B6_pure]}"
[[ "$RUN_B7"       == true ]] && run_batch "B7_mixed"        "${PIPELINE_DIRS[B7_mixed]}"
[[ "$RUN_B6_K11W5" == true ]] && run_batch "B6_pure_k11w5"  "${PIPELINE_DIRS[B6_pure_k11w5]}"
[[ "$RUN_B7_K11W5" == true ]] && run_batch "B7_mixed_k11w5" "${PIPELINE_DIRS[B7_mixed_k11w5]}"
[[ "$RUN_B6_V20260924" == true ]] && run_batch "B6_pure_v20260924"  "${PIPELINE_DIRS[B6_pure_v20260924]}"
[[ "$RUN_B7_V20260924" == true ]] && run_batch "B7_mixed_v20260924" "${PIPELINE_DIRS[B7_mixed_v20260924]}"
[[ "$RUN_B6_V20260924_REFPRIM" == true ]] && run_batch "B6_pure_v20260924_refprim"  "${PIPELINE_DIRS[B6_pure_v20260924_refprim]}"
[[ "$RUN_B7_V20260924_REFPRIM" == true ]] && run_batch "B7_mixed_v20260924_refprim" "${PIPELINE_DIRS[B7_mixed_v20260924_refprim]}"
[[ "$RUN_B6_V20260924_K11W5" == true ]] && run_batch "B6_pure_v20260924_k11w5"  "${PIPELINE_DIRS[B6_pure_v20260924_k11w5]}"
[[ "$RUN_B7_V20260924_K11W5" == true ]] && run_batch "B7_mixed_v20260924_k11w5" "${PIPELINE_DIRS[B7_mixed_v20260924_k11w5]}"
[[ "$RUN_B6_V20260924_CHAIN2" == true ]] && run_batch "B6_pure_v20260924_chain2"  "${PIPELINE_DIRS[B6_pure_v20260924_chain2]}"
[[ "$RUN_B7_V20260924_CHAIN2" == true ]] && run_batch "B7_mixed_v20260924_chain2" "${PIPELINE_DIRS[B7_mixed_v20260924_chain2]}"

echo "metrics written to: $METRICS_DIR"
