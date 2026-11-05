#!/usr/bin/env bash
# ============================================================================
# Framework B — Assembly quality metrics for EMtrace SRA consensus results
# ============================================================================
#
# For each Framework B SRA results subdir, builds a manifest of
#   sample_id, consensus, reference, stats_div
# and scores it with the ngs-benchmark Module B CLI:
#   ngs-benchmark assembly-batch --manifest M --outdir O
#
# Metrics (ROADMAP Module B): genome coverage (%), identity (%), and
# positions_called (%) from consensus.fa aligned to the fetched reference
# (global pairwise MAFFT alignment; identity excludes consensus N positions).
# Reference: assembly/refaccs.fetched.fa cached by emtrace in each sample dir.
#   coverage_pct          = covered_bp / ref_len    × 100   (breadth)
#   identity_pct          = matching_bp / aligned_bp × 100  (accuracy within aligned)
#   positions_called_pct  = called_bp / ref_len × 100
#                           called_bp = non-N bases in consensus.fa (ref-length).
#                           Read directly from the delivered consensus (the
#                           authoritative output), not from an intermediate mask.
# viral_reads is the read_count for division=Viruses in
# results/reports/stats_div.csv ("no hits" when the file is absent).
#
# Pipeline map — edit PIPELINE_DIRS to add/rename batches:
#   sra_petersen_2026  Petersen 2026 SRA (ZIKV, WNV, SARS-CoV-2, DENV-2)
#   sra_gauthier_2021  Gauthier 2021 SRA (SARS-CoV-2, Ct 14.7–30.6 + 5 NC)
#   sra_petersen_2026_v20260924_{refauto,chain2}
#                      Petersen 2026, v20260924 re-run with --ref_auto
#                      (framework_B.emtrace.v20260924.bash); scored against
#                      the fixed references in REF_DIRS
#
# Results layout (three flat TSVs per batch <label>, no per-sample dirs):
#   (<label> gets a .d<N> suffix when MIN_COVERAGE is set — see below)
#   metrics/framework_B/<label>.manifest.tsv  input manifest (kept for debugging)
#   metrics/framework_B/<label>.summary.tsv   one row per sample
#   metrics/framework_B/<label>.coverage.tsv  one row per segment per sample
#
# <label>.summary.tsv columns:
#   sample, n_segments, total_ref_len, mean_coverage_pct, mean_identity_pct,
#   weighted_coverage_pct, weighted_identity_pct,
#   mean_positions_called_pct, weighted_positions_called_pct, viral_reads
# <label>.coverage.tsv columns:
#   sample, segment, ref_len, covered_bp, aligned_bp, matching_bp,
#   coverage_pct, identity_pct, called_bp, positions_called_pct
#
# Samples whose refaccs.fetched.fa is absent are skipped (run not yet complete).
# Missing consensus.fa → NA alignment metrics and positions_called_pct.
# Missing stats_div.csv → viral_reads = "no hits".

# =============================================================================
# Pipeline map — label → results subdir
# =============================================================================
declare -A PIPELINE_DIRS=(
  [sra_petersen_2026]="results/framework_b/sra_petersen_2026"
  [sra_gauthier_2021]="results/framework_b/sra_gauthier_2021"
  [sra_petersen_2026_k11w5]="results/framework_b/sra_petersen_2026_k11w5s100"
  [sra_petersen_2026_v20260924_refauto]="results/framework_b/sra_petersen_2026_v20260924_refauto"
  [sra_petersen_2026_v20260924_chain2]="results/framework_b/sra_petersen_2026_v20260924_chain2"
)

# Fixed scoring reference per label (optional): <dir>/<sample>.fasta. The
# v20260924 batches choose their own reference (--ref_auto), so scoring against
# their own refaccs.fetched.fa would be circular. Petersen is scored against the
# isolate-matched complete genomes in data/sra/refs/petersen_2026/, the same
# files ViMOP and TELEVIR are scored against (rationale: Supplementary Table 5). Consensus records are matched by accession, else the sole record is used.
declare -A REF_DIRS=(
  [sra_petersen_2026_v20260924_refauto]="data/sra/refs/petersen_2026"
  [sra_petersen_2026_v20260924_chain2]="data/sra/refs/petersen_2026"
)

# Score only the primary consensus record (labels listed here). ref_auto may
# assemble two accessions of one virus (Petersen Zika: NC_035889.1 with 62.8k
# reads and NC_012532.1 with 353), while the fixed reference has one record.
# Name matching would then score the minor record that happens to share the
# reference's accession; instead the selected accession with the most reads in
# assembly/ref_selection.tsv is extracted to <label>.primary/<sample>.fa and
# scored alone, as ViMOP's and TELEVIR's single consensus are.
declare -A PRIMARY_RECORD=(
  [sra_petersen_2026_v20260924_refauto]=1
  [sra_petersen_2026_v20260924_chain2]=1
)

# =============================================================================
# Control switches — set to true/false before running
# =============================================================================
RUN_PETERSEN=false  # Petersen 2026: ZIKV, WNV, SARS-CoV-2, DENV-2
RUN_GAUTHIER=false   # Gauthier 2021: SARS-CoV-2 Ct 14.7–30.6 + 5 NC
RUN_PETERSEN_K11W5=false  # Petersen 2026, k11/w5 sensitive annotation config
RUN_PETERSEN_V20260924_REFAUTO=true  # Petersen 2026, v20260924 re-run, CELL refauto
RUN_PETERSEN_V20260924_CHAIN2=true   # Petersen 2026, v20260924 re-run, CELL chain2
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
  local ref_dir="${REF_DIRS[$1]:-}"
  local primary="${PRIMARY_RECORD[$1]:-}"
  [[ -n "$primary" ]] && mkdir -p "$METRICS_DIR/${label}.primary"

  if [[ ! -d "$results_dir" ]]; then
    echo "[SKIP] $label: results dir not found: $results_dir" >&2
    return 1
  fi

  local manifest="$METRICS_DIR/${label}.manifest.tsv"
  printf 'sample_id\tconsensus\treference\tstats_div\n' > "$manifest"

  local n=0 n_skip=0
  for sample_dir in "$results_dir"/*/; do
    [[ -d "$sample_dir" ]] || continue
    local sample
    sample=$(basename "$sample_dir")

    local reference="$sample_dir/$ASSEMBLY_DIR/refaccs.fetched.fa"
    [[ -n "$ref_dir" ]] && reference="$ref_dir/$sample.fasta"

    local consensus="$sample_dir/$ASSEMBLY_DIR/consensus.fa"
    local stats_div="$sample_dir/results/reports/stats_div.csv"

    local selection="$sample_dir/$ASSEMBLY_DIR/ref_selection.tsv"
    if [[ -n "$primary" && -f "$consensus" && -f "$selection" ]]; then
      local acc
      acc=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}
                         $c["selected"]=="True" && $c["reads"]>m{m=$c["reads"];a=$c["accession"]}
                         END{print a}' "$selection")
      local primary_fa="$METRICS_DIR/${label}.primary/${sample}.fa"
      awk -v a="$acc" '/^>/{p=(substr($1,2)==a)}p' "$consensus" > "$primary_fa"
      echo "  [primary] $sample: $acc"
      consensus="$primary_fa"
    fi

    if [[ ! -f "$reference" ]]; then
      echo "  [skip] $sample: refaccs.fetched.fa not found (run incomplete?)" >&2
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
[[ "$RUN_PETERSEN" == true ]] && run_batch "sra_petersen_2026" "${PIPELINE_DIRS[sra_petersen_2026]}"
[[ "$RUN_GAUTHIER" == true ]] && run_batch "sra_gauthier_2021" "${PIPELINE_DIRS[sra_gauthier_2021]}"
[[ "$RUN_PETERSEN_K11W5" == true ]] && run_batch "sra_petersen_2026_k11w5" "${PIPELINE_DIRS[sra_petersen_2026_k11w5]}"
[[ "$RUN_PETERSEN_V20260924_REFAUTO" == true ]] && run_batch "sra_petersen_2026_v20260924_refauto" "${PIPELINE_DIRS[sra_petersen_2026_v20260924_refauto]}"
[[ "$RUN_PETERSEN_V20260924_CHAIN2" == true ]] && run_batch "sra_petersen_2026_v20260924_chain2" "${PIPELINE_DIRS[sra_petersen_2026_v20260924_chain2]}"

echo "metrics written to: $METRICS_DIR"
