#!/usr/bin/env bash
# ============================================================================
# Framework B — score every Supplementary Table 3 / Supplementary Fig. 1a row with ngs-benchmark
# ============================================================================
# One scorer for all rows, `ngs-benchmark assembly-batch` (Module B: MAFFT
# G-INS-i, --adjustdirection, one pair at a time, one thread), which now also
# counts insertions and deletions (ngs_benchmark/assembly/mafft.py::indel_counts).
# (Replaced a standalone scorer, framework_B.medaka.score.py, removed 29 Sep 2026.)
#
# Rows (labels fwB.<set>.<row>):
#   EMtrace             v20260924 refauto, medaka-polished   <sample>/assembly.medaka/consensus.fa
#   EMtrace.nopolish    v20260924 refauto, unpolished        <sample>/assembly/consensus.fa
#   EMtrace.chain2      v20260924 chain2,  medaka-polished
#   EMtrace.chain2.nopolish
#   ViMOP, TELEVIR         the consensus files of their metrics manifests
# Sets: pure, mixed (simulated mutants; reference = the mutant truth) and
# petersen (field libraries; reference = the isolate-matched genome, RESULTS
# note [10]; EMtrace scored on its primary record, as framework_B.sra.metrics.bash)
# and gauthier (16 SARS-CoV-2 patient libraries; reference = Wuhan-Hu-1
# NC_045512.2 for every library, 9 Oct 2026: the libraries span several lineages
# (2020 early lineages, 2021 Alpha and N501Y+E484K) and none has an isolate genome
# in GenBank; identity is identity to Wuhan-Hu-1. The earlier OX637002.1 refs and
# scores are kept in refs/gauthier_2021.OX637002.1 and metrics/fwB.gauthier.OX637002.1).
# A sample without a consensus is listed with a missing path: ngs-benchmark
# reports NA metrics, which the tables count as coverage 0.
#
# Output:  metrics/framework_B/fwB.<set>.<row>.{summary,coverage}.tsv (+ .manifest.tsv)
# Read by: rlib/bmark_io.R::load_rescored / load_rescored_petersen
#
# Memory: G-INS-i on a 30 kb genome needs ~6.5 GB; assembly-batch aligns one
# pair at a time, so never run two of these scripts at once on a 16 GB machine.
#
# Run from the repository root:
#   bash jobs/framework_B/framework_B.rescore.bash
#   RUN_MUTANTS=false RUN_PETERSEN=false bash ...   (Gauthier only)
# ============================================================================

set -euo pipefail
[[ -d jobs ]] || { echo "error: run from the repository root" >&2; exit 1; }

M=metrics/framework_B
FB=results/framework_b
TRUTH=data/petersen_2026/mutants
RUN_MUTANTS="${RUN_MUTANTS:-true}"
RUN_PETERSEN="${RUN_PETERSEN:-true}"
RUN_GAUTHIER="${RUN_GAUTHIER:-true}"
ROWS_RE="${ROWS_RE:-.}"   # only score labels matching this regex, e.g. ROWS_RE='ViMOP|TELEVIR'

# score <label> <manifest>: assembly-batch, then name its outputs after the label
score() {
  local label="$1" manifest="$2" out="$M/.fwB.$1"
  [[ "$label" =~ $ROWS_RE ]] || return 0
  echo "[score] $label ($(( $(wc -l < "$manifest") - 1 )) samples)"
  ngs-benchmark assembly-batch --manifest "$manifest" --outdir "$out" > "$M/$label.log" 2>&1
  mv "$out/summary.tsv" "$M/$label.summary.tsv"
  mv "$out/coverage.tsv" "$M/$label.coverage.tsv"
  cp "$manifest" "$M/$label.manifest.tsv"
  rmdir "$out" 2>/dev/null || true
}

# mutant_manifest <set> <cell> <subdir>: the 28 mutants of one EMtrace cell
mutant_manifest() {
  local set="$1" cell="$2" sub="$3" out; out=$(mktemp)
  printf 'sample_id\tconsensus\treference\n' > "$out"
  local v dev s dir
  for v in COVID DENV LASV_L LASV_S; do for dev in 00 05 10 15 20 25 30; do
    s="${v}_${dev}"; [[ "$set" == mixed ]] && s="${s}_mixed"
    dir="$FB/${set}_v20260924_${cell}/$s"
    printf '%s\t%s\t%s\n' "$s" "$dir/$sub/consensus.fa" "$TRUTH/$v/${v}_${dev}.fasta" >> "$out"
  done; done
  echo "$out"
}

# Gauthier 2021 patient libraries (Ct 14.7-30.6); the 5 negative controls are not scored
GAUTHIER_PATIENTS=(
  SRR15356303 SRR15356328 SRR15356305 SRR15356306 SRR15356327 SRR15356326
  SRR15356325 SRR15356329 SRR15356307 SRR15356298 SRR15356308 SRR15356322
  SRR15356310 SRR15356311 SRR15356321 SRR15356320
)

# gauthier_manifest <subdir>: the v20260924 refauto run. A library without
# <subdir>/consensus.fa falls back to assembly/: medaka skips a library with no
# position at min_coverage (P34), so the pipeline's output is the unpolished one.
gauthier_manifest() {
  local sub="$1" out s dir cons; out=$(mktemp)
  printf 'sample_id\tconsensus\treference\n' > "$out"
  for s in "${GAUTHIER_PATIENTS[@]}"; do
    dir="$FB/sra_gauthier_2021_v20260924_refauto/$s"
    cons="$dir/$sub/consensus.fa"; [[ -f "$cons" ]] || cons="$dir/assembly/consensus.fa"
    printf '%s\t%s\t%s\n' "$s" "$cons" "data/sra/refs/gauthier_2021/$s.fasta" >> "$out"
  done
  echo "$out"
}

# patients_only <manifest>: keep the header and the GAUTHIER_PATIENTS rows, with
# the reference set to refs/gauthier_2021 (the ViMOP manifest named the old
# run's fetched OX637002.1, so every row is scored against the same genome)
patients_only() {
  local out; out=$(mktemp)
  awk -F'\t' -v OFS='\t' -v ids="${GAUTHIER_PATIENTS[*]}" \
    'BEGIN{n=split(ids,a," ");for(i=1;i<=n;i++)k[a[i]]=1}
     NR==1{print;next} ($1 in k){$3="data/sra/refs/gauthier_2021/"$1".fasta";print}' "$1" > "$out"
  echo "$out"
}

if [[ "$RUN_MUTANTS" == true ]]; then
  for set in pure mixed; do
    score "fwB.$set.EMtrace"                 "$(mutant_manifest "$set" refauto assembly.medaka)"
    score "fwB.$set.EMtrace.nopolish"        "$(mutant_manifest "$set" refauto assembly)"
    score "fwB.$set.EMtrace.chain2"          "$(mutant_manifest "$set" chain2 assembly.medaka)"
    score "fwB.$set.EMtrace.chain2.nopolish" "$(mutant_manifest "$set" chain2 assembly)"
    score "fwB.$set.ViMOP"                      "$M/vimop_${set}.manifest.tsv"
    score "fwB.$set.TELEVIR"                    "$M/televir_${set}.manifest.tsv"
  done
fi

if [[ "$RUN_PETERSEN" == true ]]; then
  # primary-record manifests are built by framework_B.medaka.petersen.score.bash
  # (polished) and framework_B.sra.metrics.bash (unpolished)
  score fwB.petersen.EMtrace          "$M/sra_petersen_2026_v20260924_refauto.medaka.manifest.tsv"
  score fwB.petersen.EMtrace.chain2   "$M/sra_petersen_2026_v20260924_chain2.medaka.manifest.tsv"
  score fwB.petersen.EMtrace.nopolish "$M/sra_petersen_2026_v20260924_refauto.manifest.tsv"
  score fwB.petersen.ViMOP               "$M/vimop_peterson_2026.manifest.tsv"
  score fwB.petersen.TELEVIR             "$M/televir_peterson_2026.manifest.tsv"
fi

if [[ "$RUN_GAUTHIER" == true ]]; then
  # 16 patient libraries; every row against Wuhan-Hu-1 NC_045512.2
  # (refs/gauthier_2021; see the header). ref_auto picked one accession
  # (NC_045512.2) per library, so no primary-record step.
  score fwB.gauthier.EMtrace          "$(gauthier_manifest assembly.medaka)"
  score fwB.gauthier.EMtrace.nopolish "$(gauthier_manifest assembly)"
  score fwB.gauthier.ViMOP               "$(patients_only "$M/vimop_gauthier_2021.manifest.tsv")"
  score fwB.gauthier.TELEVIR             "$(patients_only "$M/televir_gauthier_2021.manifest.tsv")"
fi
echo "done: $M/fwB.*.{summary,coverage}.tsv"
