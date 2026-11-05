#!/usr/bin/env bash
# ============================================================================
# Framework B — Petersen 2026: primary-record manifests for the polished rows
# ============================================================================
# Builds the manifests of the medaka-polished EMtrace rows of Supplementary Fig. 1a (the
# unpolished row, ViMOP and TELEVIR already have theirs), scored against the
# fixed isolate-matched references (data/sra/refs/petersen_2026/<acc>.fasta;
# RESULTS note [10]). Scoring is done by framework_B.rescore.bash
# (ngs-benchmark assembly-batch), for all five rows:
#   EMtrace             sra_petersen_2026_v20260924_refauto/<acc>/assembly.medaka
#   EMtrace (chain2)    sra_petersen_2026_v20260924_chain2/<acc>/assembly.medaka
#   EMtrace (nopolish)  the refauto run's assembly/ (framework_B.sra.metrics.bash manifest)
#   ViMOP, TELEVIR         their metrics manifests (vimop_/televir_peterson_2026)
#
# Where ref_auto assembled two accessions of one virus (Zika, DENV), only the
# record with the most reads is scored, the rule of framework_B.sra.metrics.bash
# (PRIMARY_RECORD): the selected accession with the most reads in
# assembly/ref_selection.tsv. Polished records keep their accession names.
#
# Output: metrics/framework_B/sra_petersen_2026_v20260924_<cell>.medaka.{primary/,manifest.tsv}
# Then:   RUN_MUTANTS=false bash jobs/framework_B/framework_B.rescore.bash
# Run from the repository root:
#   bash jobs/framework_B/framework_B.medaka.petersen.score.bash
# ============================================================================

set -euo pipefail
[[ -d jobs ]] || { echo "error: run this script from the repository root"; exit 1; }

M=metrics/framework_B
FB=results/framework_b
REFS=data/sra/refs/petersen_2026

# medaka_manifest <cell>: primary records of the polished consensus -> manifest
medaka_manifest() {
  local cell="$1" label="sra_petersen_2026_v20260924_$1.medaka"
  mkdir -p "$M/$label.primary"
  printf 'sample_id\tconsensus\treference\n' > "$M/$label.manifest.tsv"
  local d acc sample
  for d in "$FB/sra_petersen_2026_v20260924_$cell"/*/; do
    sample=$(basename "$d")
    [[ -f "$d/assembly.medaka/consensus.fa" ]] || { echo "[skip] $cell/$sample: no polished consensus"; continue; }
    acc=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}
                      $c["selected"]=="True" && $c["reads"]>m{m=$c["reads"];a=$c["accession"]}
                      END{print a}' "$d/assembly/ref_selection.tsv")
    awk -v a="$acc" '/^>/{p=(substr($1,2)==a)}p' "$d/assembly.medaka/consensus.fa" \
      > "$M/$label.primary/$sample.fa"
    [[ -s "$M/$label.primary/$sample.fa" ]] || { echo "error: $acc not in $d/assembly.medaka/consensus.fa"; exit 1; }
    printf '%s\t%s\t%s\n' "$sample" "$M/$label.primary/$sample.fa" "$REFS/$sample.fasta" \
      >> "$M/$label.manifest.tsv"
  done
}
medaka_manifest refauto
medaka_manifest chain2

echo "manifests written; score with: RUN_MUTANTS=false bash jobs/framework_B/framework_B.rescore.bash"
