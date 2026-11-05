#!/usr/bin/env bash
# ============================================================================
# Framework C — reads aligned in the final reference-based assembly (Figure b legend)
# ============================================================================
# framework_C.coverage_over_time.bash keeps only consensus.fa and the low-coverage
# BED of each point, not its BAM. This script repeats the final point's alignment
# with the assembly's own settings (assembly.minimap in emtrace.config.yaml:
# minimap2 -ax map-ont -k 11 -w 5, no secondary alignments) on the host-filtered
# reads of all batches, and counts primary mapped reads per reference segment
# (samtools view -F 0x904). The count is an upper bound on the reads the consensus
# used: markdup is not repeated (for single long reads it removes few).
#
# Output: metrics/framework_C/coverage_mapped_reads.tsv
#   library, segment, mapped_reads, library_reads (Module C's final reads_cum)
#
# Run from the repository root:
#   bash jobs/framework_C/framework_C.coverage_mapped_reads.bash
# ============================================================================

set -euo pipefail
[[ -d jobs ]] || { echo "error: run this script from the repository root"; exit 1; }

FC="results/framework_C"
OUT="metrics/framework_C/coverage_mapped_reads.tsv"
THREADS=8
# library|reference, as in framework_C.coverage_over_time.bash
LIBS=(
  "Orbivirus_kevo.refseq_abv1g_k11w5|data/orbivirus.kevo/refs/kevoov_PQ165789-PQ165798.fasta"
  "TBEV_FBC92293.refseq_abv1g|data/tbev.inhouse/refs/tbev_NC_001672.fasta"
  "ticks_35-33.refseq_abv1g_k11w5|data/ticks.ergunay/refs/cchfv_aigai_CAP12_MN811031-MN811037-MN811034.fasta"
  "ticks_58-37.refseq_abv1g_k11w5|data/ticks.ergunay/refs/ticks_58-37_cchfv_CAP12+jmtv_tulcea1.fasta"
)

printf "library\tsegment\tmapped_reads\tlibrary_reads\n" > "$OUT.tmp"
for entry in "${LIBS[@]}"; do
  lib=${entry%%|*} ref=${entry#*|}
  total=$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "reads_cum") c = i; next }
                      { last = $c } END { print last }' "$FC/$lib/timeseries.tsv")
  echo "[map] $lib ($total reads) -> $(basename "$ref")"
  cat "$FC/$lib/pipeline/$lib"/batches/batch_*/reads/reads.trim.hflt.fq.gz |
    minimap2 -ax map-ont -k 11 -w 5 --secondary=no -t "$THREADS" "$ref" - 2> /dev/null |
    samtools view -F 0x904 - | cut -f3 | sort | uniq -c |
    awk -v lib="$lib" -v total="$total" 'BEGIN { OFS = "\t" } { print lib, $2, $1, total }' >> "$OUT.tmp"
done
mv "$OUT.tmp" "$OUT"
echo "wrote $OUT"
cat "$OUT"
