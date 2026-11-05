#!/bin/bash
# =============================================================================
# Framework A — metaMix DIAMOND search: Buddle 2024 DNA viruses (Roihu)
# =============================================================================
#
# USAGE: sbatch framework_A/metamix/preprocess/buddle_dna_taxprofiler.sbatch  (once, whole dataset)
#        sbatch framework_A/metamix/preprocess/buddle_dna_diamond.sbatch   (this script, array)
# Requires framework_A/db_refseq_abvh_protein.{dload,build}.sbatch AND
# framework_A/metamix/preprocess/buddle_dna_taxprofiler.sbatch to have completed first.
#
# Output uses the custom 10-column outfmt (qseqid qlen sseqid slen stitle bitscore
# length pident evalue staxids) metaMix 0.3's generative.prob(blast.default=FALSE)
# expects — not the classic 12-col default. staxids comes straight from DIAMOND since
# $DIAMOND_DB was built with --taxonmap, so no separate read-length file or NCBI
# accession2taxid lookup is needed downstream. See framework_A/metamix/run.R and
# ROIHU.metamix_setup.md's "Which metaMix version" note for why blast.default=FALSE.
#
# REQUIRES: $databases set (Roihu reference-DB root) before submission
# =============================================================================
#SBATCH --job-name=metamix_diamond_buddle_dna
#SBATCH --account=Project_2003755
#SBATCH --partition=small
#SBATCH --array=1-5                         # 5 accs (buddle_2024_dna, 4 samples + 1 NC)
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/metamix_diamond_buddle_dna_%A_%a.out
#SBATCH --error=logs/metamix_diamond_buddle_dna_%A_%a.err

set -euo pipefail

PROJ="/scratch/project_XXXXXXX/emtrace-benchmark"
PROJAPPL="/projappl/project_XXXXXXX"
REPO="$PROJAPPL/EMtrace-benchmark"
: "${databases:?set \$databases (Roihu reference-DB root)}"

# Roihu's diamond module is broken — use the manually-installed binary instead.
export PATH="$PROJAPPL/install/diamond_dir:$PATH"

# Cross-pipeline-consistent DB (viral + archaea/bacteria 1/genus + human, RefSeq
# release 234) built by jobs/framework_A/db_refseq_abvh_protein.{dload,build}.sbatch
DIAMOND_DB="$databases/2026/NCBI/refseq.abvh.1g.234/diamond/diamond.refseq.abvh.1g.dmnd"
if [[ ! -f "$DIAMOND_DB" ]]; then
  echo "[error] $DIAMOND_DB not found — run jobs/framework_A/db_refseq_abvh_protein.dload.sbatch" \
       "then db_refseq_abvh_protein.build.sbatch first" >&2
  exit 1
fi

ACCLIST="$REPO/jobs/assets/buddle_dna_a_accs.txt"

# --- resolve this task's accession -----------------------------------------
ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCLIST")
if [[ -z "${ACC:-}" ]]; then
  echo "ERROR: no accession at line $SLURM_ARRAY_TASK_ID of $ACCLIST" >&2
  exit 1
fi

RESDIR="$PROJ/results/framework_a/metamix/buddle_2024_dna/$ACC"
mkdir -p "$RESDIR" logs

DIAMOND_OUT="$RESDIR/${ACC}_diamond.tab"
if [[ -f "$DIAMOND_OUT" ]]; then
  echo "[skip] $DIAMOND_OUT already present"
  exit 0
fi

# Host-depleted reads from framework_A/metamix/preprocess/buddle_dna_taxprofiler.sbatch
CLEAN_FASTQ="$PROJ/results/framework_a/taxprofiler/buddle_2024_dna/samtools/fastq/${ACC}_${ACC}.unmapped_other.fastq.gz"
if [[ ! -f "$CLEAN_FASTQ" ]]; then
  echo "[error] $CLEAN_FASTQ missing — run framework_A/metamix/preprocess/buddle_dna_taxprofiler.sbatch first" >&2
  exit 1
fi

echo "diamond blastx (custom 10-col outfmt for metaMix 0.3's blast.default=FALSE): $ACC"
diamond blastx \
  -d "$DIAMOND_DB" \
  -q "$CLEAN_FASTQ" \
  -o "$DIAMOND_OUT.tmp" \
  -f 6 qseqid qlen sseqid slen stitle bitscore length pident evalue staxids \
  --max-target-seqs 10 \
  --evalue 1 \
  --sensitive \
  --threads "$SLURM_CPUS_PER_TASK"

# Drop hits with no staxids (col 10 empty) — subject accessions the --taxonmap DB build
# doesn't have a taxid for. metaMix 0.3's generative.prob(blast.default=FALSE) only
# filters NA taxids when its taxonID column happens to parse as character (i.e. only if
# some OTHER row in the file has a semicolon-separated multi-taxid hit); on a file with
# zero multi-taxid hits, fread() reads the column as pure integer and that NA-filter
# branch never runs, so a single empty-staxid row reaches sparseMatrix() unfiltered and
# crashes the whole job with "'i' and 'j' must not contain NA" (confirmed 2026-09-23,
# ERR13488553_diamond.tab, 1/288 rows). Filtering here rather than in metaMix's own code.
awk -F'\t' -v OFS='\t' '$10 != ""' "$DIAMOND_OUT.tmp" > "$DIAMOND_OUT"
rm -f "$DIAMOND_OUT.tmp"

echo "[done] $ACC -> $DIAMOND_OUT"

