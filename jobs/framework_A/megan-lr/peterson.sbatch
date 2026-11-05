#!/bin/bash
# =============================================================================
# Framework A — MEGAN-LR comparator: LAST align -> MAF2DAA -> meganize -> export
# Peterson 2026 (Roihu)
# =============================================================================
# One array task per accession, same accession list
# framework_A/metamix/preprocess/peterson_taxprofiler.sbatch / framework_A/metamix/peterson.sbatch
# use. Implements the full LAST+MEGAN-LR pipeline from Huson et al. 2018 Fig. 4:
#   1. lastal (frameshift-aware, protein DB)  reads -> MAF
#   2. maf2daa                                MAF   -> DAA
#   3. daa-meganizer                          DAA   -> meganized DAA (taxonomic+functional binning)
#   4. daa2info                               meganized DAA -> taxid<TAB>count TSV
#      (this is the file ngs-benchmark's `megan_lr` adapter
#      (ngs-benchmark/ngs_benchmark/io/adapters/megan_lr.py) consumes)
#
# Reuses the SAME host-depleted reads nf-core/taxprofiler already produced for
# metaMix (framework_A/metamix/preprocess/peterson_taxprofiler.sbatch). This means
# framework_A/metamix/preprocess/peterson_taxprofiler.sbatch MUST have completed for the
# dataset before this script runs (neither framework_A/metamix/preprocess/peterson_diamond.sbatch
# nor framework_A/metamix/peterson.sbatch need to have run).
#
# Requires framework_A/megan-lr/builddb.sbatch to have completed first.
#
# USAGE: sbatch framework_A/megan-lr/peterson.sbatch
# Test one: sbatch --array=1 framework_A/megan-lr/peterson.sbatch
# =============================================================================
#SBATCH --job-name=megan_lr_peterson
#SBATCH --account=Project_2003755
#SBATCH --partition=small
#SBATCH --array=1-4                         # 4 accs (peterson_2026, no NC) — match framework_A/metamix/preprocess/peterson_taxprofiler.sbatch's accession list
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/megan_lr_peterson_%A_%a.out
#SBATCH --error=logs/megan_lr_peterson_%A_%a.err

set -euo pipefail

PROJ="/scratch/project_XXXXXXX/emtrace-benchmark"
PROJAPPL="/projappl/project_XXXXXXX"
REPO="$PROJAPPL/EMtrace-benchmark"
CONDA_ENV="$PROJAPPL/venv/megan-lr"          # bioconda: last + megan (megan-lr_setup.md item 1)

: "${databases:?set \$databases (Roihu reference-DB root)}"

export PATH="$CONDA_ENV/bin:$PATH"
# seqkit (fastq->fasta)
export PATH="$PROJAPPL/install/seqkit_dir/bin:$PATH"

RELEASE="234"                               # MUST match framework_A/megan-lr/builddb.sbatch
LASTDB_INDEX="$databases/megan_lr/refseq.abvh.1g.${RELEASE}"
MEGAN_MAPDB="$databases/megan-map-Jul2020-2-X.db"   # megan-lr_setup.md item 2

if [[ ! -s "${LASTDB_INDEX}.prj" ]]; then
  echo "[error] ${LASTDB_INDEX}.prj not found — run framework_A/megan-lr/builddb.sbatch first" >&2
  exit 1
fi
if [[ ! -f "$MEGAN_MAPDB" ]]; then
  echo "[error] $MEGAN_MAPDB not found (megan-map mapping DB, needed by daa-meganizer -mdb)" >&2
  exit 1
fi

ACCLIST="$REPO/jobs/assets/peterson_accs.txt"

ACC=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCLIST")
if [[ -z "$ACC" ]]; then
  echo "[error] no accession at line $SLURM_ARRAY_TASK_ID of $ACCLIST" >&2
  exit 1
fi

CLEAN_FASTQ="$PROJ/results/framework_a/taxprofiler/peterson_2026/samtools/fastq/${ACC}_${ACC}.unmapped_other.fastq.gz"
if [[ ! -f "$CLEAN_FASTQ" ]]; then
  echo "[error] $CLEAN_FASTQ missing — run framework_A/metamix/preprocess/peterson_taxprofiler.sbatch first" >&2
  exit 1
fi

RESDIR="$PROJ/results/framework_a/megan_lr/peterson_2026/$ACC"
mkdir -p "$RESDIR" logs

READS_FASTA="$RESDIR/${ACC}.reads.fa"
MAF="$RESDIR/${ACC}.maf"
DAA="$RESDIR/${ACC}.daa"
TAXCOUNT="$RESDIR/${ACC}_megan_taxcount.tsv"
MEGANIZED="$RESDIR/.meganized_readcount"   # new marker: forces re-meganize of runs done in default alignedBases mode

if [[ -f "$TAXCOUNT" && -f "$MEGANIZED" ]]; then
  echo "[skip] $TAXCOUNT already present"
  exit 0
fi

echo "[1/5] fastq -> fasta: $ACC"
if [[ ! -s "$READS_FASTA" ]]; then
  seqkit fq2fa "$CLEAN_FASTQ" -o "$READS_FASTA"
fi

echo "[2/5] lastal (frameshift, protein DB): $ACC"
if [[ ! -s "$MAF" ]]; then
  lastal -F15 -P "$SLURM_CPUS_PER_TASK" "$LASTDB_INDEX" "$READS_FASTA" > "$MAF.tmp"
  mv "$MAF.tmp" "$MAF"   # atomic — a time-limit kill mid-write must not leave a truncated
                          # $MAF that the `-s "$MAF"` resume check above would mistake for done
fi

echo "[3/5] maf2daa: $ACC"
if [[ ! -s "$DAA" ]]; then
  maf2daa -i "$MAF" -r "$READS_FASTA" -o "$DAA.tmp.daa" -p "$SLURM_CPUS_PER_TASK"
  mv "$DAA.tmp.daa" "$DAA"   # atomic — same resume-safety reason as the lastal step above
fi

echo "[4/5] daa-meganizer: $ACC"
if [[ ! -f "$MEGANIZED" ]]; then
  # -ram readCount: -lg defaults to alignedBases, which would make counts aligned bases, not reads
  daa-meganizer -i "$DAA" -lg -ram readCount -class -mdb "$MEGAN_MAPDB" -t "$SLURM_CPUS_PER_TASK"
  touch "$MEGANIZED"
fi

echo "[5/5] daa2info (taxid -> count export): $ACC"
{
  printf 'taxid\tcount\n'   # daa2info's own output has no header row
  daa2info -i "$DAA" -c2c Taxonomy
} > "$TAXCOUNT"

echo "[done] $ACC -> $TAXCOUNT"
