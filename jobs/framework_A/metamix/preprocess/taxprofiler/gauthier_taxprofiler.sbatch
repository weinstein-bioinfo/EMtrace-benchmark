#!/bin/bash
# =============================================================================
# Framework A — nf-core/taxprofiler preprocessing: Gauthier 2021 (Roihu)
# =============================================================================
#
# USAGE: sbatch framework_A/metamix/preprocess/gauthier_taxprofiler.sbatch
# Then: sbatch framework_A/metamix/preprocess/gauthier_diamond.sbatch (DIAMOND step, array)
#
# REQUIRES: $hostgenomes set (Roihu host-genome root) before submission
# =============================================================================
#SBATCH --job-name=taxprofiler_prep_gauthier
#SBATCH --account=Project_2003755
#SBATCH --partition=small
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G                          # placeholder — retune
#SBATCH --time=24:00:00
#SBATCH --output=logs/taxprofiler_prep_gauthier_%j.out
#SBATCH --error=logs/taxprofiler_prep_gauthier_%j.err

set -euo pipefail

PROJ="/scratch/project_XXXXXXX/emtrace-benchmark"
PROJAPPL="/projappl/project_XXXXXXX"
REPO="$PROJAPPL/EMtrace-benchmark"
: "${hostgenomes:?set \$hostgenomes (Roihu host-genome root)}"

# Manually-installed nextflow needs a working JDK — module load first
module load biojava
export PATH="$PROJAPPL/install/nextflow:$PATH"


export NXF_SINGULARITY_CACHEDIR="$PROJ/singularity_cache"
export APPTAINER_TMPDIR="${LOCAL_SCRATCH:-${TMPDIR:-/tmp}}"
export APPTAINER_CACHEDIR="${LOCAL_SCRATCH:-${TMPDIR:-/tmp}}"
mkdir -p "$NXF_SINGULARITY_CACHEDIR"

# Same GRCh38 build emtrace's own Homo_sapiens bgfilter uses
# (jobs/emtrace.config.yaml) — keeps host removal consistent across
# every pipeline in this project.
HOST_REF="$hostgenomes/GCF_000001405.40_genomic+rna.fna.gz"
if [[ ! -f "$HOST_REF" ]]; then
  echo "[error] $HOST_REF not found" >&2
  exit 1
fi

ACCLIST="$REPO/jobs/assets/gauthier_accs.txt"
SRA_DIR="$PROJ/sra"
OUTDIR="$PROJ/results/framework_a/taxprofiler/gauthier_2021"
WORKDIR="$PROJ/nf-work/taxprofiler_gauthier"
mkdir -p "$OUTDIR" "$WORKDIR" logs

SAMPLESHEET="$OUTDIR/samplesheet.csv"
DATABASES="$OUTDIR/databases.csv"   # header-only — no profiler enabled

echo "sample,run_accession,instrument_platform,fastq_1,fastq_2,fasta" > "$SAMPLESHEET"
while IFS= read -r ACC || [[ -n "$ACC" ]]; do
  [[ -z "$ACC" ]] && continue
  FASTQ="$SRA_DIR/${ACC}.fastq.gz"
  if [[ ! -f "$FASTQ" ]]; then
    echo "[error] $FASTQ not found" >&2
    exit 1
  fi
  echo "${ACC},${ACC},OXFORD_NANOPORE,${FASTQ},," >> "$SAMPLESHEET"
done < "$ACCLIST"

echo "tool,db_name,db_params,db_type,db_path" > "$DATABASES"

nextflow run nf-core/taxprofiler -r 2.0.1 \
  -profile singularity \
  -work-dir "$WORKDIR" \
  --input "$SAMPLESHEET" \
  --databases "$DATABASES" \
  --outdir "$OUTDIR" \
  --perform_longread_qc \
  --perform_longread_hostremoval \
  --hostremoval_reference "$HOST_REF" \
  --save_hostremoval_unmapped \
  --max_cpus "$SLURM_CPUS_PER_TASK" \
  --max_memory '128.GB' \
  --max_time '24.h'

echo "[done] taxprofiler host-removal: gauthier_2021 -> $OUTDIR/samtools/fastq/<ACC>_<ACC>.unmapped_other.fastq.gz per accession"
