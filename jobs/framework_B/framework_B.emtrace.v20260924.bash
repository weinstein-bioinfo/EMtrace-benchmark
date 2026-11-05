#!/usr/bin/env bash
# ============================================================================
# Framework B — EMtrace full re-run on the mutant matrix (v20260924)
# ============================================================================
#
# Re-runs the WHOLE pipeline (`--pipe main,assembly`) on the pure and mixed
# mutant sets with the EMtrace code as of 2026-09-24, and with the
# reference chosen by the pipeline itself (--ref_auto) instead of handed to it
# (--refaccs).  Changes since the original Framework B runs (May/June 2026):
#
#   annotation  all hits kept per read (hit-mode change, 17 Sep)
#   assembly    reference: automatic selection (ref_auto, 5848bb6) —
#                 primary: reads mapped to assembly.database (RefSeq viral),
#                 secondary: annotation hits; see assembly/ref_selection.tsv
#               reads: all host-filtered reads (assembly.input_reads)
#               mapping: minimap2 -k 11 -w 5, default -m/-s (assembly.minimap)
#               masking: samtools depth -J (deleted positions count as covered)
#               calling: bcftools mpileup -X ont-sup (ONT indels), AF tiers,
#                        norm -m -any (463a0ff, 5b5266d)
#               min_coverage: MIN_COVERAGE below (10; config default is 20)
#
#   CELL            which variant to run; sets --ann1 and the ref_auto source,
#                   and names the output batch {pure,mixed}_v20260924_<CELL>:
#
#     refauto   ORIGINAL.  Annotation index minimap.refseq.abv.1g (k15/w10),
#               references from annotation hits only (use_assembly_db false).
#               One database governs detection, selection and the reference
#               sequence -- the cleanest setup, and it selected correctly for
#               all 14 COVID and all 14 DENV samples.  But a virus annotation
#               misses cannot be selected: LASV_L 10-30 % and LASV_S 25-30 %
#               ended with no reference and no consensus.
#
#     k11w5     SAME single-database design, sensitive annotation index
#               minimap.refseq.abv.1g.k11w5 (-k 11 -w 5, prebuilt .mmi).
#               Fixes the detection step itself, so selection needs no second
#               database.  Probed 25 Sep 2026 against the Lassa RefSeq pair:
#               LASV_L 25/30 % recruits 32/19 reads where the default index
#               recruits 0, and LASV_S 30 % recruits 70 against 22 -- all past
#               ref_auto.min_reads=10.  Runs the two Lassa viruses only (the
#               RUN_* switches below); COVID and DENV already select correctly
#               on the default index, so they would only cost runtime.
#
#     chain2    Two-stage annotation, then a nucleotide reference.  Round 1 is a
#               CHAIN: minimap.refseq.abv.1g (nucleotide, archaea+bacteria+
#               viruses) runs on all reads, then diamond.uniref.vi (DIAMOND
#               blastx, UniRef100 viral proteins: 2,787,313 seqs / 1.08 G aa)
#               runs on whatever round 1 left UNANNOTATED -- the chain rule in
#               steps/annotate.py, "subsequent searches target unannotated
#               reads".  Protein search survives divergence that nucleotide
#               seeding cannot: probed 25 Sep 2026 on LASV_L_20 (261 reads),
#               DIAMOND recovers 142 Lassa reads at bitscore >= 70 and 119 at
#               >= 80, with ZERO non-Lassa hits at either, against 9 reads from
#               nucleotide annotation alone.  4 s for all 261 reads.
#
#               Reference selection is then split by design, because a protein
#               hit cannot name a genome:
#                 identification -> the chain (DIAMOND taxids; the .dmnd carries
#                                  --taxonmap, staxids are populated)
#                 reference      -> ref_auto's PRIMARY search, i.e. the viral
#                                  read partition mapped onto assembly.database
#                                  (minimap.refseq.vi, 18,795 nucleotide seqs)
#                                  at -k 11 -w 5.
#               Hence --ref_auto_assembly_db together with
#               --no-ref_auto_annotation: annotation is deliberately barred from
#               supplying accessions, since its hits here are UniRef100_* protein
#               IDs that _fetch_accessions() cannot resolve to a genome (it would
#               either fail or hand assembly a protein sequence).
#
#               This also turns the caveat in the refprim cell into the
#               mechanism: the chain POPULATES reads.ann1.vi.fq, so
#               _find_viral_reads() (assembly.py:153) hands the primary search a
#               real viral partition instead of one stray read.
#
#               Run via the named strategy `--anns vi.chain2`, redefined
#               25 Sep 2026 (in BOTH emtrace.config.yaml and the package
#               default.yaml) from minimap.refseq.abv to minimap.refseq.abv.1g.
#               The old definition could not run: only the .acc2taxid and .stats
#               sidecars of minimap.refseq.abv are on this machine, no FASTA and
#               no .mmi, while `databases --list` still reports it installed.
#               refseq.abv is also the full 28.2 Gbp / 30,639-sequence set
#               against 9.1 Gbp for .1g, so keeping .abv would have changed the
#               nucleotide database and the annotation method at once.  The
#               strategy was unreferenced anywhere else, so redefining it was
#               safe.  NOTE `--anns` overrides `--ann1`/`--ann2`.
#
#     refprim   Annotation index unchanged; ref_auto ALSO searches
#               assembly.database (RefSeq viral, 18,795 seqs) at -k 11 -w 5
#               from assembly.ref_auto.minimap.  Requires use_assembly_db: true
#               and minimap k/w 11/5 in jobs/emtrace.config.yaml (the
#               CLI flag toggles only the source; k/w are config-only).
#               NOTE: for background-embedded (mixed) samples the selection
#               search is handed the annotation viral partition whenever that
#               file is non-empty (steps/assembly.py:153), so a single stray
#               read there suppresses this cell's effect.  Pure samples, whose
#               partition is empty, are unaffected.
#
# Petersen 2026 SRA (RUN_PETERSEN_REFAUTO / RUN_PETERSEN_CHAIN2): the same
# v20260924 pipeline on the 4 real Petersen libraries (data/sra/<ACC>.fastq[.gz]),
# once per cell, refauto and chain2, regardless of CELL.  The earlier
# sra_petersen_2026 runs were given their reference (--refaccs); here ref_auto
# picks it, as for the mutants.
#
# Gauthier 2021 SRA (RUN_GAUTHIER_REFAUTO / RUN_GAUTHIER_CHAIN2): the same run on
# the 21 Gauthier SARS-CoV-2 libraries (16 patients, Ct 14.7-30.6, + 5 negative
# controls), replacing the June 2026 sra_gauthier_2021 run (--refaccs
# OX637002.1, min_coverage 5, pre-5b5266d calling); the project notes.
# Off by default: 21 libraries, several hours.
#
# Nothing existing is overwritten: results go to NEW batch directories,
#   results/framework_b/pure_v20260924_<CELL>/<SAMPLE>/
#   results/framework_b/mixed_v20260924_<CELL>/<SAMPLE>/
#   results/framework_b/sra_petersen_2026_v20260924_{refauto,chain2}/<ACC>/
#   results/framework_b/sra_gauthier_2021_v20260924_{refauto,chain2}/<ACC>/
# e.g. pure_v20260924_k11w5/ for the default cell, with the same <SAMPLE> names
# as pure/ and mixed/, so the metrics scripts score them unchanged (the three
# cells are registered in PIPELINE_DIRS in framework_B.metrics.bash).
# The script refuses to write into pure/ or mixed/.
#
# Scoring note: ref_auto may choose another accession of the same virus than
# the one the mutants were made from (e.g. NC_045512.2 for COVID, RefSeq
# Lassa NC_004296/7.1 instead of MG812630/1: the annotation and assembly DBs
# are RefSeq-only).  ngs-benchmark aligns each consensus to the mutant truth
# with global MAFFT, so a different accession of the same virus still scores.
#
# Skip logic:
#   A sample is skipped once <SAMPLE>/.done.v20260924 exists (written after a
#   successful run).  Not keyed on consensus.fa: with ref_auto a sample may
#   legitimately end without one (no reference selected).  To redo a sample:
#     rm -rf results/framework_b/pure_v20260924_<CELL>/<SAMPLE>
#
# Timing:   each run is wrapped by bench_time.py, appending one row
#           sample,user,system,cpu,total,max_rss_mb,exit_status to
#           metrics/framework_B/<batch>.timings.csv.
#
# Log layout:
#   logs/emtrace/<batch>/<SAMPLE>/   pipeline log dir (--logs)
#   logs/<batch>/<SAMPLE>.log           bash stdout/stderr (tee)
#
# Run from the repository root:
#   bash jobs/framework_B/framework_B.emtrace.v20260924.bash
# ============================================================================

set -o pipefail

[[ -d jobs ]] || { echo "error: run this script from the repository root"; exit 1; }

VERSION="v20260924"
CELL="${CELL:-chain2}"   # chain2 | k11w5 | refauto | refprim  (see header)
TAG="${VERSION}_${CELL}"
RESULTS_ROOT="results/framework_b"
LOG_DIR="logs"
PIPE_LOG_DIR="logs/emtrace"
METRICS_DIR="metrics/framework_B"
CONFIG="jobs/emtrace.config.yaml"
THREADS=8
MIN_DIAMOND_BITS=80      # chain2 only; config default is 200, which keeps 11 of 196 hits
# ANN_FLAGS, REF_DB_FLAG and CELL_FLAGS are set from CELL below.

# =============================================================================
# Control switches — set to true/false before running
# =============================================================================
RUN_PURE=false            # pure/  mutants  → pure_${TAG}/   (4 viruses × 7 deviations)
RUN_MIXED=false           # mixed/ mutants  → mixed_${TAG}/  (same, + 500 000 background reads)

# Petersen 2026 SRA samples (4 real libraries: ZIKV, WNV, SARS-CoV-2, DENV),
# one switch per cell, independent of CELL above → sra_petersen_2026_v20260924_<cell>/
RUN_PETERSEN_REFAUTO=false   # default annotation, references from annotation hits
RUN_PETERSEN_CHAIN2=false    # nucleotide + DIAMOND chain, references from the RefSeq search
PETERSEN_ACCS=(SRR31174237 ERR14941561 SRR15356294 SRR32419728)

# Gauthier 2021 SARS-CoV-2 libraries (BioProject PRJNA752146), one switch per
# cell → sra_gauthier_2021_v20260924_<cell>/. Same 21 accessions as
# framework_B.sra.emtrace.bash: 16 patients, then 5 negative controls.
RUN_GAUTHIER_REFAUTO=true  # default annotation, references from annotation hits
RUN_GAUTHIER_CHAIN2=false   # nucleotide + DIAMOND chain, references from the RefSeq search
GAUTHIER_ACCS=(
  SRR15356303 SRR15356328 SRR15356305 SRR15356306 SRR15356327 SRR15356326
  SRR15356325 SRR15356329 SRR15356307 SRR15356298 SRR15356308 SRR15356322
  SRR15356310 SRR15356311 SRR15356321 SRR15356320                           # patients
  SRR15356319 SRR15356312 SRR15356295 SRR15356318 SRR15356317               # negative controls
)
SRA_DIR="data/sra"

# Set for CELL=k11w5, whose question is the two divergent Lassa segments.
# Turn COVID and DENV back on when running the refauto or refprim cells.
RUN_COVID=true
RUN_DENV=true
RUN_LASV_L=true
RUN_LASV_S=true

# set_cell <cell>: sets ANN_FLAGS, REF_DB_FLAG, CELL_FLAGS and STAMP for one cell
set_cell() {
  CELL_FLAGS=()
  case "$1" in
    refauto) ANN_FLAGS=(--ann1 minimap.refseq.abv.1g)       ; REF_DB_FLAG="--no-ref_auto_assembly_db" ;;
    k11w5)   ANN_FLAGS=(--ann1 minimap.refseq.abv.1g.k11w5) ; REF_DB_FLAG="--no-ref_auto_assembly_db" ;;
    refprim) ANN_FLAGS=(--ann1 minimap.refseq.abv.1g)       ; REF_DB_FLAG="--ref_auto_assembly_db" ;;
    chain2)  ANN_FLAGS=(--anns vi.chain2)
             REF_DB_FLAG="--ref_auto_assembly_db"
             CELL_FLAGS=(--no-ref_auto_annotation --min_diamond_bits "$MIN_DIAMOND_BITS") ;;
    *) echo "error: CELL must be chain2, k11w5, refauto or refprim"; exit 1 ;;
  esac
  STAMP=".done.${VERSION}.$1"
}
set_cell "$CELL"

DEVS=(00 05 10 15 20 25 30)
MIN_COVERAGE=10          # assembly.min_coverage (--min_coverage)
DRY_RUN=false            # list what would run; touch nothing
# =============================================================================

VIRUSES=()
[[ "$RUN_COVID"  == true ]] && VIRUSES+=(COVID)
[[ "$RUN_DENV"   == true ]] && VIRUSES+=(DENV)
[[ "$RUN_LASV_L" == true ]] && VIRUSES+=(LASV_L)
[[ "$RUN_LASV_S" == true ]] && VIRUSES+=(LASV_S)

RAN=0; SKIPPED=0; FAILED=0; NOCONS=0
FAILED_SAMPLES=(); NOCONS_SAMPLES=()

# ---------------------------------------------------------------------------
# run_sample <batch> <sample> <fastq>
# ---------------------------------------------------------------------------
run_sample() {
  local batch="$1" sample="$2" fastq="$3"
  local rdir="$RESULTS_ROOT/$batch/$sample"

  case "$batch" in pure|mixed) echo "error: refusing to write into $batch/"; exit 1;; esac

  if [[ -f "$rdir/$STAMP" ]]; then
    echo "[skip] $batch/$sample — done ($STAMP)"; SKIPPED=$((SKIPPED + 1)); return 0
  fi
  if [[ ! -f "$fastq" ]]; then
    echo "[skip] $batch/$sample — input FASTQ not found: $fastq"; SKIPPED=$((SKIPPED + 1)); return 0
  fi
  if [[ "$DRY_RUN" == true ]]; then
    echo "[dry-run] $batch/$sample  -1 $fastq  --pipe main,assembly --ref_auto $REF_DB_FLAG ${CELL_FLAGS[*]} --min_coverage $MIN_COVERAGE"
    return 0
  fi

  mkdir -p "$LOG_DIR/$batch" "$PIPE_LOG_DIR/$batch"
  echo "[run]  $batch/$sample"
  local rc=0
  python3 jobs/bench_time.py --csv "$METRICS_DIR/${batch}.timings.csv" --label "$sample" -- \
    emtrace run \
    -1 "$fastq" \
    -s "$sample" \
    -r "$RESULTS_ROOT/$batch" \
    -t "$THREADS" \
    --ref_auto \
    "$REF_DB_FLAG" \
    ${CELL_FLAGS[@]+"${CELL_FLAGS[@]}"} \
    --min_coverage "$MIN_COVERAGE" \
    --bgfilter Homo_sapiens \
    "${ANN_FLAGS[@]}" \
    --config "$CONFIG" \
    --pipe main,assembly \
    --logs "$PIPE_LOG_DIR/$batch/$sample" \
    -v \
    2>&1 | tee "$LOG_DIR/$batch/$sample.log" || rc=$?

  if [[ $rc -ne 0 ]]; then
    echo "[FAIL] $batch/$sample — exit $rc (see $LOG_DIR/$batch/$sample.log)"
    FAILED=$((FAILED + 1)); FAILED_SAMPLES+=("$batch/$sample"); return 0
  fi
  date -u +"%Y-%m-%dT%H:%M:%SZ" > "$rdir/$STAMP"
  RAN=$((RAN + 1))
  if [[ ! -f "$rdir/assembly/consensus.fa" ]]; then
    echo "[note] $batch/$sample — no consensus.fa (no reference selected? see assembly/ref_selection.tsv)"
    NOCONS=$((NOCONS + 1)); NOCONS_SAMPLES+=("$batch/$sample")
  fi
}

mkdir -p "$METRICS_DIR"
echo "Framework B re-run $VERSION [$CELL]: ref_auto $REF_DB_FLAG ${CELL_FLAGS[*]}, min_coverage $MIN_COVERAGE"
echo "annotation: ${ANN_FLAGS[*]}"
echo "EMtrace: $(emtrace --version 2>/dev/null || echo '?')"

if [[ "$RUN_PURE" == true ]]; then
  echo "=== pure_${TAG} ==="
  for v in "${VIRUSES[@]}"; do for dev in "${DEVS[@]}"; do
    s="${v}_${dev}"
    run_sample "pure_${TAG}" "$s" "data/petersen_2026/fastq/${s}/${s}.fastq.gz"
  done; done
fi

if [[ "$RUN_MIXED" == true ]]; then
  echo "=== mixed_${TAG} ==="
  for v in "${VIRUSES[@]}"; do for dev in "${DEVS[@]}"; do
    s="${v}_${dev}_mixed"
    run_sample "mixed_${TAG}" "$s" "data/petersen_2026/mixed/${s}.fastq.gz"
  done; done
fi

# SRA sets (Petersen 2026, Gauthier 2021): real libraries, so no mutant truth; the
# reference is chosen by ref_auto (the older sra_* runs were given it with --refaccs).
# run_sra <set> <cell> <acc>...   batch = sra_<set>_v20260924_<cell>
run_sra() {
  local set="$1" cell="$2"; shift 2
  local batch="sra_${set}_${VERSION}_${cell}" acc fq
  set_cell "$cell"
  echo "=== ${batch}  (annotation: ${ANN_FLAGS[*]}; ref_auto $REF_DB_FLAG ${CELL_FLAGS[*]}) ==="
  for acc in "$@"; do
    fq="$SRA_DIR/${acc}.fastq"; [[ -f "$fq" ]] || fq="${fq}.gz"
    run_sample "$batch" "$acc" "$fq"
  done
}
[[ "$RUN_PETERSEN_REFAUTO" == true ]] && run_sra petersen_2026 refauto "${PETERSEN_ACCS[@]}"
[[ "$RUN_PETERSEN_CHAIN2"  == true ]] && run_sra petersen_2026 chain2  "${PETERSEN_ACCS[@]}"
[[ "$RUN_GAUTHIER_REFAUTO" == true ]] && run_sra gauthier_2021 refauto "${GAUTHIER_ACCS[@]}"
[[ "$RUN_GAUTHIER_CHAIN2"  == true ]] && run_sra gauthier_2021 chain2  "${GAUTHIER_ACCS[@]}"

echo
echo "=== summary ==="
echo "ran: $RAN   skipped: $SKIPPED   failed: $FAILED   without consensus: $NOCONS"
[[ $FAILED -gt 0 ]] && printf '  failed: %s\n' "${FAILED_SAMPLES[@]}"
[[ $NOCONS -gt 0 ]] && printf '  no consensus: %s\n' "${NOCONS_SAMPLES[@]}"
echo "timings: $METRICS_DIR/<batch>.timings.csv"
echo
echo "Next: score with framework_B.metrics.bash — batches pure_${TAG} / mixed_${TAG}"
echo "are registered in PIPELINE_DIRS (leave MIN_COVERAGE unset: assembly/ holds"
echo "the d${MIN_COVERAGE} run)."
