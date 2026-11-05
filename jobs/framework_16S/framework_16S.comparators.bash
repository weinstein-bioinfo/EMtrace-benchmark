#!/usr/bin/env bash
# ============================================================================
# Framework 16S — comparators: Emu and Spaghetti, every library  (local/Mac)
# ============================================================================
#
# The one comparator runner for the 16S benchmark. EMtrace cells are
# jobs/framework_16S/framework_16S.emtrace.bash (Roihu); these two run here, because
# neither is worth a Slurm allocation (45 min and 18 min on barcode79).
#
# Panel B.1: each tool on its OWN database, at its published defaults. Those
# defaults are the operating point the tool was evaluated at, so overriding one
# turns the comparator into a different tool. See roadmap section 7.1.
#   Emu        its own 49 k-sequence DB, -N 50, 10-read threshold
#   Spaghetti  SILVA, primary alignments only, GENUS RANK ONLY by design
#
# Panel B.2 for Spaghetti (--spag-db emu): the same pipeline against Emu's
# sequences. Its database is only dna-sequences.mmi plus a
# "Feature ID<TAB>Taxon" table, and taxonomyTable.py reads that generically,
# so nothing in it is SILVA-specific -- the feasibility concern in roadmap
# section 7.2 does not hold. Built here by build_spag_emu_db().
#
# Spaghetti's yacrd chimera removal stays SKIPPED: the all-against-all step
# grows with the square of the read count (~2.5 h and ~490 GB of alignments at
# 88 k reads). Both libraries skip it, so the setup is at least consistent.
#
# Usage (from the repo root):
#   bash jobs/framework_16S/framework_16S.comparators.bash                  # both tools, both libraries
#   bash jobs/framework_16S/framework_16S.comparators.bash --tool emu
#   bash jobs/framework_16S/framework_16S.comparators.bash --library synthgut
#   bash jobs/framework_16S/framework_16S.comparators.bash \
#       --tool spaghetti --spag-db emu          # Panel B.2 only
#   bash jobs/framework_16S/framework_16S.comparators.bash \
#       --library zymo_winand --reads data/sra/16S/SRR10391201.fastq.gz
#
# Outputs (roadmap section 9):
#   results/framework_16S/<library>/emu/<library>_rel-abundance.tsv
#   results/framework_16S/<library>/spaghetti.tsv        (scorer input, SILVA)
#   results/framework_16S/<library>/spaghetti.emu.tsv    (scorer input, Emu ref)
#   results/framework_16S/<library>/spaghetti[.emu]/     (full run trees)
#   each with a *.provenance.txt recording tool version, DB build and thresholds
#
# Helpers kept separate on purpose:
#   framework_16S.spaghetti.sh        vendored upstream runner, not ours to edit
#   framework_16S.spaghetti2table.py  its output -> the scorer's `table` format
#
# Prints the manifest rows to paste into jobs/framework_16S/framework_16S.manifest.tsv.
# ============================================================================

set -euo pipefail

# Reads per library. Adding a library is one line here.
# Reads per library on this Mac. Several candidates per library, because
# barcode79 is .fq.gz while synthgut is a plain .fastq -- hard-coding one
# suffix silently skips the other library.
reads_for() {
  local lib="$1" c
  local -a cands
  case "$lib" in
    barcode79) cands=("data/mock_community.16S/barcode79.fq.gz") ;;
    synthgut)  cands=("data/sra/16S/SRR14326315.fastq"
                      "data/sra/16S/SRR14326315.fastq.gz") ;;
    *)         cands=("data/sra/16S/$lib.fastq"
                      "data/sra/16S/$lib.fastq.gz"
                      "data/sra/16S/$lib.fq.gz") ;;
  esac
  for c in "${cands[@]}"; do
    [[ -f "$c" ]] && { echo "$c"; return 0; }
  done
  echo ""
}

SAMPLES=""          # empty = every library in reads_for()
READS_OVERRIDE=""
TOOL="both"
# Which reference(s) to run Spaghetti against:
#   silva  its own, as published        -> Panel B.1
#   emu    Emu's sequences, unchanged   -> Panel B.2
SPAG_DBS="${SPAG_DBS:-silva emu}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --sample|--library) SAMPLES="$2"; shift 2 ;;
    --reads)  READS_OVERRIDE="$2"; shift 2 ;;
    --tool)   TOOL="$2";   shift 2 ;;
    --spag-db) SPAG_DBS="$2"; shift 2 ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
case "$TOOL" in emu|spaghetti|both) ;; *) echo "--tool must be emu|spaghetti|both" >&2; exit 2 ;; esac
for d in $SPAG_DBS; do
  case "$d" in silva|emu) ;; *) echo "--spag-db must be silva and/or emu" >&2; exit 2 ;; esac
done
[[ -n "$SAMPLES" ]] || SAMPLES="barcode79 synthgut"

# Build a Spaghetti-format database from Emu's sequences. Spaghetti needs
# exactly two files -- dna-sequences.mmi and a "Feature ID<TAB>Taxon" table --
# and taxonomyTable.py reads the second generically, so nothing here is
# SILVA-specific. Emu names every sequence "<taxid>:emu_db:<n>" and ships a
# taxonomy.tsv keyed by that taxid, so the QIIME-style lineage is a join with
# no external lookup. Probed on barcode79.1ksample: 49,243/49,243 sequences
# resolve, 916/1000 reads map, all eight mock genera recovered.
build_spag_emu_db() {
  local db="$1"
  if [[ -f "$db/dna-sequences.mmi" && -s "$db/taxonomy.tsv" ]]; then
    echo "[skip] Spaghetti emu DB exists ($db)"; return 0
  fi
  local fa="data/comparators/emu_db/species_taxid.fasta"
  local tx="data/comparators/emu_db/taxonomy.tsv"
  [[ -f "$fa" && -f "$tx" ]] || { echo "Missing Emu DB ($fa)" >&2; return 1; }
  mkdir -p "$db"
  echo "[build] Spaghetti emu DB -> $db"
  python3 - "$fa" "$tx" "$db/taxonomy.tsv" <<'TAXPY'
import csv, sys
fa, tx, out = sys.argv[1:4]
tax = {r["tax_id"]: r for r in csv.DictReader(open(tx), delimiter="\t")}
RANKS = [("superkingdom","d"),("phylum","p"),("class","c"),("order","o"),
         ("family","f"),("genus","g"),("species","s")]
n = miss = 0
with open(out, "w") as fh:
    fh.write("Feature ID\tTaxon\n")
    for line in open(fa):
        if not line.startswith(">"):
            continue
        sid = line[1:].split()[0]
        t = tax.get(sid.split(":")[0])
        if t is None:
            miss += 1
            continue
        fh.write(sid + "\t" + "; ".join(
            f"{p}__{(t.get(k) or '').strip().replace(' ', '_')}"
            for k, p in RANKS) + "\n")
        n += 1
print(f"  taxonomy.tsv: {n} features, {miss} unresolved")
if miss:
    raise SystemExit("unresolved sequence taxids — refusing a partial database")
TAXPY
  [[ $? -eq 0 ]] || return 1
  minimap2 -x map-ont -d "$db/dna-sequences.mmi" "$fa" 2>/dev/null \
    || { echo "minimap2 index failed" >&2; return 1; }
  cp "$fa" "$db/dna-sequences.fasta"
}

REPO="$PWD"
LOG_DIR="logs/framework_16S"
TIMINGS="metrics/framework_16S/timings/framework_16S.comparators.timings.csv"
THREADS="${THREADS:-8}"
mkdir -p "$LOG_DIR" "$(dirname "$TIMINGS")"

# --- one library at a time ---------------------------------------------------
for SAMPLE in $SAMPLES; do

READS="${READS_OVERRIDE:-$(reads_for "$SAMPLE")}"
if [[ ! -f "$READS" ]]; then
  echo "[WARN] $SAMPLE: reads not found (${READS:-<unknown library>}) — skipping."
  continue
fi

# One tree keyed library -> tool (roadmap section 9).
OUTROOT="results/framework_16S/$SAMPLE"
mkdir -p "$OUTROOT"

echo
echo "=== $SAMPLE ==="
echo "reads:   $READS"
echo "tools:   $TOOL   (threads $THREADS)"

# ---------------------------------------------------------------------------
# Emu
# ---------------------------------------------------------------------------
if [[ "$TOOL" == "emu" || "$TOOL" == "both" ]]; then
  EMU_PY="${EMU_PY:-data/comparators/emu-venv/bin/python}"
  EMU_SCRIPT="${EMU_SCRIPT:-data/comparators/emu-src/emu}"
  EMU_DB="${EMU_DB:-data/comparators/emu_db}"
  EMU_OUT="$OUTROOT/emu"

  [[ -x "$EMU_PY" ]] || { echo "Missing Emu venv python: $EMU_PY (set EMU_PY)" >&2; exit 1; }
  [[ -f "$EMU_SCRIPT" ]] || { echo "Missing Emu script: $EMU_SCRIPT (set EMU_SCRIPT)" >&2; exit 1; }
  [[ -f "$EMU_DB/species_taxid.fasta" && -f "$EMU_DB/taxonomy.tsv" ]] || {
    echo "Emu database incomplete: $EMU_DB (want species_taxid.fasta + taxonomy.tsv)" >&2; exit 1; }
  command -v minimap2 > /dev/null || { echo "minimap2 not on PATH — Emu shells out to it" >&2; exit 1; }

  mkdir -p "$EMU_OUT"

  # Emu holds a read x candidate-taxon likelihood matrix in memory (up to 50
  # candidates per read) and peaked at 2.6 GB on barcode79's 100 k reads. This
  # library is 5.2x that, so expect ~13 GB against the Mac's 16 GB: roadmap
  # "Emu's memory may not fit on library 2". Read the max_rss_mb column in the
  # timings CSV afterwards — it decides whether library 2 can run here at all.
  echo "Emu: ~4 h expected; watch memory (~13 GB projected on 16 GB)."

  python3 jobs/bench_time.py \
    --csv "$TIMINGS" --label "emu.$SAMPLE" -- \
    "$EMU_PY" "$EMU_SCRIPT" abundance \
      --db "$EMU_DB" \
      --type map-ont \
      --output-basename "$SAMPLE" \
      --threads "$THREADS" \
      --keep-counts \
      --output-dir "$EMU_OUT" \
      "$READS" \
    2>&1 | tee "$LOG_DIR/emu.$SAMPLE.log"

  {
    echo "date:        $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "emu:         $("$EMU_PY" "$EMU_SCRIPT" --version 2>&1 | head -1)"
    echo "emu_script:  $EMU_SCRIPT"
    echo "python:      $("$EMU_PY" --version 2>&1)"
    echo "minimap2:    $(minimap2 --version 2>&1 | head -1)"
    echo "database:    $EMU_DB"
    echo "db_files:"
    ls -la "$EMU_DB" | sed 's/^/  /'
    echo "sample:      $SAMPLE"
    echo "reads:       $READS"
    echo "reads_md5:   $(md5 -q "$READS" 2>/dev/null || md5sum "$READS" | cut -d' ' -f1)"
    echo "threads:     $THREADS"
    echo "thresholds:  Emu defaults (10-read minimum; no {genus} sp. exclusion)"
  } > "$EMU_OUT/emu.provenance.txt"

  echo "Emu done -> $EMU_OUT/${SAMPLE}_rel-abundance.tsv"
fi

# ---------------------------------------------------------------------------
# Spaghetti
# ---------------------------------------------------------------------------
if [[ "$TOOL" == "spaghetti" || "$TOOL" == "both" ]]; then
 for SPAG_WHICH in $SPAG_DBS; do
  # silva keeps the established paths so Part 1 results stay addressable;
  # emu is suffixed.
  if [[ "$SPAG_WHICH" == silva ]]; then
    SPAG_RUN="$REPO/$OUTROOT/spaghetti"
    SPAG_OUT="$OUTROOT/spaghetti.tsv"
    export SPAG_DB="${SPAG_DB_SILVA:-$REPO/data/comparators/spaghetti_db}"
  else
    SPAG_RUN="$REPO/$OUTROOT/spaghetti.emu"
    SPAG_OUT="$OUTROOT/spaghetti.emu.tsv"
    export SPAG_DB="$REPO/data/comparators/spaghetti_db.emu"
    build_spag_emu_db "data/comparators/spaghetti_db.emu" || { echo "[WARN] emu DB build failed — skipping"; continue; }
  fi
  echo "-- Spaghetti on $SPAG_WHICH ($SPAG_DB)"
  export SPAG_THREADS="$THREADS"
  # Part 1 skipped chimera removal and Part 2 keeps the setup unchanged: the
  # all-vs-all step grows with the square of the read count and was already
  # ~2.5 h at 88 k reads (the project notes).
  export SKIP_YACRD="${SKIP_YACRD:-true}"

  [[ -f "$SPAG_DB/dna-sequences.mmi" && -f "$SPAG_DB/taxonomy.tsv" ]] || {
    echo "Spaghetti DB incomplete: $SPAG_DB" >&2; exit 1; }

  export PATH="$HOME/miniconda3/envs/yacrd/bin:\
$REPO/data/comparators/spaghetti-venv/bin:\
$REPO/data/comparators/spaghetti-src/module1:$PATH"
  for t in porechop NanoFilt NanoStat minimap2 filterPAF.py merfePAF.py taxonomyTable.py; do
    p="$(command -v "$t" || true)"
    printf "  %-18s %s\n" "$t" "${p:-MISSING}"
    [[ -n "$p" ]] || { echo "Spaghetti needs $t on PATH" >&2; exit 1; }
  done
  [[ "$SKIP_YACRD" == true ]] || command -v yacrd > /dev/null || {
    echo "SKIP_YACRD=false but yacrd is not on PATH" >&2; exit 1; }

  # The runner globs "$1*fastq" for PLAIN uncompressed input and writes every
  # intermediate beside it, so the library gets a directory of its own. The
  # decompressed reads are ~1.6 GB here and porechop/NanoFilt copy them twice
  # more: budget ~5 GB.
  mkdir -p "$SPAG_RUN"
  if [[ ! -s "$SPAG_RUN/$SAMPLE.fastq" ]]; then
    # synthgut is already plain .fastq; gzip -dc on it leaves an empty file.
    if [[ "$READS" == *.gz ]]; then
      echo "Decompressing reads into $SPAG_RUN (~1.6 GB)"
      gzip -dc "$READS" > "$SPAG_RUN/$SAMPLE.fastq"
    else
      echo "Copying reads into $SPAG_RUN"
      cp "$READS" "$SPAG_RUN/$SAMPLE.fastq"
    fi
  else
    echo "Reusing $SPAG_RUN/$SAMPLE.fastq"
  fi

  # The trailing slash is required — the runner concatenates $1 with a glob.
  python3 jobs/bench_time.py \
    --csv "$TIMINGS" --label "spaghetti.$SAMPLE" -- \
    bash jobs/framework_16S/framework_16S.spaghetti.sh "$SPAG_RUN/" \
    2>&1 | tee "$LOG_DIR/spaghetti.$SAMPLE.log"

  # Join the phyloseq pair into the genus x species x readn table the scorer
  # reads. Part 1 wrote this to data/comparators/spaghetti.tsv, which is
  # gitignored and would be overwritten per library; Part 2 keeps it beside the
  # run instead.
  python3 jobs/framework_16S/framework_16S.spaghetti2table.py \
    "$SPAG_RUN" "$SPAG_OUT"

  {
    echo "date:         $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "runner:       jobs/framework_16S/framework_16S.spaghetti.sh (upstream ac9d2ec, vendored)"
    echo "porechop:     $(porechop --version 2>&1 | head -1)"
    echo "NanoFilt:     $(NanoFilt --version 2>&1 | head -1)"
    echo "minimap2:     $(minimap2 --version 2>&1 | head -1)"
    if [[ "$SPAG_WHICH" == silva ]]; then
      echo "database:     $SPAG_DB (SILVA 138, QIIME2 release) -- Panel B.1"
    else
      echo "database:     $SPAG_DB (Emu species_taxid.fasta, unchanged) -- Panel B.2"
      echo "db build:     build_spag_emu_db in this script; lineage from Emu taxonomy.tsv"
    fi
    echo "sample:       $SAMPLE"
    echo "reads:        $READS"
    echo "reads_md5:    $(md5 -q "$READS" 2>/dev/null || md5sum "$READS" | cut -d' ' -f1)"
    echo "threads:      $SPAG_THREADS"
    echo "skip_yacrd:   $SKIP_YACRD"
    echo "length_filter: 1200-1800 bp (NanoFilt, upstream default)"
  } > "$SPAG_RUN/spaghetti.provenance.txt"

  # The chopper/NanoFilt window is frozen, so record what it cut rather than
  # retuning it (roadmap step 1: barcode79 lost 12 %).
  if [[ -s "$SPAG_RUN/$SAMPLE-porechop.fastq" && -s "$SPAG_RUN/$SAMPLE-porechop-nanofilt.fastq" ]]; then
    pre=$(( $(wc -l < "$SPAG_RUN/$SAMPLE-porechop.fastq") / 4 ))
    post=$(( $(wc -l < "$SPAG_RUN/$SAMPLE-porechop-nanofilt.fastq") / 4 ))
    printf "reads_porechop\t%d\nreads_nanofilt\t%d\npct_kept\t%.1f\n" \
      "$pre" "$post" "$(echo "scale=4; 100*$post/$pre" | bc)" \
      | tee "$SPAG_RUN/length_filter.tsv"
  fi

  echo "Spaghetti done -> $SPAG_OUT"
 done   # SPAG_WHICH
fi

# --- manifest rows for this library ------------------------------------------
# Columns match jobs/framework_16S/framework_16S.manifest.tsv:
#   label, library, panel, format, taxonomy, path, notes
echo ""
echo "Manifest rows for $SAMPLE:"
if [[ "$TOOL" == "emu" || "$TOOL" == "both" ]]; then
  printf "Emu\t%s\tB.1\temu\tncbi\t%s/emu/%s_rel-abundance.tsv\tEmu defaults, own DB\n" \
    "$SAMPLE" "$OUTROOT" "$SAMPLE"
fi
if [[ "$TOOL" == "spaghetti" || "$TOOL" == "both" ]]; then
  for d in $SPAG_DBS; do
    if [[ "$d" == silva ]]; then
      printf "Spaghetti\t%s\tB.1\ttable\tsilva\t%s/spaghetti.tsv\tgenus only; chimera removal skipped\n" \
        "$SAMPLE" "$OUTROOT"
    else
      printf "Spaghetti.emu\t%s\tB.2\ttable\tncbi\t%s/spaghetti.emu.tsv\tgenus only; Emu reference; chimera removal skipped\n" \
        "$SAMPLE" "$OUTROOT"
    fi
  done
fi

done   # SAMPLE

echo
echo "Results: results/framework_16S/<library>/{emu,spaghetti}/"
echo "Timings: $TIMINGS"
echo "Then:    python3 jobs/framework_16S/framework_16S.score.py --manifest jobs/framework_16S/framework_16S.manifest.tsv"
