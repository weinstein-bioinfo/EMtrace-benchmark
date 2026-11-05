#!/usr/bin/env bash
# ============================================================================
# Framework 16S — Kraken 2 + Bracken, both databases, both libraries  (Mac)
# ============================================================================
#
# Roadmap status rows 11 (Panel B.1) and 14 (Panel B.2). Small enough to run
# here: both databases are 16S-only, so they build in minutes and load in well
# under 2 GB -- nothing like a standard `nt` build.
#
# WHY THIS TOOL IS IN THE BENCHMARK. Kraken 2 is the no-EM baseline. With the
# Emu database it gives four points on one reference:
#     Kraken 2   k-mer LCA, no abundance model
#     Bracken    Bayesian re-estimation on top of the same classification
#     Emu        EM
#     EMtrace EM
# so the EM's contribution is isolated rather than inferred. Expect Kraken 2 to
# report many false positives on a 16S-only database -- that is the documented
# result, not a misconfiguration, and it is why Bracken exists.
#
# DATABASES
#   emu    built from Emu's species_taxid.fasta, unchanged -- the Panel B.2
#          reference. Headers already carry the taxid as their first
#          colon-field, so relabelling for Kraken 2 is mechanical.
#   std16s Kraken 2's own prebuilt 16S database (Greengenes/SILVA/RDP), the
#          Panel B.1 "its own database" arm. Downloaded, not built.
#
# TWO BUILD DETAILS THAT FAIL SILENTLY IF MISSED
#   - Bracken must be built with -l 1500 (full-length 16S). The default -l 100
#     builds a k-mer distribution for short reads and produces plausible,
#     wrong abundances with no error.
#   - Kraken 2 needs |kraken:taxid| in the FASTA headers, and the taxonomy it
#     was built against must be the one the scorer uses. The emu build reuses
#     the taxdump EMtrace is built against (NCBI_TAXDUMP) rather than
#     downloading a fresh one, so "same database" does not quietly become
#     "same sequences, different nomenclature".
#   - kraken2-build --download-taxonomy is NOT used: it rsyncs to NCBI, which
#     now refuses connections, and it pulls multi-GB accession2taxid maps this
#     database has no use for (every header carries an explicit taxid).
#
# Usage (from the repo root):
#   bash jobs/framework_16S/framework_16S.kraken.bash --pipe all
#
# NCBI_TAXDUMP=<dir with nodes.dmp/names.dmp> overrides the taxdump; if that
# directory is absent the script fetches taxdump.tar.gz over HTTPS.
#   bash ... --pipe build            # databases only
#   bash ... --pipe run --db emu     # classify only, one database
#   bash ... --library barcode79
#
# Outputs (roadmap section 9):
#   results/framework_16S/<library>/kraken2.<db>/…kreport
#   results/framework_16S/<library>/bracken.<db>/…kreport
#   results/framework_16S/<library>/{kraken2,bracken}.<db>.tsv   scorer input
#   plus a provenance file per tool recording versions and thresholds
#
# Prints the manifest rows to paste into framework_16S.manifest.tsv.
#
# Install (bioconda has no osx-arm64 build for kraken2; force osx-64):
#   CONDA_SUBDIR=osx-64 conda create -n kraken2 -c bioconda -c conda-forge \
#       kraken2 bracken
#   conda activate kraken2
# ============================================================================

set -uo pipefail

REPO="$PWD"
DB_ROOT="${DB_ROOT:-data/comparators/kraken_db}"
EMU_FASTA="${EMU_FASTA:-data/comparators/emu_db/species_taxid.fasta}"
# Taxdump reused for the emu build; the same one EMtrace is built against.
NCBI_TAXDUMP="${NCBI_TAXDUMP:-${databases:-$HOME/databases}/emtrace/taxonomy/ncbi}"
LOG_DIR="logs/framework_16S"
TIMINGS="metrics/framework_16S/timings/framework_16S.kraken.timings.csv"
THREADS="${THREADS:-8}"
READLEN="${READLEN:-1500}"     # Bracken k-mer distribution; see note above
BRACKEN_MIN="${BRACKEN_MIN:-10}"   # match the 10-read cutoff every tool is scored at

SAMPLES=""
DBS="emu std16s"
PIPE="all"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --library|--sample) SAMPLES="$2"; shift 2 ;;
    --db)      DBS="$2";     shift 2 ;;
    --pipe)    PIPE="$2";    shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) sed -n '2,50p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
[[ -n "$SAMPLES" ]] || SAMPLES="barcode79 synthgut"
case "$PIPE" in all|build|run) ;; *) echo "--pipe must be all|build|run" >&2; exit 2 ;; esac

for t in kraken2 kraken2-build bracken bracken-build; do
  command -v "$t" >/dev/null || {
    echo "Missing $t. Install the env named in this script's header." >&2; exit 1; }
done
mkdir -p "$DB_ROOT" "$LOG_DIR" "$(dirname "$TIMINGS")"

reads_for() {
  local lib="$1" c
  local -a cands
  case "$lib" in
    barcode79) cands=("data/mock_community.16S/barcode79.fq.gz") ;;
    synthgut)  cands=("data/sra/16S/SRR14326315.fastq"
                      "data/sra/16S/SRR14326315.fastq.gz") ;;
    *)         cands=("data/sra/16S/$lib.fastq"
                      "data/sra/16S/$lib.fastq.gz") ;;
  esac
  for c in "${cands[@]}"; do [[ -f "$c" ]] && { echo "$c"; return 0; }; done
  echo ""
}

# ---------------------------------------------------------------------------
# build
# ---------------------------------------------------------------------------
# Run a build step and stop on failure.
# `cmd | tee log` returns tee's status, so a failed kraken2-build used to pass
# unnoticed and the next step ran against a half-built database. This checks
# PIPESTATUS[0] instead.
run_step() {
  local log="$1"; shift
  "$@" 2>&1 | tee "$log"
  local st=${PIPESTATUS[0]}
  [[ $st -eq 0 ]] || { echo "[FATAL] failed (exit $st): $* — see $log" >&2; return "$st"; }
}

# Kraken 2 needs names.dmp + nodes.dmp. Its own --download-taxonomy uses rsync
# to NCBI, which now refuses connections, and it also pulls the multi-GB
# accession2taxid maps that this database does not need: every sequence carries
# an explicit |kraken:taxid|, so no accession lookup happens.
#
# Prefer the taxdump EMtrace is already built against. Pinning one taxdump
# across Kraken 2 and the EMtrace cells is what keeps "same database" from
# quietly meaning "same sequences, different nomenclature" (roadmap 4.2).
ensure_taxonomy() {
  local db="$1" tax="$1/taxonomy"
  if [[ -f "$tax/nodes.dmp" && -f "$tax/names.dmp" ]]; then
    echo "[skip] taxonomy present ($tax)"; return 0
  fi
  mkdir -p "$tax"
  if [[ -f "$NCBI_TAXDUMP/nodes.dmp" && -f "$NCBI_TAXDUMP/names.dmp" ]]; then
    echo "[taxonomy] copying from $NCBI_TAXDUMP (the EMtrace taxdump)"
    cp "$NCBI_TAXDUMP"/{nodes,names,merged,delnodes}.dmp "$tax/" 2>/dev/null
  else
    echo "[taxonomy] fetching taxdump.tar.gz over HTTPS (rsync to NCBI is refused)"
    curl -fsSL -o "$tax/taxdump.tar.gz" \
      https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz \
      || { echo "[FATAL] taxdump download failed" >&2; return 1; }
    tar -xzf "$tax/taxdump.tar.gz" -C "$tax" nodes.dmp names.dmp merged.dmp delnodes.dmp \
      || { echo "[FATAL] taxdump unpack failed" >&2; return 1; }
  fi
  [[ -f "$tax/nodes.dmp" && -f "$tax/names.dmp" ]] \
    || { echo "[FATAL] no nodes.dmp/names.dmp in $tax" >&2; return 1; }
  echo "[taxonomy] $(wc -l < "$tax/nodes.dmp") nodes"
}

build_emu() {
  local db="$DB_ROOT/emu"
  if [[ -f "$db/hash.k2d" && -f "$db/database${READLEN}mers.kmer_distrib" ]]; then
    echo "[skip] emu database exists ($db)"; return 0
  fi
  [[ -f "$EMU_FASTA" ]] || { echo "Missing $EMU_FASTA" >&2; return 1; }
  echo "[build] emu database from $(grep -c '^>' "$EMU_FASTA") sequences"
  mkdir -p "$db"

  # Emu names sequences "<taxid>:emu_db:<n>"; Kraken 2 wants |kraken:taxid|<id>.
  # Rewriting the header is the whole conversion -- the sequences are untouched,
  # which is the point of Panel B.2.
  if [[ ! -s "$db/emu.kraken.fna" ]]; then
    awk '/^>/ {split(substr($1,2), p, ":"); print ">" substr($1,2) "|kraken:taxid|" p[1]; next} {print}' \
      "$EMU_FASTA" > "$db/emu.kraken.fna"
  fi

  ensure_taxonomy "$db" || return 1

  if [[ ! -s "$db/seqid2taxid.map" ]]; then
    run_step "$LOG_DIR/kraken2-build.emu.add.log" \
      kraken2-build --add-to-library "$db/emu.kraken.fna" --db "$db" --threads "$THREADS" \
      || return 1
  else
    echo "[skip] library already added ($db/seqid2taxid.map)"
  fi

  run_step "$LOG_DIR/kraken2-build.emu.build.log" \
    kraken2-build --build --db "$db" --threads "$THREADS" || return 1

  # -l 1500: full-length 16S. The default 100 would silently mis-weight Bracken.
  run_step "$LOG_DIR/bracken-build.emu.log" \
    bracken-build -d "$db" -t "$THREADS" -k 35 -l "$READLEN" || return 1
}

build_std16s() {
  local db="$DB_ROOT/std16s"
  if [[ -f "$db/hash.k2d" && -f "$db/database${READLEN}mers.kmer_distrib" ]]; then
    echo "[skip] std16s database exists ($db)"; return 0
  fi
  echo "[build] std16s: Kraken 2's own 16S database (Panel B.1 arm)"
  mkdir -p "$db"
  # Kraken 2 ships three 16S builds; Greengenes is the one its docs default to.
  kraken2-build --db "$db" --special greengenes --threads "$THREADS" \
    2>&1 | tee "$LOG_DIR/kraken2-build.std16s.log" || {
      echo "[WARN] --special greengenes failed. Try 'silva' or 'rdp':" >&2
      echo "       kraken2-build --db $db --special silva --threads $THREADS" >&2
      return 1; }
  bracken-build -d "$db" -t "$THREADS" -k 35 -l "$READLEN" \
    2>&1 | tee "$LOG_DIR/bracken-build.std16s.log"
}

if [[ "$PIPE" == all || "$PIPE" == build ]]; then
  for d in $DBS; do
    case "$d" in
      emu)    build_emu    || exit 1 ;;
      std16s) build_std16s || echo "[WARN] std16s unavailable — Panel B.1 arm skipped" ;;
      *) echo "unknown --db '$d'" >&2; exit 2 ;;
    esac
  done
fi
[[ "$PIPE" == build ]] && { echo "Databases under $DB_ROOT/"; exit 0; }

# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------
for SAMPLE in $SAMPLES; do
  READS="$(reads_for "$SAMPLE")"
  if [[ ! -f "$READS" ]]; then
    echo "[WARN] $SAMPLE: reads not found — skipping."; continue
  fi
  OUTROOT="results/framework_16S/$SAMPLE"
  echo; echo "=== $SAMPLE  ($(basename "$READS")) ==="

  for d in $DBS; do
    db="$DB_ROOT/$d"
    [[ -f "$db/hash.k2d" ]] || { echo "[skip] $d: not built"; continue; }

    kout="$OUTROOT/kraken2.$d"; bout="$OUTROOT/bracken.$d"
    mkdir -p "$kout" "$bout"
    kreport="$kout/$SAMPLE.kreport"

    if [[ -f "$kreport" ]]; then
      echo "[skip] kraken2.$d (kreport exists)"
    else
      echo "[run ] kraken2 --db $d"
      python3 jobs/bench_time.py --csv "$TIMINGS" --label "kraken2.$d.$SAMPLE" -- \
        kraken2 --db "$db" --threads "$THREADS" \
                --report "$kreport" --output "$kout/$SAMPLE.kraken" \
                "$READS" 2>&1 | tee "$LOG_DIR/kraken2.$d.$SAMPLE.log"
    fi
    [[ -f "$kreport" ]] || { echo "[WARN] no kreport for $d/$SAMPLE"; continue; }

    # Bracken rewrites the same report with its corrected abundances (-w), so
    # both tools go through one converter and differ only in the input file.
    breport="$bout/$SAMPLE.bracken.kreport"
    if [[ -f "$breport" ]]; then
      echo "[skip] bracken.$d (kreport exists)"
    else
      echo "[run ] bracken --db $d -r $READLEN -l S -t $BRACKEN_MIN"
      python3 jobs/bench_time.py --csv "$TIMINGS" --label "bracken.$d.$SAMPLE" -- \
        bracken -d "$db" -i "$kreport" -o "$bout/$SAMPLE.bracken.tsv" \
                -w "$breport" -r "$READLEN" -l S -t "$BRACKEN_MIN" \
        2>&1 | tee "$LOG_DIR/bracken.$d.$SAMPLE.log"
    fi

    python3 jobs/framework_16S/framework_16S.kreport2table.py \
      "$kreport" "$OUTROOT/kraken2.$d.tsv"
    [[ -f "$breport" ]] && python3 jobs/framework_16S/framework_16S.kreport2table.py \
      "$breport" "$OUTROOT/bracken.$d.tsv"

    # Provenance is not optional: without the DB build and the thresholds, the
    # bar is not defensible in a manuscript.
    {
      echo "date:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
      echo "kraken2:   $(kraken2 --version 2>&1 | head -1)"
      echo "bracken:   $(bracken -v 2>&1 | head -1)"
      echo "database:  $d  ($db)"
      [[ "$d" == emu ]] && echo "db source: $EMU_FASTA (Emu's own, unchanged)"
      echo "reads:     $READS"
      echo "bracken:   -r $READLEN -l S -t $BRACKEN_MIN"
      echo "threads:   $THREADS"
    } > "$OUTROOT/kraken.$d.provenance.txt"
  done

  echo
  echo "Manifest rows for $SAMPLE:"
  for d in $DBS; do
    [[ -f "$OUTROOT/kraken2.$d.tsv" ]] || continue
    panel=$([[ "$d" == emu ]] && echo "B.2" || echo "B.1")
    printf "Kraken2.%s\t%s\t%s\ttable\tncbi\t%s/kraken2.%s.tsv\tKraken 2 defaults, %s DB\n" \
      "$d" "$SAMPLE" "$panel" "$OUTROOT" "$d" "$d"
    printf "Bracken.%s\t%s\t%s\ttable\tncbi\t%s/bracken.%s.tsv\tBracken -r %s -l S -t %s, %s DB\n" \
      "$d" "$SAMPLE" "$panel" "$OUTROOT" "$d" "$READLEN" "$BRACKEN_MIN" "$d"
  done
done

echo
echo "Databases: $DB_ROOT/"
echo "Timings:   $TIMINGS"
echo "Then:      python3 jobs/framework_16S/framework_16S.score.py \\"
echo "             --manifest jobs/framework_16S/framework_16S.manifest.tsv"
