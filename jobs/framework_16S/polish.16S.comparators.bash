#!/usr/bin/env bash
# ============================================================================
# 16S polishing — R10.4.1 benchmark: data staging + comparators   [local Mac]
# ============================================================================
# Plan: the project notes. ONT open data zymo_16s_2025.09
# (https://epi2me.nanoporetech.com/zymo_16s_2025.09/), R10.4.1, Dorado 1.1.1,
# kit SQK-MAB114.24. One library per replicate:
#
#   community  bucket folder  replicates  library names
#   gutmock    GutMock        4           gutmock_<model>_r1 .. r4   ZymoBIOMICS Gut Microbiome Standard (D6331)
#   msplus     MSPlus         6           msplus_<model>_r1 .. r6    Zymo DNA standard + 4 ATCC species
#
# <model> is sup (default) or hac. Main run: gutmock, sup (§2.4 step 5).
#
# PARTS (--part, comma-separated; default fetch,emu,spaghetti,kraken,wf16s)
#   fetch      download the replicate BAMs (public S3 over HTTPS, no aws CLI)
#              and write one FASTQ per replicate to
#              data/sra/16S/<library>.fastq.gz -- the fallback path both
#              framework_16S.comparators.bash and framework_16S.kraken.bash
#              resolve, so they run unchanged with --library <library>.
#              The FASTQ header comment keeps Dorado's st (start time) and qs
#              tags: st is what a replay batches on.
#   emu        Emu v3.6.2 defaults, own DB                      (Panel B.1 = B.2)
#   spaghetti  Spaghetti on SILVA (B.1) and on Emu's sequences (B.2)
#   kraken     Kraken 2 + Bracken on std16s (B.1) and emu (B.2); databases
#              must already be built (framework_16S.kraken.bash --pipe build)
#   wf16s      ONT wf-16s, minimap2 classifier, its default ncbi_16s_18s
#              database (B.1; the vendor workflow), via Nextflow + Docker.
#              Converted to the scorer's `table` format (wf16s.tsv).
#   manifest   print the manifest rows only
#
# Usage (from the repo root):
#   bash jobs/framework_16S/polish.16S.comparators.bash                       # gutmock sup, everything
#   bash jobs/framework_16S/polish.16S.comparators.bash --part fetch
#   bash jobs/framework_16S/polish.16S.comparators.bash --community msplus --part fetch,emu
#   bash jobs/framework_16S/polish.16S.comparators.bash --reps "1 2" --part emu
#
# Before scoring:
#   - build the truth tables (the project notes step 1) and register the
#     libraries in framework_16S.score.py TRUTHS (exact-key lookup today);
#   - msplus only after the DNA standard's catalogue number (D6305/D6306) is
#     confirmed.
# Then paste the printed rows into framework_16S.manifest.tsv (panel diag
# until the section is written up) and score.
#
# Cost (Mac, 8 threads, ~73-110 k reads per replicate): Emu ~30-45 min and
# ~3 GB, Spaghetti ~45 min on the Emu DB, Kraken 2/Bracken minutes, wf-16s
# ~10 min. Download: ~0.4 GB of BAMs per community and model.
# ============================================================================

set -uo pipefail
[[ -d jobs ]] || { echo "error: run from the repository root" >&2; exit 1; }

COMMUNITY=gutmock
MODEL=sup
REPS=""
PARTS="fetch,emu,spaghetti,kraken,wf16s"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --community) COMMUNITY="$2"; shift 2 ;;
    --model)     MODEL="$2";     shift 2 ;;
    --reps)      REPS="$2";      shift 2 ;;
    --part)      PARTS="$2";     shift 2 ;;
    -h|--help)   sed -n '2,50p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done
case "$COMMUNITY" in
  gutmock) BUCKET_DIR=GutMock; NREPS=4 ;;
  msplus)  BUCKET_DIR=MSPlus;  NREPS=6 ;;
  *) echo "--community must be gutmock|msplus" >&2; exit 2 ;;
esac
case "$MODEL" in sup|hac) ;; *) echo "--model must be sup|hac" >&2; exit 2 ;; esac
[[ -n "$REPS" ]] || REPS="$(seq 1 "$NREPS" | tr '\n' ' ')"
has_part() { [[ ",$PARTS," == *",$1,"* ]]; }

S3="https://ont-open-data.s3.amazonaws.com"
PREFIX="zymo_16s_2025.09/basecalls/$BUCKET_DIR/MAB114/$MODEL"
BAM_DIR="data/ont_open/zymo_16s_2025.09/$BUCKET_DIR/MAB114/$MODEL"
FQ_DIR="data/sra/16S"
LOG_DIR="logs/framework_16S/polish.16S"
TIMINGS="metrics/framework_16S/timings/polish.16S.comparators.timings.csv"
THREADS="${THREADS:-8}"
WF16S_REV="${WF16S_REV:-}"   # pin a wf-16s release, e.g. WF16S_REV=v1.5.0; empty = latest
KRAKEN_ENV="${KRAKEN_ENV:-$HOME/miniconda3/envs/kraken2}"          # conda env with kraken2 + bracken (osx-64)
KRAKEN_DB_ROOT="${KRAKEN_DB_ROOT:-data/comparators/kraken_db}"  # emu/ and std16s/, built by framework_16S.kraken.bash
mkdir -p "$BAM_DIR" "$FQ_DIR" "$LOG_DIR" "$(dirname "$TIMINGS")"

lib_for() { echo "${COMMUNITY}_${MODEL}_r$1"; }

# ---------------------------------------------------------------------------
# fetch: BAM -> FASTQ per replicate
# ---------------------------------------------------------------------------
# One BAM per replicate folder; its name carries the flow cell and run id, so
# it is listed rather than hard-coded.
bam_key() {
  python3 - "$S3" "$PREFIX/rep$1/" <<'LISTPY'
import re, sys, urllib.request
s3, prefix = sys.argv[1:]
x = urllib.request.urlopen(f"{s3}/?list-type=2&prefix={prefix}").read().decode()
keys = [k for k in re.findall(r"<Key>(.*?)</Key>", x) if k.endswith(".bam")]
if len(keys) != 1:
    raise SystemExit(f"expected one BAM under {prefix}, found {len(keys)}")
print(keys[0])
LISTPY
}

fetch_rep() {
  local r="$1" lib key bam fq
  lib="$(lib_for "$r")"
  fq="$FQ_DIR/$lib.fastq.gz"
  if [[ -s "$fq" ]]; then echo "[skip] $lib: $fq exists"; return 0; fi
  key="$(bam_key "$r")" || return 1
  bam="$BAM_DIR/rep$r.bam"
  if [[ ! -s "$bam" ]]; then
    echo "[get ] $key"
    curl -sSf -o "$bam.part" "$S3/$key" && mv "$bam.part" "$bam" || { echo "[FAIL] download $key" >&2; return 1; }
  fi
  echo "[fq  ] $lib <- $bam"
  samtools fastq -T st,qs "$bam" 2>/dev/null | gzip -c > "$fq.part" && mv "$fq.part" "$fq" \
    || { echo "[FAIL] samtools fastq $bam" >&2; return 1; }
  {
    echo "date:       $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "source:     $S3/$key"
    echo "library:    $lib"
    echo "bam_md5:    $(md5 -q "$bam" 2>/dev/null || md5sum "$bam" | cut -d' ' -f1)"
    # The header lists every barcode's read group; the reads carry their own.
    echo "barcode:    $(samtools view "$bam" | head -1 | grep -o 'BC:Z:[^[:space:]]*' | cut -d: -f3)"
    echo "model:      $(samtools view "$bam" | head -1 | grep -o 'RG:Z:[^[:space:]]*' | grep -o 'dna_r10[^_]*_[^_]*_[^_]*_[a-z]*@v[0-9.]*')"
    echo "basecaller: $(samtools view -H "$bam" | awk -F'\t' '/^@PG/ && /dorado/ {print $4, $5; exit}')"
    echo "n_reads:    $(samtools view -c "$bam")"
  } > "$FQ_DIR/$lib.provenance.txt"
  sed 's/^/  /' "$FQ_DIR/$lib.provenance.txt"
}

# ---------------------------------------------------------------------------
# wf16s: ONT's workflow, then convert to the scorer's `table` format
# ---------------------------------------------------------------------------
# abundance_table_species.tsv has one row per lineage: `tax` holds eight
# ;-separated levels (superkingdom;kingdom;phylum;class;order;family;genus;
# species) and the next column the read count. "Unknown"/"<rank>_none" levels
# become empty cells, so an unresolved species scores at genus only.
wf16s_to_table() {
  python3 - "$1" "$2" <<'CONVPY'
import csv, sys
src, dst = sys.argv[1:]
LEVELS = ["division", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
rows = list(csv.DictReader(open(src), delimiter="\t"))
if not rows or "tax" not in rows[0]:
    raise SystemExit(f"{src}: no 'tax' column")
count_col = [c for c in rows[0] if c not in ("tax", "total")][0]
with open(dst, "w") as fh:
    fh.write("division\tspecies\tgenus\tfamily\treadn\n")
    n = 0
    for r in rows:
        lv = dict(zip(LEVELS, r["tax"].split(";")))
        if lv.get("division") in (None, "", "Unclassified"):
            continue
        clean = {k: ("" if (v in ("Unknown", "") or v.endswith("_none")) else v) for k, v in lv.items()}
        fh.write(f"{clean['division']}\t{clean.get('species','')}\t{clean.get('genus','')}\t"
                 f"{clean.get('family','')}\t{int(float(r[count_col]))}\n")
        n += 1
print(f"  wf16s table: {n} lineages -> {dst}")
CONVPY
}

# Per-read assignments -> scorer table. Newer wf-16s releases default to
# taxonomic_rank G and then write only abundance_table_genus.tsv; the per-read
# file (--include_read_assignments) still carries each read's full lineage,
# `C|U  read_id  taxid  ...  division|phylum|class|order|family|genus|species`,
# and counting reads per lineage is what the species table would contain.
wf16s_assign_to_table() {
  python3 - "$1" "$2" <<'CONVPY'
import collections, sys
src, dst = sys.argv[1:]
n = collections.Counter(); unclassified = 0
for line in open(src):
    f = line.rstrip("\n").split("\t")
    if not f or f[0] != "C":
        unclassified += 1; continue
    lv = f[-1].split("|") + [""] * 7
    div, fam, gen, sp = lv[0], lv[4], lv[5], lv[6]
    if not div:   # taxid 0, empty lineage: failed wf-16s's identity/coverage filter
        unclassified += 1; continue
    n[(div, sp, gen, fam)] += 1
with open(dst, "w") as fh:
    fh.write("division\tspecies\tgenus\tfamily\treadn\n")
    for (div, sp, gen, fam), c in n.most_common():
        fh.write(f"{div}\t{sp}\t{gen}\t{fam}\t{c}\n")
print(f"  wf16s table: {len(n)} lineages, {sum(n.values())} classified reads, "
      f"{unclassified} unclassified -> {dst}")
CONVPY
}

run_wf16s() {
  local lib="$1" fq="$FQ_DIR/$1.fastq.gz" out="results/framework_16S/$1/wf16s"
  local table="results/framework_16S/$1/wf16s.tsv"
  command -v nextflow > /dev/null || { echo "[WARN] nextflow not on PATH — skipping wf16s" >&2; return 0; }
  docker info > /dev/null 2>&1 || { echo "[WARN] Docker not running — skipping wf16s" >&2; return 0; }
  if [[ -s "$table" ]]; then echo "[skip] wf16s $lib ($table exists)"; return 0; fi
  local assign; assign="$(find "$out" -name '*lineages.minimap2.assignments.tsv' 2>/dev/null | head -1)"
  if [[ -s "$assign" ]]; then
    # a finished run whose species table is missing: convert, do not re-run
    echo "[conv] wf16s $lib from per-read assignments ($assign)"
    wf16s_assign_to_table "$assign" "$table" || return 1
    echo "table:     $table built from $assign (per-read lineages)" >> "$out/wf16s.provenance.txt"
    return 0
  fi
  mkdir -p "$out"
  local rev=(); [[ -n "$WF16S_REV" ]] && rev=(-r "$WF16S_REV")
  echo "[run ] wf-16s $lib  (minimap2, ncbi_16s_18s${WF16S_REV:+, $WF16S_REV})"
  python3 jobs/bench_time.py --csv "$TIMINGS" --label "wf16s.$lib" -- \
    nextflow run epi2me-labs/wf-16s ${rev[@]+"${rev[@]}"} \
      --fastq "$fq" \
      --classifier minimap2 \
      --database_set ncbi_16s_18s \
      --taxonomic_rank S \
      --include_read_assignments \
      --threads "$THREADS" \
      --out_dir "$out" \
      -profile standard \
    > "$LOG_DIR/wf16s.$lib.log" 2>&1
  local abund; abund="$(find "$out" -name 'abundance_table_species.tsv' | head -1)"
  [[ -s "$abund" ]] || { echo "[FAIL] wf16s $lib: no abundance_table_species.tsv — see $LOG_DIR/wf16s.$lib.log" >&2; return 1; }
  wf16s_to_table "$abund" "$table" || return 1
  {
    echo "date:      $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "workflow:  epi2me-labs/wf-16s ${WF16S_REV:-<latest>} ($(nextflow info epi2me-labs/wf-16s 2>/dev/null | awk '/revision/ {print $NF; exit}'))"
    echo "classifier: minimap2; database_set ncbi_16s_18s; taxonomic_rank S; workflow defaults otherwise"
    echo "reads:     $fq"
  } > "$out/wf16s.provenance.txt"
}

# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
if has_part fetch; then
  command -v samtools > /dev/null || { echo "samtools not on PATH" >&2; exit 1; }
  for r in $REPS; do fetch_rep "$r" || exit 1; done
fi

for r in $REPS; do
  LIB="$(lib_for "$r")"
  if [[ ! -s "$FQ_DIR/$LIB.fastq.gz" ]] && ! has_part manifest; then
    echo "[WARN] $LIB: $FQ_DIR/$LIB.fastq.gz missing — run --part fetch first"; continue
  fi
  echo; echo "=== $LIB ==="
  if has_part emu; then
    THREADS="$THREADS" bash jobs/framework_16S/framework_16S.comparators.bash \
      --library "$LIB" --tool emu || echo "[WARN] Emu failed on $LIB"
  fi
  if has_part spaghetti; then
    THREADS="$THREADS" bash jobs/framework_16S/framework_16S.comparators.bash \
      --library "$LIB" --tool spaghetti --spag-db "silva emu" || echo "[WARN] Spaghetti failed on $LIB"
  fi
  if has_part kraken; then
    # framework_16S.kraken.bash exits unless kraken2, kraken2-build, bracken and
    # bracken-build are on PATH. They live in the osx-64 conda env of that
    # script's header; append (not prepend) its bin so the env's python does
    # not shadow the python3 that bench_time.py and kreport2table.py use.
    if [[ ! -x "$KRAKEN_ENV/bin/kraken2" || ! -x "$KRAKEN_ENV/bin/bracken" ]]; then
      echo "[WARN] kraken2/bracken not found in $KRAKEN_ENV (set KRAKEN_ENV) — skipping Kraken 2 on $LIB"
    elif [[ ! -f "$KRAKEN_DB_ROOT/emu/hash.k2d" || ! -f "$KRAKEN_DB_ROOT/std16s/hash.k2d" ]]; then
      echo "[WARN] Kraken 2 databases not built under $KRAKEN_DB_ROOT — run framework_16S.kraken.bash --pipe build"
    else
      PATH="$PATH:$KRAKEN_ENV/bin" DB_ROOT="$KRAKEN_DB_ROOT" \
        bash jobs/framework_16S/framework_16S.kraken.bash \
          --library "$LIB" --pipe run --db "emu std16s" --threads "$THREADS" \
        || echo "[WARN] Kraken 2 failed on $LIB — see logs/framework_16S/kraken2.*.$LIB.log"
    fi
  fi
  if has_part wf16s; then
    run_wf16s "$LIB" || echo "[WARN] wf-16s failed on $LIB"
  fi
done

# --- manifest rows -----------------------------------------------------------
# Columns: label, library, panel, format, taxonomy, path, notes. Panel diag
# until the section is written up (the project notes).
echo
echo "Manifest rows (panel diag):"
for r in $REPS; do
  LIB="$(lib_for "$r")"; R="results/framework_16S/$LIB"
  printf "Emu\t%s\tdiag\temu\tncbi\t%s/emu/%s_rel-abundance.tsv\tR10.4.1 (the project notes): Emu defaults, own DB\n" "$LIB" "$R" "$LIB"
  printf "Spaghetti\t%s\tdiag\ttable\tsilva\t%s/spaghetti.tsv\tR10.4.1: genus only; chimera removal skipped\n" "$LIB" "$R"
  printf "Spaghetti.emu\t%s\tdiag\ttable\tncbi\t%s/spaghetti.emu.tsv\tR10.4.1: Emu reference; chimera removal skipped\n" "$LIB" "$R"
  printf "Kraken2.std16s\t%s\tdiag\ttable\tncbi\t%s/kraken2.std16s.tsv\tR10.4.1: Kraken 2 defaults, std16s DB\n" "$LIB" "$R"
  printf "Bracken.std16s\t%s\tdiag\ttable\tncbi\t%s/bracken.std16s.tsv\tR10.4.1: Bracken -r 1500 -l S -t 10, std16s DB\n" "$LIB" "$R"
  printf "Kraken2.emu\t%s\tdiag\ttable\tncbi\t%s/kraken2.emu.tsv\tR10.4.1: Kraken 2 defaults, emu DB\n" "$LIB" "$R"
  printf "Bracken.emu\t%s\tdiag\ttable\tncbi\t%s/bracken.emu.tsv\tR10.4.1: Bracken -r 1500 -l S -t 10, emu DB\n" "$LIB" "$R"
  printf "wf16s\t%s\tdiag\ttable\tncbi\t%s/wf16s.tsv\tR10.4.1: ONT wf-16s, minimap2, ncbi_16s_18s\n" "$LIB" "$R"
done
echo
echo "EMtrace cells: sbatch jobs/framework_16S/polish.16S.emtrace.bash (Roihu)."
echo "Stage the reads there first:"
echo "  rsync -avP $FQ_DIR/${COMMUNITY}_${MODEL}_r*.fastq.gz roihu:/scratch/project_XXXXXXX/emtrace-benchmark/sra/"
