#!/usr/bin/env bash
#
# Step 3 — Score CZ ID (Framework A comparator) predictions against ground
# truth.
#
# Usage:
#   bash jobs/framework_A/czid/03_ngs_benchmark.sh [config.yaml]
#
# Dependencies: python3 + PyYAML, ngs-benchmark (on PATH), emtrace
# Inputs:  results/framework_A/czid/<dataset>/<sample>/<sample>[_N]_<czid_id>_taxon_report.csv
#          data/truth/<dataset>_truth.tsv                        (provided by user)
# Outputs: detection mode — metrics/<dataset>/czid/<sample>/...
#          buddle mode     — metrics/<dataset>/{manifest.tsv,combined_metrics.tsv,
#                             sensitivity_vs_gcml.tsv,summary.json}
# Logs:    logs/03_bmark_czid_<dataset>_<sample>.log
#
# NOTE: `ngs_benchmark.io.adapters.czid.load()` reads the CZ ID Taxon Report
# CSV directly (comma-separated, tax_id renamed to taxid) and resolves the
# missing species/genus/family/order/class/division lineage from `taxid` via
# the shared NCBI taxonomy DB. Set NGS_BMARK_NCBI_TAXONOMY_DB if it isn't at
# the adapter's default path.

set -euo pipefail

# <<< CHANGE THIS to your own EMtrace checkout's src/ dir.
EMTRACE_SRC="${EMTRACE_SRC:-../EMtrace/src}"
export PYTHONPATH="$EMTRACE_SRC${PYTHONPATH:+:$PYTHONPATH}"

# <<< CHANGE THIS if your NCBI taxonomy DB lives elsewhere. Needed by the
# czid adapter's tax_id -> lineage resolution (ngs_benchmark/io/taxonomy.py).
export NGS_BMARK_NCBI_TAXONOMY_DB="${NGS_BMARK_NCBI_TAXONOMY_DB:-$HOME/databases/emtrace/taxonomy/ncbi/taxonomy.db}"

PIPELINE="czid"
NB_PRED_FMT="czid"

CONFIG="${1:-jobs/config.yaml}"
RES_DIR="results/framework_A/$PIPELINE"
MET_ROOT="metrics"
LOG_DIR="logs"

mkdir -p "$MET_ROOT" "$LOG_DIR"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config not found: $CONFIG" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Read ngs_benchmark scoring parameters
# ---------------------------------------------------------------------------
eval "$(python3 - "$CONFIG" <<'PY'
import sys, yaml, shlex
with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f) or {}
nb = cfg.get("ngs_benchmark", {}) or {}
def emit(name, val):
    if val is not None:
        print(f"export {name}={shlex.quote(str(val))}")
emit("NB_RANKS",     nb.get("ranks", "species,genus,family"))
emit("NB_TAXGROUP",  nb.get("taxgroup", "Viruses"))
emit("NB_SCORE_COL", nb.get("score_col", "readn"))
emit("NB_READ_MIN",  nb.get("read_min", 0))
emit("NB_TAIL",      nb.get("tail", 0.0))
emit("NB_TAIL_COL",  nb.get("tail_col", "readn"))
PY
)"

# ---------------------------------------------------------------------------
# Emit dataset rows: name<TAB>mode<TAB>truth_file
# ---------------------------------------------------------------------------
read_datasets() {
    python3 - "$CONFIG" <<'PY'
import sys, yaml
with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)
for name, ds in (cfg.get("datasets") or {}).items():
    mode = ds.get("mode", "detection")
    truth = ds.get("truth_file", "")
    print(f"{name}\t{mode}\t{truth}")
PY
}

# ---------------------------------------------------------------------------
# Emit sample rows for one dataset: sra<TAB>gc_per_ml<TAB>nc_sra
# Skips PLACEHOLDER. Missing gc_per_ml/nc_sra emit as empty fields.
# ---------------------------------------------------------------------------
read_dataset_samples() {
    local ds_name="$1"
    python3 - "$CONFIG" "$ds_name" <<'PY'
import sys, yaml
cfg_path, ds_name = sys.argv[1], sys.argv[2]
with open(cfg_path) as f:
    cfg = yaml.safe_load(f)
ds = (cfg.get("datasets") or {}).get(ds_name, {}) or {}
for sample in (ds.get("samples") or []):
    sra = sample.get("sra")
    if not sra or sra == "PLACEHOLDER":
        continue
    gc = sample.get("gc_per_ml", "")
    nc = sample.get("nc_sra", "")
    if nc == "PLACEHOLDER":
        nc = ""
    print(f"{sra}\t{gc}\t{nc}")
PY
}

# ---------------------------------------------------------------------------
# Resolve a sample's CZ ID Taxon Report CSV. CZ ID appends its own
# sample-id (and an occasional dedup "_N") to the filename, so the exact
# name isn't derivable from config — glob for it instead.
# Prints the path on stdout; prints nothing (and warns on stderr) if zero
# or more than one file matches.
# ---------------------------------------------------------------------------
resolve_pred() {
    local dataset="$1" sra="$2"
    local dir="$RES_DIR/$dataset/$sra"
    local matches=("$dir/${sra}"*_taxon_report.csv)

    if [[ ! -e "${matches[0]}" ]]; then
        echo "WARN: pred missing for $dataset/$sra (no ${sra}*_taxon_report.csv under $dir)" >&2
        return
    fi
    if [[ "${#matches[@]}" -gt 1 ]]; then
        echo "WARN: multiple taxon_report.csv matches for $dataset/$sra; skipping: ${matches[*]}" >&2
        return
    fi
    echo "${matches[0]}"
}

# ===========================================================================
# Mode: detection — per-sample `ngs-benchmark detection`
# ===========================================================================
run_detection_mode() {
    local dataset="$1"
    local truth="$2"

    if [[ ! -f "$truth" ]]; then
        echo "WARN: truth file missing for $dataset ($truth); skipping dataset" >&2
        return
    fi

    while IFS=$'\t' read -r sra _gc _nc; do
        [[ -z "$sra" ]] && continue
        local pred
        pred="$(resolve_pred "$dataset" "$sra")"
        [[ -z "$pred" ]] && continue

        local outdir="$MET_ROOT/$dataset/$PIPELINE/$sra"
        local logfile="$LOG_DIR/03_bmark_${PIPELINE}_${dataset}_${sra}.log"

        mkdir -p "$outdir"
        echo "[$(date +%H:%M:%S)] detection: $dataset/$sra" | tee -a "$logfile"
        ngs-benchmark detection \
            --pred      "$pred" \
            --truth     "$truth" \
            --outdir    "$outdir" \
            --pred-fmt  "$NB_PRED_FMT" \
            --ranks     "$NB_RANKS" \
            --taxgroup  "$NB_TAXGROUP" \
            --score-col "$NB_SCORE_COL" \
            --read-min  "$NB_READ_MIN" \
            --tail      "$NB_TAIL" \
            --tail-col  "$NB_TAIL_COL" \
            2>&1 | tee -a "$logfile"
    done < <(read_dataset_samples "$dataset")
}

# ===========================================================================
# Mode: buddle
# Per-sample --buddle PMR/RPMR threshold mode is a pending feature in
# runner.run_batch (line 220), so it's not yet wired in here — same caveat as
# jobs/misc/03_ngs_benchmark.sh's emtrace version.
# ===========================================================================
run_buddle_mode() {
    local dataset="$1"
    local truth="$2"

    if [[ ! -f "$truth" ]]; then
        echo "WARN: truth file missing for $dataset ($truth); skipping dataset" >&2
        return
    fi

    local outdir="$MET_ROOT/$dataset"
    local manifest="$outdir/manifest.tsv"
    local logfile="$LOG_DIR/03_bmark_${PIPELINE}_${dataset}.log"
    mkdir -p "$outdir"

    # Build this pipeline's rows. Skip samples whose pred CSV is missing.
    local new_rows
    new_rows="$(mktemp)"
    {
        while IFS=$'\t' read -r sra gc nc; do
            [[ -z "$sra" ]] && continue
            local pred
            pred="$(resolve_pred "$dataset" "$sra")"
            [[ -z "$pred" ]] && continue
            printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$PIPELINE" "$sra" "$pred" "$NB_PRED_FMT" "$gc" "$nc"
        done < <(read_dataset_samples "$dataset")
    } > "$new_rows"

    if [[ ! -s "$new_rows" ]]; then
        echo "WARN: no usable $PIPELINE samples for $dataset; leaving shared manifest/combined_metrics untouched" >&2
        rm -f "$new_rows"
        return
    fi

    # Merge into the dataset's shared manifest: drop this pipeline's stale
    # rows (reruns), keep every other pipeline's rows, append the fresh ones.
    {
        printf "pipeline\tsample_id\tpath\tfmt\tgc_per_ml\tnc_sample_id\n"
        if [[ -f "$manifest" ]]; then
            awk -F'\t' -v p="$PIPELINE" 'NR>1 && $1!=p' "$manifest"
        fi
        cat "$new_rows"
    } > "${manifest}.tmp"
    mv "${manifest}.tmp" "$manifest"
    rm -f "$new_rows"

    echo "[$(date +%H:%M:%S)] detection-batch: $dataset ($(($(wc -l < "$manifest") - 1)) samples across all pipelines)" | tee -a "$logfile"
    ngs-benchmark detection-batch \
        --pred-list "$manifest" \
        --truth     "$truth" \
        --outdir    "$outdir" \
        --pred-fmt  "$NB_PRED_FMT" \
        --ranks     "$NB_RANKS" \
        --taxgroup  "$NB_TAXGROUP" \
        --score-col "$NB_SCORE_COL" \
        --read-min  "$NB_READ_MIN" \
        --tail      "$NB_TAIL" \
        --tail-col  "$NB_TAIL_COL" \
        --dilution-series \
        2>&1 | tee -a "$logfile"
}

# ---------------------------------------------------------------------------
# Main loop.
# ---------------------------------------------------------------------------
while IFS=$'\t' read -r dataset mode truth; do
    [[ -z "$dataset" ]] && continue
    case "$mode" in
        buddle)    run_buddle_mode    "$dataset" "$truth" ;;
        detection) run_detection_mode "$dataset" "$truth" ;;
        *)         echo "WARN: unknown mode '$mode' for dataset $dataset; skipping" >&2 ;;
    esac
done < <(read_datasets)

echo "done — metrics under $MET_ROOT/*/ ($PIPELINE rows merged into each dataset's shared manifest/combined_metrics.tsv)"
