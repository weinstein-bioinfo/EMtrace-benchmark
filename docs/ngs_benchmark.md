# ngs-benchmark

Scoring library of the EMtrace benchmark. Four modules, one CLI (`ngs-benchmark`):

| Module | Commands | Scores |
|:--|:--|:--|
| A — detection accuracy | `detection`, `detection-batch` | Taxon calls against a truth table: precision, recall, F1 per rank; dilution series (Buddle) |
| B — assembly quality | `assembly`, `assembly-batch` | Consensus against its reference: genome coverage, identity, indels |
| C — real-time performance | `realtime` | A MinKNOW run replayed through `emtrace watch`: per-batch accuracy, time to a correct and stable call, latency |
| D — background filtering | `bgfilter`, `bgfilter-batch`, `bgfilter-divbatch` | Host-filter settings: fraction of host and target reads removed |

<br>

Modules A and B take any tool's output through adapters (EMtrace, EPI2ME, CZ ID, MEGAN-LR, metaMix, Kraken 2/Bracken); C and D drive or read EMtrace.

## Install

```bash
pip install -e .
```

## CLI

### Module A — detection accuracy

```bash
ngs-benchmark detection \
  --pred       abundance.tsv \
  --truth      truth.tsv \
  --outdir     results/
```

See `ngs-benchmark detection --help` for the full flag set.

### Module B — assembly quality

Genome coverage (%) and identity (%) of a consensus aligned to its reference.

```bash
ngs-benchmark assembly \
  --consensus  consensus.fasta \
  --reference  reference.fasta \
  --outdir     results/ \
  --candidates annotation/assembly_candidates.tsv   # optional (EMtrace)

# Many consensus/reference pairs from a TSV manifest
#   cols: sample_id, consensus, reference, [candidates]
ngs-benchmark assembly-batch \
  --manifest   manifest.tsv \
  --outdir     results/
```

See `ngs-benchmark assembly --help` for the full flag set. Requires `mafft` on `PATH`.

## Contents

### Module A — detection accuracy

- §3.1 normalisation (phage / retrovirus filtering, tail-trim, group-by sum)
- §3.2 set metrics (TP, FP, FN, Pr, Rc, F1, mean_accuracy) per `(taxgroup, rank)`;
  `--drop-unknown` leaves "unknown" labels (reads a tool places above a rank, e.g. MEGAN-LR's
  LCA) out of the predicted set at that rank, so they are no taxon call. Read counts and the
  `--min-rel` denominator are unchanged; the Excel summary still lists them
- §3.3 PR curve + AUPR + F-max operating point
- §3.5 Buddle PMR/RPMR threshold mode
- §3.6 dilution-series detection metrics
- Adapters: `emtrace`, `czid`, `epi2me` (Bracken format), `kraken_bracken`
  (Kraken2 kreport + auto-detected Bracken), `megan_lr` (tax-to-count export),
  `metamix` (species-assignments table)
- Truth-table loader (membership + dilution-series variants)
- TSV / JSON writers + legacy-R `summary.xlsx` (SUMMARY tab + per-taxgroup
  `{vi,ba,…}_{true,pred,TP,FN,FP}` detail tabs with FP/FN red-fill highlighting
  and percent-formatted Pr/Rc/F1)

### Module B — assembly quality

- MAFFT G-INS-i alignment of consensus → reference (`--globalpair --maxiterate 1000
  --adjustdirection`); coverage, identity and indels (gap runs) from the alignment
- per-segment `coverage.tsv` (segment, ref_len, covered_bp, aligned_bp,
  matching_bp, coverage_pct, identity_pct) — multi-segment genomes get one row
  per record; unaligned segments report zeros
- `summary.tsv`: unweighted mean and length-weighted coverage/identity across
  segments, plus optional `viral_reads`
- `coverage_pct = covered_bp / ref_len` (breadth over the reference);
  `identity_pct = matching_bp / aligned_bp` (accuracy within aligned blocks)
- edge cases: missing/empty consensus → `NA` metrics; missing `--candidates`
  file → `viral_reads = "no hits"` (annotation selected zero viral reads)
- `assembly-batch` writes per-sample dirs + a combined `combined_summary.tsv`

### Taxonomy DB for comparator adapters

The four comparator adapters that ingest taxid-only output (`epi2me`,
`kraken_bracken`, `megan_lr`, `metamix`) call into
[emtrace.models.taxonomy](https://github.com/weinstein-bioinfo/EMtrace/blob/main/src/emtrace/models/taxonomy.py)
to fill in the canonical rank columns. Install `emtrace` once via
`pip install -e /path/to/EMtrace/` and point the package at a prebuilt
taxonomy SQLite DB:

```bash
export NGS_BMARK_NCBI_TAXONOMY_DB=/path/to/taxonomy.db
```

The default lookup path is
`$EMTRACE_DB_DIR/taxonomy/ncbi/taxonomy.db`.

## Tests

```bash
pip install -e ".[dev]"
pytest -q
```
