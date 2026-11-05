# Running metaMix on Roihu — Framework A

Running metaMix (Morfopoulou & Plagnol 2015).

**metaMix version:** CRAN's last release, **0.3**
(archived 2022-12-14, "issues were not corrected in time"; still downloadable from
CRAN's Archive at
`https://cran.r-project.org/src/contrib/Archive/metaMix/metaMix_0.3.tar.gz`)

## Environment setup (Roihu)

**Layout:**

*Already available if emtrace environment is set:*

| what | path |
|---|---|
| code (this repo) | `/projappl/project_XXXXXXX/EMtrace-benchmark` (`$REPO`) |
| manually-installed CLI tools | `/projappl/project_XXXXXXX/install/<tool>` |
| working root (inputs+results) | `/scratch/project_XXXXXXX/emtrace-benchmark` |
| host genomes | `$hostgenomes` (env var) |

*Needs setup:*

| R packages | `/projappl/project_XXXXXXX/project_rpackages` (`$R_LIBS`) |
| reference DBs | `$databases` |
| metaMix pipeline (.tar) | under /projappl/project_XXXXXXX/metaMix |


### 1. Get the code

Clone this repository to `$PROJAPPL` and install EMtrace as in the README.


### 2. External CLI tools

| tool | expected path | used by |
|---|---|---|
| `seqkit` | `$PROJAPPL/install/seqkit_dir/bin/seqkit` | preprocess (`preprocess/diamond/<dataset>_diamond.sbatch`'s `seqkit fx2tab` read-length step) |
| `diamond` | `$PROJAPPL/install/diamond_dir/diamond` | metaMix DB build (`makedb`) + preprocess (`blastx`) |
| `nextflow` | `$PROJAPPL/install/nextflow/nextflow` | preprocessing (`preprocess/<dataset>_taxprofiler.sbatch` → nf-core/taxprofiler) |

Preprocessing (host removal) goes through nf-core/taxprofiler.

Roihu's `diamond` module is broken (`module load diamond` — default or `/2.2.5` — reports
success but leaves a nonexistent, typo'd path on PATH; `command -v diamond` then finds
nothing). Use the manually-installed binary instead — already wired into the sbatch
scripts via `export PATH="$PROJAPPL/install/diamond_dir:$PATH"`, no module load needed.

The manually-installed `nextflow` needs a working JDK — `module load biojava` first
(also already wired into `preprocess/<dataset>_taxprofiler.sbatch`).


### 4. R environment for metaMix
```bash
module load r-env/461                           
mkdir -p /projappl/project_XXXXXXX/project_rpackages
echo "R_LIBS=/projappl/project_XXXXXXX/project_rpackages" >> ~/.Renviron
```

- Inside R (`R --no-save` after `module load r-env/461`):
  ```r
  install.packages(c("data.table", "Matrix", "gtools", "ggplot2"))
  # Do NOT install.packages("Rmpi") — CSC already bundles a working Rmpi, compiled
  # against r-env's internal MPI
  install.packages("/projappl/project_XXXXXXX/metaMix/metaMix_0.3.tar.gz",
                  repos = NULL, type = "source")
  library(metaMix)   # sanity check
  ```

### 5. Build the database

Default (2026-08): cross-pipeline-consistent DB shared with EMtrace/MEGAN-LR
(viral + archaea/bacteria 1/genus + human, RefSeq release 234), *not* built by
anything in this directory — see `framework_A/db_refseq_abvh_protein.dload.sbatch`
and `.build.sbatch` (one directory up). Run those two (dload, then build) before
step 1 below.

## Pipeline

One file per dataset for the full metaMix run (`buddle_dna` | `buddle_rna` | `gauthier` | `peterson`):

| step | script |
|---|---|---|
| 0 (once) | `framework_A/db_refseq_abvh_protein.dload.sbatch` → `.build.sbatch` |
| 1 (once, per dataset) | `framework_A/metamix/preprocess/<dataset>_taxprofiler.sbatch` — nf-core/taxprofiler qc & host removal, whole dataset in one job |
| 2 (array, per dataset) | `framework_A/metamix/preprocess/diamond/<dataset>_diamond.sbatch` — DIAMOND blastx, custom 10-col outfmt (`qseqid qlen sseqid slen stitle bitscore length pident evalue staxids`) for metaMix 0.3's `blast.default=FALSE` mode |
| 3 (array, per dataset) | `framework_A/metamix/<dataset>.sbatch` → `framework_A/metamix/run.R` — currently run locally, see "Running metaMix: local, not Roihu" |
| 4 (once, all samples; buddle datasets only) | `framework_A/metamix/nc_filter.py` → `.../metamix/<dataset>/<ACC>/<ACC>_metamix_final.tsv` |
| 5 | `ngs-benchmark detection` / `detection-batch` |


```bash
sbatch framework_A/db_refseq_abvh_protein.dload.sbatch           # once, resubmit if it times out (resumable)
sbatch framework_A/db_refseq_abvh_protein.build.sbatch
# wait for it to finish, then, for each dataset (buddle_dna shown):
sbatch framework_A/metamix/preprocess/buddle_dna_taxprofiler.sbatch         # once, whole dataset (nf-core/taxprofiler host removal)
# wait for it to finish, then:
sbatch framework_A/metamix/preprocess/buddle_dna_diamond.sbatch          # array 1-x (DIAMOND blastx)
# wait for the whole array, then:
sbatch framework_A/metamix/buddle_dna.sbatch                     # array 1-x

# once every array task in step 3 has finished (buddle_dna/buddle_rna only
# NC comparison per Buddle et al. 2024:
python3 framework_A/metamix/nc_filter.py \
  --config jobs/config.yaml \
  --results-dir /scratch/project_XXXXXXX/emtrace-benchmark/results/framework_a/metamix

# then, run ngs-benchmark detection mode
```


## Design decisions

- **Custom 10-column DIAMOND outfmt for `blast.default=FALSE`** —
  `qseqid qlen sseqid slen stitle bitscore length pident evalue staxids`; read length
  comes from `qlen` and taxid from `staxids`, so no separate read-length or
  accession2taxid file is needed (see "BLOCKING — resolved 2026-09-23" below).

- **Preprocessing uses nf-core/taxprofiler, matching Buddle et al. 2024's own
  methods.** This preprocessing step is shared with MEGAN-LR.



### BLOCKING — resolved 2026-09-23: use `blast.default=FALSE`

Two problems in 0.3's own `generative.prob()` meant it would not get through a real
DIAMOND-against-abvh.1g sample using `blast.default=TRUE` (the mode `run.R` used to
use):

1. **`cbind.data.frame()` → `unused argument (with = FALSE)` crash in the
   `blast.default=TRUE` branch** — `blast.out.gi<-cbind.data.frame(blast.output, protID)`
   (`R/step1_preprocess.R`) demotes the `fread()` data.tables to a plain `data.frame`,
   so a later data.table-only `[, ..., with=FALSE]` subset fails.
2. **0.3's `blast.default=TRUE` subject-ID parsing assumes the legacy
   `gi|1234|ref|ACCESSION.1|...` pipe-delimited NCBI header format**, not a bare
   accession — for both `gi.or.prot="gi"` and `"prot"`, it re-reads the *entire* raw
   DIAMOND output file with `fread(..., sep="|", select=c(2))` / `select=c(4))`. This
   pipeline's DIAMOND DB (`diamond.refseq.abvh.1g.faa`) uses bare
   `>WP_133517669.1 description` headers (no `gi|...|` wrapper — deprecated by NCBI in
   2016) — DIAMOND's `sseqid` output has zero `|` characters, so this would error
   outright.

**Fix: switched to `generative.prob(blast.default=FALSE, ...)`** — metaMix 0.3's other,
equally-published input mode. Verified this actually resolves
both problems by reading `blast.default=FALSE`'s code path line-by-line (not just the
`TRUE` branch) — it builds its equivalent of `blast.out.gi` (`blast.output.length`)
entirely via `fread()`/data.table row-subsetting, never `cbind.data.frame()`, so it
stays a genuine `data.table` all the way to the `with=FALSE` line; and it takes
`staxids` directly from a custom 10-column outfmt instead of parsing subject-ID
headers at all. 

Changes made:
- `framework_A/metamix/preprocess/diamond/<dataset>_diamond.sbatch` (all 4 datasets):
  DIAMOND's `-f 6` outfmt changed from the classic 12-column default to
  `qseqid qlen sseqid slen stitle bitscore length pident evalue staxids` — the exact
  10 columns, in order, `blast.default=FALSE` expects. The `seqkit fx2tab` read-length
  step was removed (no longer needed — `blast.default=FALSE` takes read length
  straight from DIAMOND's own `qlen` column, confirmed by reading the source's
  `select=c(1,2,7,8,10)` → `read,length,aln,ident,taxonID` mapping).
- `framework_A/metamix/run.R`: `generative.prob()` no longer takes
  `read.length.file`/`accession.taxon.file`/`gi.or.acc` (all specific to
  `blast.default=TRUE`) — just `blast.output.file`, `blast.default=FALSE`,
  `contig.weight.file=1`, `outDir`. `reduce.space()` no longer takes
  `taxon.name.map`/`outDir` (0.3's real signature doesn't have them — it reads
  `outDir` back out of `step1$outDir` automatically). `bayes.model.aver()` is
  unchanged (still takes `taxon.name.map`) but its `poster.prob.thr` is left at 0.3's
  own published default (`0.9`).
- All four `preprocess/diamond/<dataset>_diamond.sbatch` scripts drop DIAMOND hits
  with an empty `staxids` field (`awk -F'\t' '$10 != ""'`) — otherwise they crash
  `generative.prob()` with `'i' and 'j' must not contain NA`.

## Running metaMix: local, not Roihu

metaMix is currently not compatible with Roihu's R environment (`r-env`), so the
metaMix run itself (step 3) is done locally.
Preprocessing (steps 0–2) runs on Roihu as described above.
