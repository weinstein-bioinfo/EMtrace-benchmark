# bmark_io.R — shared loaders/constants for EMtrace Framework B tables & figures.
#
# Sourced by rtables/*.R and rplots/*.R so the three-pipeline
# comparison (EMtrace, ViMOP, TELEVIR) is defined in ONE place. All paths are
# relative to the repo root; run R from the repo root.
#
# Requires: readr, dplyr, tidyr, stringr

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
})

BMARK_METRICS_DIR <- "metrics/framework_B"

# --- pipeline identity (order + colours shared across all figures/tables) -----
# "EMtrace (k11/w5)" is the sensitive-annotation configuration of EMtrace,
# not a separate tool — kept as a distinct series so the default and sensitive
# configs are compared transparently. "EMtrace (chain2)" is the nucleotide +
# DIAMOND chain annotation (vi.chain2); it exists only as a v20260924 batch.
PIPELINE_LEVELS  <- c("EMtrace", "EMtrace (k11/w5)", "EMtrace (chain2)",
                      "EMtrace (nopolish)", "ViMOP", "TELEVIR")
PIPELINE_COLOURS <- c(EMtrace = "#1b9e77", "EMtrace (k11/w5)" = "#66c2a5",
                      "EMtrace (chain2)" = "#e7298a",
                      "EMtrace (nopolish)" = "#a6a6a6",
                      ViMOP = "#d95f02", TELEVIR = "#7570b3")

# --- dataset maps -------------------------------------------------------------
MUTANT_VIRUSES <- c("COVID", "DENV", "LASV_L", "LASV_S")
# Facet labels for the virus codes, as the manuscript text abbreviates them
VIRUS_LABELS <- c(COVID = "SARS-CoV-2", DENV = "DENV-2", LASV_L = "LASV L", LASV_S = "LASV S")
TELEVIR_RATES  <- c(0, 5, 10, 15, 20, 25, 30)   # full mutant rate matrix
PETERSEN_VIRUS <- c(SRR31174237 = "Zika", ERR14941561 = "WNV",
                    SRR15356294 = "SARS-CoV-2", SRR32419728 = "DENV")

# Summary file (in BMARK_METRICS_DIR) for each pipeline × dataset.
SUMMARY_FILES <- list(
  pure     = c(EMtrace = "B6_pure.summary.tsv",
               "EMtrace (k11/w5)" = "B6_pure_k11w5.summary.tsv",
               "EMtrace (chain2)" = "B6_pure_chain2.summary.tsv",
               ViMOP = "vimop_pure.summary.tsv",
               TELEVIR = "televir_pure.summary.tsv"),
  mixed    = c(EMtrace = "B7_mixed.summary.tsv",
               "EMtrace (k11/w5)" = "B7_mixed_k11w5.summary.tsv",
               "EMtrace (chain2)" = "B7_mixed_chain2.summary.tsv",
               ViMOP = "vimop_mixed.summary.tsv",
               TELEVIR = "televir_mixed.summary.tsv"),
  # Petersen: the v20260924 re-run (ref_auto), CELL refauto and chain2, scored
  # against the fixed references (framework_B.sra.metrics.bash). The earlier
  # default and k11/w5 runs (sra_petersen_2026[_k11w5], given --refaccs,
  # pre-5b5266d calling) are superseded.
  petersen = c(EMtrace = "sra_petersen_2026_v20260924_refauto.summary.tsv",
               "EMtrace (chain2)" = "sra_petersen_2026_v20260924_chain2.summary.tsv",
               ViMOP = "vimop_peterson_2026.summary.tsv",
               TELEVIR = "televir_peterson_2026.summary.tsv"),
  # Gauthier: the v20260924 re-run (ref_auto, medaka-polished) and the comparators,
  # all scored against Wuhan-Hu-1 NC_045512.2 by framework_B.rescore.bash. The
  # earlier run (sra_gauthier_2021, given reference, 5x) and the OX637002.1
  # scoring are superseded.
  gauthier = c(EMtrace = "fwB.gauthier.EMtrace.summary.tsv",
               ViMOP = "fwB.gauthier.ViMOP.summary.tsv",
               TELEVIR = "fwB.gauthier.TELEVIR.summary.tsv")
)

# read one summary -> tibble(sample, coverage, identity, pipeline) ------------
read_pipeline_summary <- function(file, pipeline) {
  path <- file.path(BMARK_METRICS_DIR, file)
  if (!file.exists(path)) { warning("missing summary, skipped: ", path); return(NULL) }
  read_tsv(path, show_col_types = FALSE) |>
    transmute(sample   = as.character(sample),
              coverage = as.numeric(mean_coverage_pct),
              identity = as.numeric(mean_identity_pct),
              pipeline = pipeline) |>
    # identity is undefined without called bases; the scorer writes 0.0 there
    # (e.g. ViMOP Gauthier, TELEVIR at 10x). NA keeps such samples out of the
    # identity means, as for samples with no consensus at all (Methods, Assembly).
    mutate(identity = if_else(!is.na(coverage) & coverage == 0, NA_real_, identity))
}

# Full re-runs, scored as separate batches <label>_<version>.summary.tsv
# (framework_B.metrics.bash B6_pure_v20260924 etc.); value = facet label.
RUN_VERSIONS <- c(v20260924 = "v20260924 (10×)")

# EMtrace summary scored at one assembly.min_coverage (framework_B.metrics.bash
# with MIN_COVERAGE=N writes <label>.dN.summary.tsv); NULL = the unsuffixed file.
# A RUN_VERSIONS name ("v20260924") selects <label>_v20260924.summary.tsv instead;
# a config without a batch at that version is then missing and skipped (chain2
# has only the v20260924 batch, so it drops from the unsuffixed and dN figures).
mc_file <- function(file, mc) {
  if (is.null(mc)) return(file)
  # A run version is inserted after the batch stem, so the config tag stays last:
  # B6_pure.summary.tsv -> B6_pure_v20260924.summary.tsv
  # B6_pure_k11w5.summary.tsv -> B6_pure_v20260924_k11w5.summary.tsv
  # (this is how framework_B.metrics.bash labels the batches)
  if (mc %in% names(RUN_VERSIONS)) return(sub("^(B[67]_(pure|mixed))", paste0("\\1_", mc), file))
  sub("\\.summary\\.tsv$", paste0(".d", mc, ".summary.tsv"), file)
}

# facet label for a threshold argument: "10" -> "min_coverage 10×",
# "10.ontsup" (a tagged scoring, METRICS_TAG) -> "min_coverage 10× (ontsup)",
# "v20260924" -> RUN_VERSIONS label
mc_label <- function(mc) {
  n <- sub("\\..*$", "", mc); tag <- ifelse(grepl("\\.", mc), sub("^[^.]*\\.", "", mc), "")
  lab <- paste0("min_coverage ", n, "×", ifelse(tag == "", "", paste0(" (", tag, ")")))
  unname(ifelse(mc %in% names(RUN_VERSIONS), RUN_VERSIONS[mc], lab))
}

# output-file tag for a threshold argument: "10" -> "d10", "v20260924" -> "v20260924"
mc_tag <- function(mc) ifelse(mc %in% names(RUN_VERSIONS), mc, paste0("d", mc))

# summary files for a dataset, EMtrace's switched to min_coverage mc ------
pipeline_files <- function(dataset, mc = NULL) {
  files <- SUMMARY_FILES[[dataset]]
  emtrace  <- startsWith(names(files), "EMtrace")
  files[emtrace] <- vapply(files[emtrace], mc_file, "", mc = mc)
  files
}

# bind all three pipelines for a dataset; mc selects EMtrace's min_coverage --
load_dataset <- function(dataset, mc = NULL) {
  files <- pipeline_files(dataset, mc)
  bind_rows(lapply(names(files), \(p) read_pipeline_summary(files[[p]], p))) |>
    mutate(pipeline = factor(pipeline, levels = PIPELINE_LEVELS))
}

# read one <label>.coverage.tsv, summing base counts over segments per sample --
read_pipeline_counts <- function(file, pipeline) {
  path <- file.path(BMARK_METRICS_DIR, sub("\\.summary\\.tsv$", ".coverage.tsv", file))
  if (!file.exists(path)) { warning("missing coverage table, skipped: ", path); return(NULL) }
  read_tsv(path, col_types = cols(sample = "c", .default = "c")) |>
    mutate(across(c(ref_len, covered_bp, aligned_bp, matching_bp), as.numeric)) |>
    group_by(sample) |>
    summarise(across(c(ref_len, covered_bp, aligned_bp, matching_bp), \(x) sum(x, na.rm = TRUE)),
              .groups = "drop") |>
    mutate(pipeline = pipeline)
}

# Mutants for one condition, all pipelines: per-sample base counts, coverage (%)
# and substitution errors per 10 kbp of aligned consensus. The scorer excludes
# gap and N columns (ngs_benchmark/assembly/mafft.py), so indels are NOT
# counted. Samples with no aligned bases (no consensus) get errors NA.
load_mutant_counts <- function(condition = c("pure", "mixed"), mc = NULL) {
  condition <- match.arg(condition)
  files <- pipeline_files(condition, mc)
  bind_rows(lapply(names(files), \(p) read_pipeline_counts(files[[p]], p))) |>
    mutate(pipeline = factor(pipeline, levels = PIPELINE_LEVELS)) |>
    parse_mutant() |>
    filter(!is.na(rate)) |>
    mutate(virus    = factor(virus, levels = MUTANT_VIRUSES),
           coverage = 100 * covered_bp / ref_len,
           errors   = if_else(aligned_bp > 0,
                              1e4 * (aligned_bp - matching_bp) / aligned_bp, NA_real_))
}

# sample id -> virus + mutation rate (strip "_mixed", trailing NN) ------------
parse_mutant <- function(df) {
  df |> mutate(base  = str_remove(sample, "_mixed$"),
               rate  = as.integer(str_extract(base, "\\d{2}$")),
               virus = str_remove(base, "_\\d{2}$"))
}

# Mutants for one condition, all pipelines, with TELEVIR's no-consensus
# samples 0-filled (coverage 0, identity NA) across the full 4-virus x 7-rate
# grid so failures are represented rather than silently dropped from the mean.
load_mutants <- function(condition = c("pure", "mixed"), mc = NULL) {
  condition <- match.arg(condition)
  dat <- load_dataset(condition, mc) |> parse_mutant() |>
    filter(!is.na(rate)) |>
    select(pipeline, virus, rate, coverage, identity)

  tel <- dat |> filter(pipeline == "TELEVIR") |>
    complete(virus = MUTANT_VIRUSES, rate = TELEVIR_RATES,
             fill = list(coverage = 0)) |>          # identity stays NA
    mutate(pipeline = factor("TELEVIR", levels = PIPELINE_LEVELS))

  # EMtrace with ref_auto emits no consensus when no reference is selected
  # (coverage NA in the summary) -> coverage 0, like TELEVIR's failures
  bind_rows(dat |> filter(pipeline != "TELEVIR"), tel) |>
    mutate(coverage = coalesce(coverage, 0)) |>
    mutate(virus = factor(virus, levels = MUTANT_VIRUSES))
}

# 16 Gauthier patient accessions (Ct present); joins Ct for figures ------------
gauthier_ct <- function() {
  read_tsv(file.path(BMARK_METRICS_DIR, "gauthier_ct_vs_coverage.tsv"),
           show_col_types = FALSE) |>
    transmute(sample = as.character(sample), ct = suppressWarnings(as.numeric(ct)))
}
gauthier_patients <- function() gauthier_ct() |> filter(!is.na(ct)) |> pull(sample)

# --- markdown pipe-table writer (no caption, no trailing <br>) ---------------
# align: character vector "l"/"r"/"c" per column.
md_table <- function(df, align) {
  cell <- function(v) ifelse(is.na(v), "", as.character(v))
  hdr  <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep  <- paste0("|", paste0(vapply(align, \(a)
           switch(a, l = ":--", r = "--:", c = ":-:", "---"), ""), collapse = "|"), "|")
  rows <- apply(df, 1, \(r) paste0("| ", paste(cell(r), collapse = " | "), " |"))
  c(hdr, sep, rows)
}

# write a table as both .csv and .md fragment under tables/ --------
write_table <- function(df, id, align, tables_dir = "tables") {
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(df, file.path(tables_dir, paste0(id, ".csv")))
  writeLines(md_table(df, align), file.path(tables_dir, paste0(id, ".md")))
  message("wrote ", id, ".{csv,md}")
}

# --- rescored Framework B (Tables 4-6, STables 4-6, Figures 2-4 "medaka") ----
# Every row is scored by ngs-benchmark assembly-batch (MAFFT G-INS-i,
# --adjustdirection; per-segment coverage, identity, and insertions/deletions as
# gap runs), run by jobs/framework_B/framework_B.rescore.bash into
#   BMARK_METRICS_DIR/fwB.<set>.<row>.coverage.tsv   (set: pure | mixed | petersen)
# Segments are summed per sample: coverage = covered_bp / ref_len, identity =
# matching_bp / aligned_bp, substitutions = aligned_bp - matching_bp. Every
# pipeline is placed on the full sample grid: a sample without consensus counts
# as coverage 0, identity/indels/substitutions NA (Methods, Assembly).
MUT_ROWS <- tibble::tribble(
  ~Pipeline,               ~label,
  "EMtrace",            "EMtrace",
  "EMtrace (chain2)",   "EMtrace.chain2",
  "EMtrace (nopolish)", "EMtrace.nopolish",
  "ViMOP",                 "ViMOP",
  "TELEVIR",               "TELEVIR"
)

# per-sample sums of one fwB.<set>.<label>.coverage.tsv
read_fwb <- function(set, label) {
  path <- file.path(BMARK_METRICS_DIR, sprintf("fwB.%s.%s.coverage.tsv", set, label))
  if (!file.exists(path)) stop("missing ", path, " -- run jobs/framework_B/framework_B.rescore.bash")
  read_tsv(path, col_types = cols(sample = "c", segment = "c", .default = "d")) |>
    group_by(sample) |>
    summarise(ref_len = sum(ref_len), covered_bp = sum(covered_bp),
              aligned_bp = sum(aligned_bp), matching_bp = sum(matching_bp),
              ins = sum(insertions), del = sum(deletions),
              has_cons = any(!is.na(coverage_pct)), .groups = "drop") |>
    transmute(sample,
              coverage   = 100 * covered_bp / ref_len,
              identity   = if_else(aligned_bp > 0, 100 * matching_bp / aligned_bp, NA_real_),
              indels     = if_else(has_cons & covered_bp > 0, ins + del, NA_real_),
              subs       = if_else(has_cons & covered_bp > 0, aligned_bp - matching_bp, NA_real_),
              covered_bp)
}

load_rescored <- function(condition) {
  scored <- purrr::map2_dfr(MUT_ROWS$Pipeline, MUT_ROWS$label,
                            \(p, l) read_fwb(condition, l) |> mutate(Pipeline = p)) |>
    parse_mutant() |>
    select(Pipeline, virus, rate, coverage, identity, indels, subs, covered_bp)
  expand_grid(Pipeline = MUT_ROWS$Pipeline, virus = MUTANT_VIRUSES,
              rate = c(0, 5, 10, 15, 20, 25, 30)) |>
    left_join(scored, by = c("Pipeline", "virus", "rate")) |>
    mutate(coverage = coalesce(coverage, 0),
           Pipeline = factor(Pipeline, levels = MUT_ROWS$Pipeline))
}

# Petersen 2026 (Supplementary Fig. 1a), same rows and scorer
load_rescored_petersen <- function() {
  scored <- purrr::map2_dfr(MUT_ROWS$Pipeline, MUT_ROWS$label,
                            \(p, l) read_fwb("petersen", l) |> mutate(Pipeline = p))
  expand_grid(Pipeline = MUT_ROWS$Pipeline, sample = names(PETERSEN_VIRUS)) |>
    left_join(scored, by = c("Pipeline", "sample")) |>
    mutate(coverage = coalesce(coverage, 0),
           virus    = unname(PETERSEN_VIRUS[sample]),
           Pipeline = factor(Pipeline, levels = MUT_ROWS$Pipeline))
}
