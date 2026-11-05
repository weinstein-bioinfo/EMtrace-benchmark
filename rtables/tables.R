#!/usr/bin/env Rscript
# tables.R — Framework B tables (EMtrace / ViMOP / TELEVIR) from the metrics TSVs,
# written as tables/<id>.{csv,md}.
#
# Tables (ids):
#   fwB_mut_pure_cov     Supplementary Table 3, coverage % by divergence (mean of 4 viruses)
#   fwB_mut_mixed_cov    Supplementary Table 3, "Mean, mixed" column
#   fwB_mut_pure_indels  Supplementary Table 3, indels per 10 kbp by divergence
#   fwB_mut_mixed_indels the same, background-embedded reads
#   fwB_petersen         Supplementary Fig. 1a, per-library coverage/identity
#   fwB_petersen_indels  indels per 10 kbp per field library (Results)
#   fwB_gauthier         Gauthier mean coverage/identity over the 16 patients (Results)
#
# Run from the repo root:  Rscript rtables/tables.R

source("rlib/bmark_io.R")

fmt1 <- function(x) ifelse(is.na(x), "—", formatC(round(x, 1), format = "f", digits = 1))

# ============================================================================
# Supplementary Table 3 — mutants: coverage % and indels per 10 kbp by divergence
#
# All rows are scored by ngs-benchmark assembly-batch (MAFFT G-INS-i,
# --adjustdirection; indels counted as gap runs), run by
# jobs/framework_B/framework_B.rescore.bash -> results/metrics/fwB.<set>.<row>.coverage.tsv,
# read by bmark_io.R::load_rescored. coverage = covered_bp / ref_len;
# identity = matching_bp / aligned_bp.
# Every pipeline is scored on the full 4 viruses x 7 rates matrix: a sample
# without consensus (no reference selected, TELEVIR without alignment) counts
# as coverage 0 and is left out of the identity and indel means (note [2]).
# Indels are pooled rates, sum(ins + del) / sum(covered_bp) x 1e4, so a
# pipeline that calls more of the genome is not penalised for it.
# ============================================================================
RATE_COLS <- paste0(c(0, 5, 10, 15, 20, 25, 30), "%")

# MUT_ROWS and load_rescored() live in rlib/bmark_io.R (shared with Fig. 3)

indel_rate <- function(indels, covered_bp) {
  ok <- !is.na(indels)
  if (!any(ok)) NA_real_ else 1e4 * sum(indels[ok]) / sum(covered_bp[ok])
}

# by-rate values (over the 4 viruses) + one overall column, one row per pipeline
rate_wide <- function(m, value, overall) {
  m |> group_by(Pipeline, rate) |>
    summarise(v = {{ value }}, .groups = "drop") |>
    mutate(v = if_else(is.nan(v), NA_real_, v), rate = paste0(rate, "%")) |>
    pivot_wider(names_from = rate, values_from = v) |>
    left_join(overall, by = "Pipeline") |>
    arrange(Pipeline) |>
    mutate(Pipeline = as.character(Pipeline), across(-Pipeline, fmt1))
}

mutant_table <- function(condition, id) {
  m <- load_rescored(condition)
  overall <- m |> group_by(Pipeline) |>
    summarise(`Overall cov`    = mean(coverage),
              `Overall id`     = mean(identity, na.rm = TRUE),
              `Overall indels` = indel_rate(indels, covered_bp), .groups = "drop")
  tbl <- rate_wide(m, mean(coverage), overall)
  write_table(tbl, id, align = c("l", rep("r", ncol(tbl) - 1)))
}

indel_table <- function(condition, id) {
  m <- load_rescored(condition)
  overall <- m |> group_by(Pipeline) |>
    summarise(Overall = indel_rate(indels, covered_bp), .groups = "drop")
  tbl <- rate_wide(m, indel_rate(indels, covered_bp), overall)
  write_table(tbl, id, align = c("l", rep("r", ncol(tbl) - 1)))
}

mutant_table("pure",  "fwB_mut_pure_cov")
mutant_table("mixed", "fwB_mut_mixed_cov")
indel_table("pure",   "fwB_mut_pure_indels")
indel_table("mixed",  "fwB_mut_mixed_indels")

# ============================================================================
# Supplementary Fig. 1a — Petersen: per-sample coverage / identity, cols Zika/WNV/SARS/DENV
# ============================================================================
# Rows as Tables 4/5, all scored by ngs-benchmark (bmark_io.R::load_rescored_petersen).
# Indels per 10 kbp of called bases, per genome.
pet <- load_rescored_petersen() |>
  mutate(cell = paste0(fmt1(coverage), " / ", fmt1(identity)))

pet_overall <- pet |>
  group_by(Pipeline) |>
  summarise(cov = mean(coverage), id = mean(identity, na.rm = TRUE),
            indels = indel_rate(indels, covered_bp), .groups = "drop") |>
  mutate(overall = paste0(fmt1(cov), " / ", fmt1(id)))

virus_order <- c("Zika", "WNV", "SARS-CoV-2", "DENV")
by_virus <- function(df, col) {
  df |> select(Pipeline, virus, v = {{ col }}) |>
    tidyr::pivot_wider(names_from = virus, values_from = v) |>
    select(Pipeline, all_of(virus_order))
}
pet_tbl <- by_virus(pet, cell) |>
  left_join(pet_overall, by = "Pipeline") |>
  arrange(Pipeline) |>
  transmute(Pipeline = as.character(Pipeline), across(all_of(virus_order)),
            `Overall (cov / id)` = overall, `Overall indels` = fmt1(indels))
write_table(pet_tbl, "fwB_petersen",
            align = c("l", rep("r", ncol(pet_tbl) - 1)))

pet_indels <- pet |>
  mutate(rate = if_else(is.na(indels), NA_real_, 1e4 * indels / covered_bp)) |>
  by_virus(rate) |>
  left_join(pet_overall |> select(Pipeline, Overall = indels), by = "Pipeline") |>
  arrange(Pipeline) |>
  mutate(Pipeline = as.character(Pipeline), across(-Pipeline, fmt1))
write_table(pet_indels, "fwB_petersen_indels",
            align = c("l", rep("r", ncol(pet_indels) - 1)))

# ============================================================================
# Gauthier (Results): overall coverage % / identity % over the 16 patients
# ============================================================================
patients <- gauthier_patients()
gau_tbl <- load_dataset("gauthier") |>
  filter(sample %in% patients) |>
  group_by(pipeline) |>
  summarise(`Overall cov %` = mean(coverage),
            `Overall id %`  = mean(identity, na.rm = TRUE), .groups = "drop") |>
  arrange(pipeline) |>
  transmute(Pipeline = as.character(pipeline),
            `Overall cov %` = fmt1(`Overall cov %`),
            `Overall id %`  = fmt1(`Overall id %`))
write_table(gau_tbl, "fwB_gauthier", align = c("l", "r", "r"))

message("all Framework B tables written to tables/")
