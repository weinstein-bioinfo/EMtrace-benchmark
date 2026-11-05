#!/usr/bin/env Rscript
# Framework B — consensus substitution errors vs sequence divergence (3 pipelines)
#
# Companion to framework_B.coverage_vs_divergence.R: same panels, colours and
# arguments, with the accuracy of the recovered bases on the y axis:
#   errors per 10 kbp = (aligned_bp - matching_bp) / aligned_bp × 10^4
# summed over segments per sample, from <label>[.dN].coverage.tsv.
#
# What is and is not counted:
#   * SUBSTITUTIONS only — the scorer excludes gap and N columns, so indels
#     (ONT's homopolymer errors) do not appear here.
#   * At high divergence most errors are reference bias (uncalled variant sites
#     keep the reference base), not sequencing error.
#   * Samples with no consensus are dropped, not 0-filled: they have coverage 0
#     (shown in the coverage figure) but no accuracy.
#   * The "Pooled" column sums counts over the viruses a pipeline produced a
#     consensus for, so it is weighted by aligned bases.
# Zero errors are drawn at the 0.1 floor, labelled 0. One error in the longest
# genome (~29 kb) is ~0.34, so the floor never hides a non-zero value.
# Dashed line: 1.47 wrong bp per 10 kbp, the bound for majority calling
# at d=10, 10% per-read error, ≥50% rule.
#
# EMtrace's assembly.min_coverage is compared by passing thresholds; each
# needs framework_B.metrics.bash run with MIN_COVERAGE=N (<label>.dN.coverage.tsv).
#
# "medaka" plots the series of Fig. 3b instead (rlib/bmark_io.R::load_rescored,
# all scored by ngs-benchmark): EMtrace and EMtrace (chain2) polished with medaka,
# ViMOP and TELEVIR. That scorer does count indels; they are in STables 4/5, so
# this figure still shows substitutions only.
#
# Output: figures/framework_B.errors_vs_divergence.png         (no args)
#         figures/framework_B.errors_vs_divergence.d10_d20.png (args 10 20)
#         figures/framework_B.errors_vs_divergence.medaka.png  (arg medaka)
# Run from the repo root:  Rscript rplots/framework_B.errors_vs_divergence.R [10 20 ... | medaka]

source("rlib/bmark_io.R")
suppressPackageStartupMessages(library(ggplot2))

out_dir <- "figures"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

FLOOR     <- 0.1    # where zero-error samples are drawn
REF_BOUND <- 1.47   # majority-calling bound: d=10, 10% error, ≥50% rule

mcs <- commandArgs(trailingOnly = TRUE)   # e.g. 10 20; none = unsuffixed tables
rescored <- identical(mcs, "medaka")
if (rescored) mcs <- character(0)
mc_levels <- if (length(mcs)) mc_label(mcs) else ""

load_mc <- function(mc) {
  bind_rows(
    load_mutant_counts("pure",  mc) |> mutate(condition = "Pure"),
    load_mutant_counts("mixed", mc) |> mutate(condition = "Mixed")
  ) |> mutate(min_cov = if (is.null(mc)) "" else mc_label(mc))
}
load_medaka <- function() {
  bind_rows(load_rescored("pure")  |> mutate(condition = "Pure"),
            load_rescored("mixed") |> mutate(condition = "Mixed")) |>
    filter(Pipeline != "EMtrace (nopolish)", !is.na(subs)) |>
    transmute(pipeline = factor(as.character(Pipeline), levels = PIPELINE_LEVELS),
              condition, virus = factor(virus, levels = MUTANT_VIRUSES), rate,
              aligned_bp = covered_bp, matching_bp = covered_bp - subs,
              errors = 1e4 * subs / covered_bp, min_cov = "")
}
dat <- (if (rescored) load_medaka() else if (length(mcs)) bind_rows(lapply(mcs, load_mc)) else load_mc(NULL)) |>
  filter(!is.na(errors))

# As in the coverage figure, only a series covering all 4 viruses is pooled; the
# Lassa-only cells (k11/w5, chain2) would pool over a different set of viruses.
# A series is judged over all its rates, not per rate: samples without a consensus
# are already dropped, so a full series can lack a virus at a single rate.
full_series <- dat |>
  group_by(pipeline, condition, min_cov) |>
  filter(n_distinct(virus) == length(MUTANT_VIRUSES)) |>
  ungroup()
pooled <- full_series |>
  group_by(pipeline, condition, min_cov, rate) |>
  summarise(errors = 1e4 * (sum(aligned_bp) - sum(matching_bp)) / sum(aligned_bp),
            .groups = "drop") |>
  mutate(virus = "Pooled")

offset <- c(EMtrace = -0.8, "EMtrace (k11/w5)" = -0.4, "EMtrace (chain2)" = 0,
            ViMOP = 0.4, TELEVIR = 0.8)
plot_dat <- bind_rows(dat |> select(pipeline, condition, min_cov, virus, rate, errors), pooled) |>
  mutate(
    virus     = factor(virus, levels = c(MUTANT_VIRUSES, "Pooled"),
                       labels = c(VIRUS_LABELS[MUTANT_VIRUSES], "Pooled")),
    condition = factor(condition, levels = c("Pure", "Mixed")),
    min_cov   = factor(min_cov, levels = mc_levels),
    rate_x    = rate + offset[as.character(pipeline)],
    errors_y  = pmax(errors, FLOOR)
  )

p <- ggplot(plot_dat, aes(rate_x, errors_y, colour = pipeline)) +
  geom_hline(yintercept = REF_BOUND, linetype = "dashed", colour = "grey45") +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.6, alpha = 0.85) +
  (if (length(mcs)) facet_grid(condition + min_cov ~ virus) else facet_grid(condition ~ virus)) +
  scale_colour_manual(values = PIPELINE_COLOURS, labels = \(x) sub("^EMtrace", "EMtrace", x)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30)) +
  scale_y_log10(breaks = c(FLOOR, 1, 10, 100, 1000),
                labels = c("0", "1", "10", "100", "1000")) +
  # Dashed: 1.47 substitutions per 10 kbp, expected at 10x depth and 10% per-read error
  labs(x = "Divergence (%)", y = "Substitution errors per 10 kbp", colour = NULL,
       tag = if (rescored) "b") +   # panel label for the manuscript figure
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 16))

tag <- if (rescored) ".medaka" else if (length(mcs)) paste0(".", paste(mc_tag(mcs), collapse = "_")) else ""
out <- file.path(out_dir, paste0("framework_B.errors_vs_divergence", tag, ".png"))
ggsave(out, p, width = 12, height = 2.5 * 2 * max(1, length(mcs)), dpi = 150)
cat("wrote", out, "\n")
