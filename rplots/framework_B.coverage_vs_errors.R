#!/usr/bin/env Rscript
# Framework B — coverage vs substitution errors per sample (supplementary)
#
# The trade-off behind the min_coverage choice: one point per mutant sample,
# genome coverage (%) on x, substitution errors per 10 kbp on y (log). A
# pipeline that masks more reports only the deep, easy regions and looks more
# accurate; plotting both axes together makes that visible.
#   colour — pipeline
#   shape  — EMtrace min_coverage threshold (comparators: one shape)
#   size   — mutation rate
#   grey segments join each EMtrace sample across thresholds, so the
#   recovery gained by a lower threshold can be read against the errors it adds.
# Errors are substitutions only (the scorer excludes gap and N columns); samples
# with no consensus are dropped. Zero errors are drawn at the 0.1 floor.
# Metric definitions: framework_B.errors_vs_divergence.R.
#
# Output: figures/framework_B.coverage_vs_errors.png         (no args)
#         figures/framework_B.coverage_vs_errors.d10_d20.png (args 10 20)
# Run from the repo root:  Rscript rplots/framework_B.coverage_vs_errors.R [10 20 ...]

source("rlib/bmark_io.R")
suppressPackageStartupMessages(library(ggplot2))

out_dir <- "figures"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

FLOOR     <- 0.1
REF_BOUND <- 1.47

mcs <- commandArgs(trailingOnly = TRUE)
mc_levels <- if (length(mcs)) paste0(mcs, "×") else "EMtrace"
COMPARATOR <- "comparator"

load_mc <- function(mc) {
  bind_rows(
    load_mutant_counts("pure",  mc) |> mutate(condition = "Pure"),
    load_mutant_counts("mixed", mc) |> mutate(condition = "Mixed")
  ) |> mutate(min_cov = if (is.null(mc)) mc_levels else paste0(mc, "×"))
}
dat <- (if (length(mcs)) bind_rows(lapply(mcs, load_mc)) else load_mc(NULL)) |>
  filter(!is.na(errors)) |>
  mutate(emtrace    = startsWith(as.character(pipeline), "EMtrace"),
         min_cov = if_else(emtrace, min_cov, COMPARATOR)) |>
  distinct(pipeline, condition, min_cov, sample, .keep_all = TRUE) |>   # comparators load once per mc
  mutate(condition = factor(condition, levels = c("Pure", "Mixed")),
         min_cov   = factor(min_cov, levels = c(mc_levels, COMPARATOR)),
         errors_y  = pmax(errors, FLOOR))

p <- ggplot(dat, aes(coverage, errors_y)) +
  geom_hline(yintercept = REF_BOUND, linetype = "dashed", colour = "grey45") +
  (if (length(mcs) > 1)
     geom_path(data = dat |> filter(emtrace) |> arrange(min_cov),
               aes(group = interaction(pipeline, sample)), colour = "grey60", linewidth = 0.4)) +
  geom_point(aes(colour = pipeline, shape = min_cov, size = rate), alpha = 0.8) +
  facet_grid(condition ~ virus) +
  scale_colour_manual(values = PIPELINE_COLOURS) +
  scale_shape_manual(values = setNames(c(16, 17, 15, 18, 3)[seq_along(mc_levels)], mc_levels) |>
                       c(setNames(4, COMPARATOR))) +
  scale_size_continuous(range = c(1, 3.5), breaks = c(0, 10, 20, 30)) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_log10(breaks = c(FLOOR, 1, 10, 100, 1000),
                labels = c("0", "1", "10", "100", "1000")) +
  labs(title = "Framework B — genome coverage vs consensus substitution errors",
       subtitle = "One point per mutant sample; grey segments join a EMtrace sample across min_coverage",
       x = "Genome coverage (%)", y = "Substitution errors per 10 kbp",
       colour = NULL, shape = "min_coverage", size = "Mutation rate (%)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "top", legend.box = "vertical", panel.grid.minor = element_blank())

tag <- if (length(mcs)) paste0(".", paste(mc_tag(mcs), collapse = "_")) else ""
out <- file.path(out_dir, paste0("framework_B.coverage_vs_errors", tag, ".png"))
ggsave(out, p, width = 12, height = 6.5, dpi = 150)
cat("wrote", out, "\n")
