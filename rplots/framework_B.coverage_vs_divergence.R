#!/usr/bin/env Rscript
# Framework B — consensus genome coverage vs sequence divergence (3 pipelines)
#
# Genome coverage (%) vs divergence (mutation rate) for EMtrace / ViMOP / TELEVIR, per
# virus (COVID, DENV, LASV_L, LASV_S) and averaged, for pure and
# background-embedded ("mixed") reads. TELEVIR covers rates 0/10/20/30 only;
# its no-consensus LASV samples are 0-filled (see rlib/bmark_io.R::load_mutants).
#
# EMtrace's assembly.min_coverage can be compared by passing thresholds; each
# needs framework_B.metrics.bash run with MIN_COVERAGE=N (<label>.dN.summary.tsv).
# One facet row per condition × threshold, the comparators repeated in each row.
#
# "medaka" plots the rows of Tables 4/5 instead (rlib/bmark_io.R::load_rescored,
# all scored by ngs-benchmark): EMtrace and EMtrace (chain2) polished with medaka,
# ViMOP and TELEVIR. The unpolished row, EMtrace (nopolish), is left out:
# its coverage is nearly that of EMtrace (Tables 4/5 show both).
#
# Output: figures/framework_B.coverage_vs_divergence.png        (no args)
#         figures/framework_B.coverage_vs_divergence.d10_d20.png  (args 10 20)
#         figures/framework_B.coverage_vs_divergence.medaka.png   (arg medaka)
# Run from the repo root:  Rscript rplots/framework_B.coverage_vs_divergence.R [10 20 ... | medaka]

source("rlib/bmark_io.R")
suppressPackageStartupMessages(library(ggplot2))

out_dir <- "figures"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

mcs <- commandArgs(trailingOnly = TRUE)   # e.g. 10 20; none = unsuffixed summaries
rescored <- identical(mcs, "medaka")
if (rescored) mcs <- character(0)
mc_levels <- if (length(mcs)) mc_label(mcs) else ""

load_mc <- function(mc) {
  bind_rows(
    load_mutants("pure",  mc) |> mutate(condition = "Pure"),
    load_mutants("mixed", mc) |> mutate(condition = "Mixed")
  ) |> mutate(min_cov = if (is.null(mc)) "" else mc_label(mc))
}
load_medaka <- function() {
  bind_rows(load_rescored("pure")  |> mutate(condition = "Pure"),
            load_rescored("mixed") |> mutate(condition = "Mixed")) |>
    filter(Pipeline != "EMtrace (nopolish)") |>
    transmute(pipeline = factor(as.character(Pipeline), levels = PIPELINE_LEVELS),
              condition, virus = as.character(virus), rate, coverage, min_cov = "")
}
dat <- if (rescored) load_medaka() else if (length(mcs)) bind_rows(lapply(mcs, load_mc)) else load_mc(NULL)

# The mean facet is only meaningful for a series covering all 4 viruses; a cell
# run on a subset (e.g. the v20260924 k11/w5 cell, Lassa only) is left out of it
# rather than averaged over a different denominator.
mean_dat <- dat |>
  group_by(pipeline, condition, min_cov, rate) |>
  summarise(coverage = mean(coverage), n_virus = n_distinct(virus), .groups = "drop") |>
  filter(n_virus == length(MUTANT_VIRUSES)) |>
  select(-n_virus) |>
  mutate(virus = "Mean (4 viruses)")

offset <- c(EMtrace = -0.8, "EMtrace (k11/w5)" = -0.4, "EMtrace (chain2)" = 0,
            "EMtrace (nopolish)" = -0.4, ViMOP = 0.4, TELEVIR = 0.8)
plot_dat <- bind_rows(dat |> select(pipeline, condition, min_cov, virus, rate, coverage), mean_dat) |>
  mutate(
    virus     = factor(virus,
                       levels = c(MUTANT_VIRUSES, "Mean (4 viruses)"),
                       labels = c(VIRUS_LABELS[MUTANT_VIRUSES], "Mean (4 viruses)")),
    condition = factor(condition, levels = c("Pure", "Mixed")),
    min_cov   = factor(min_cov, levels = mc_levels),
    rate_x    = rate + offset[as.character(pipeline)]
  )

p <- ggplot(plot_dat, aes(rate_x, coverage, colour = pipeline)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 1.6, alpha = 0.85) +
  (if (length(mcs)) facet_grid(condition + min_cov ~ virus) else facet_grid(condition ~ virus)) +
  scale_colour_manual(values = PIPELINE_COLOURS, labels = \(x) sub("^EMtrace", "EMtrace", x)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Divergence (%)", y = "Genome coverage (%)", colour = NULL,
       tag = if (rescored) "a") +   # panel label for the manuscript figure
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 16))

tag <- if (rescored) ".medaka" else if (length(mcs)) paste0(".", paste(mc_tag(mcs), collapse = "_")) else ""
out <- file.path(out_dir, paste0("framework_B.coverage_vs_divergence", tag, ".png"))
ggsave(out, p, width = 12, height = 2.5 * 2 * max(1, length(mcs)), dpi = 150)
cat("wrote", out, "\n")
