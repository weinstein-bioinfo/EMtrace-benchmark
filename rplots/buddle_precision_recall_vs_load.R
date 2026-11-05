#!/usr/bin/env Rscript
# Framework A — Buddle 2024 mock community: precision & recall vs viral load
#
# Species-level Precision and Recall as a function of dilution (genome copies/mL)
# for EMtrace / EPI2ME / CZ ID / MEGAN-LR, on DNA and RNA libraries.
# Scored like Fig. 2a,b: 10-read threshold, 0.01 % floor, MEGAN-LR from its
# read-count export with reads above species not counted as a species call.
#
# Layout: rows = metric (Precision, Recall) = the two requested panels;
#         cols = nucleic acid (DNA, RNA) so the two library types don't overplot;
#         colour = pipeline. A small per-pipeline x-dodge (in log space) keeps
#         coincident points (e.g. all pipelines at 0 at low load) visible.
#
# Inputs (repo root):
#   metrics/framework_A_floor/floor_em/buddle_2024_{dna,rna}/combined_metrics.tsv
#   (+ manifest.tsv for gc_per_ml; jobs/framework_A/framework_A.rescore_floor.py)
# Output:
#   figures/buddle_precision_recall_vs_load.png
#
# Run from the repo root:
#   Rscript rplots/buddle_precision_recall_vs_load.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

metrics_dir <- "metrics/framework_A_floor/floor_em"
out_dir     <- "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_buddle <- function(sub, na_type) {
  load <- read_tsv(file.path(metrics_dir, sub, "manifest.tsv"), show_col_types = FALSE) |>
    distinct(pipeline, sample_id, gc_per_ml)
  read_tsv(file.path(metrics_dir, sub, "combined_metrics.tsv"), show_col_types = FALSE) |>
    filter(rank == "species") |>
    left_join(load, by = c("pipeline", "sample_id")) |>
    transmute(nucacid = na_type, gc_per_ml, pipeline, Precision = Pr, Recall = Rc)
}

dodge <- c(EMtrace = 0.84, EPI2ME = 0.94, "CZ ID" = 1.06, "MEGAN-LR" = 1.19)  # log-space x-offset

dat <- bind_rows(
  read_buddle("buddle_2024_dna", "DNA"),
  read_buddle("buddle_2024_rna", "RNA")
) |>
  mutate(
    pipeline = recode(pipeline, emtrace = "EMtrace", epi2me = "EPI2ME", czid = "CZ ID",
                      megan_lr = "MEGAN-LR"),
    gc_x     = gc_per_ml * dodge[pipeline]
  ) |>
  pivot_longer(c(Precision, Recall), names_to = "metric", values_to = "value") |>
  mutate(
    pipeline = factor(pipeline, levels = c("EMtrace", "EPI2ME", "CZ ID", "MEGAN-LR")),
    metric   = factor(metric,   levels = c("Precision", "Recall")),
    nucacid  = factor(nucacid,  levels = c("DNA", "RNA"))
  )

p <- ggplot(dat, aes(gc_x, value, colour = pipeline)) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 2, alpha = 0.85) +
  facet_grid(metric ~ nucacid) +
  scale_x_log10(breaks = c(60, 600, 6000, 60000), labels = c("60", "600", "6k", "60k")) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_manual(values = c(EMtrace = "#1b9e77", EPI2ME = "#d95f02", "CZ ID" = "#7570b3",
                                 "MEGAN-LR" = "#e6ab02")) +
  labs(title = "Buddle 2024 — precision & recall vs viral load",
       x = "Viral load (genome copies / mL)", y = NULL, colour = NULL,
       tag = "c") +   # panel label for the manuscript figure
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 16))

out <- file.path(out_dir, "buddle_precision_recall_vs_load.png")
ggsave(out, p, width = 8, height = 5.5, dpi = 150)
cat("wrote", out, "\n")
