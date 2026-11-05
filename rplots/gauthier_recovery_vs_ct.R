#!/usr/bin/env Rscript
# Gauthier 2021 — SARS-CoV-2 genome recovery vs RT-qPCR Ct (3 pipelines)
#
# Genome coverage vs Ct, one line per pipeline (EMtrace / ViMOP / TELEVIR), as in
# the Framework B divergence figures. Each point is one patient sample; lines join
# samples in Ct order. No x offset (unlike Framework B): the cliff lies between
# Ct 25.4 and 25.5, so any offset would move it per pipeline. Above the cliff the
# lines overlap (coverage differs by < 1 point), so they are layered instead:
# TELEVIR widest at the bottom, then ViMOP, EMtrace thinnest on top.
# Dashed line at the Ct ~25 recovery cliff. Negative controls (no Ct) are excluded.
# Consensus identity is not plotted: 99.8-100 % where the genome is recovered
# (Ct <= 25.4) and computed over <= ~5.6 kbp below that, so it is given in the text.
#
# Output: figures/gauthier_recovery_vs_ct.png
# Run from the repo root:  Rscript rplots/gauthier_recovery_vs_ct.R

source("rlib/bmark_io.R")
suppressPackageStartupMessages(library(ggplot2))

out_dir <- "figures"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ct <- gauthier_ct()

dat <- load_dataset("gauthier") |>
  left_join(ct, by = "sample") |>
  filter(!is.na(ct)) |>                                    # patients only
  mutate(layer = factor(as.character(pipeline), levels = c("TELEVIR", "ViMOP", "EMtrace"))) |>
  arrange(layer)                                           # draw order: bottom to top

WIDTH <- c(TELEVIR = 2.0, ViMOP = 1.3, EMtrace = 0.6)   # line width per layer
SIZE  <- c(TELEVIR = 3.2, ViMOP = 2.3, EMtrace = 1.4)   # point size per layer

p <- ggplot(dat, aes(ct, coverage, colour = pipeline, group = layer)) +
  geom_vline(xintercept = 25, linetype = 2, colour = "grey60") +
  geom_line(aes(linewidth = layer)) +
  geom_point(aes(size = layer)) +
  scale_linewidth_manual(values = WIDTH, guide = "none") +
  scale_size_manual(values = SIZE, guide = "none") +
  scale_colour_manual(values = PIPELINE_COLOURS, labels = \(x) sub("^EMtrace", "EMtrace", x)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "RT-qPCR Ct", y = "Genome coverage (%)", colour = NULL,
       tag = "b") +   # panel label for the manuscript figure
  theme_bw(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.tag = element_text(face = "bold", size = 16))

out <- file.path(out_dir, "gauthier_recovery_vs_ct.png")
ggsave(out, p, width = 8, height = 4.5, dpi = 150)   # same size as Fig. 4a
cat("wrote", out, "\n")
