#!/usr/bin/env Rscript
# Framework B — Petersen 2026 SRA: per-sample genome recovery (3 pipelines)
#
# Bar HEIGHT = genome coverage (%); the number printed inside each bar (black,
# vertical) = consensus identity (%). x = virus, fill = pipeline.
# Rows as Supplementary Fig. 1a, all scored by ngs-benchmark (rlib/bmark_io.R::load_rescored_petersen):
# EMtrace and EMtrace (chain2) polished with medaka, ViMOP and TELEVIR.
# The unpolished row, EMtrace (nopolish), is left out.
#
# Output: figures/petersen_coverage_identity.png
# Run from the repo root:  Rscript rplots/petersen_coverage_identity.R

source("rlib/bmark_io.R")
suppressPackageStartupMessages(library(ggplot2))

out_dir <- "figures"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- load_rescored_petersen() |>
  filter(Pipeline != "EMtrace (nopolish)") |>
  mutate(pipeline = factor(as.character(Pipeline), levels = PIPELINE_LEVELS),
         virus    = factor(virus, levels = c("Zika", "DENV", "WNV", "SARS-CoV-2")))

dodge <- position_dodge(width = 0.8)
p <- ggplot(dat, aes(virus, coverage, fill = pipeline)) +
  geom_col(position = dodge, width = 0.72) +
  geom_text(aes(y = pmax(coverage / 2, 8),
                label = if_else(is.na(identity), "", sprintf("%.1f", identity))),
            position = dodge, angle = 90, colour = "black", size = 2.8) +
  scale_fill_manual(values = PIPELINE_COLOURS, labels = \(x) sub("^EMtrace", "EMtrace", x)) +
  scale_y_continuous(limits = c(0, 100)) +
  # bar height = genome coverage %; in-bar number = consensus identity %
  labs(x = NULL, y = "Genome coverage (%)", fill = NULL,
       tag = "a") +   # panel label for the manuscript figure
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        plot.tag = element_text(face = "bold", size = 16))

out <- file.path(out_dir, "petersen_coverage_identity.png")
ggsave(out, p, width = 8, height = 4.5, dpi = 150)
cat("wrote", out, "\n")
