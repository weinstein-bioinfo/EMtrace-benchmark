#!/usr/bin/env Rscript
# Framework A — detection accuracy on field and clinical samples (manuscript Fig. 2a,b)
#
# (a) Petersen 2026, four single-virus field libraries: species F1 per library,
#     and pooled (micro-averaged) precision, recall and F1, per tool.
# (b) Gauthier 2021, 16 SARS-CoV-2 patient samples: pooled species precision,
#     recall and F1 per tool, as tiles in the style of (a)'s pooled block (left),
#     and which false-positive species each tool reported, with the number of
#     samples it reported them in (right).
# Panel (c) of the figure is buddle_precision_recall_vs_load.R.
#
# Input: <metrics root>/{peterson_2026,gauthier_2021}/<tool>/<sample>/{metrics,per_taxon}.tsv
#        (ngs-benchmark detection; species rank, 10-read threshold, 0.01 % abundance floor)
#        All four tools: metrics/framework_A_floor/floor_em
#        (jobs/framework_A/framework_A.rescore_floor.py; MEGAN-LR from its read-count
#        export, with reads placed above species not counted as a species call).
# Output: figures/frameworkA_detection.png
# Run from the repo root:  Rscript rplots/frameworkA_detection.R

source("rlib/bmark_io.R")
suppressPackageStartupMessages({ library(ggplot2); library(grid) })

MET <- c(emtrace = "metrics/framework_A_floor/floor_em",
         epi2me     = "metrics/framework_A_floor/floor_em",
         czid       = "metrics/framework_A_floor/floor_em",
         megan_lr   = "metrics/framework_A_floor/floor_em")
TOOLS <- c(emtrace = "EMtrace", epi2me = "EPI2ME", czid = "CZ ID", megan_lr = "MEGAN-LR")
COLS  <- c(EMtrace = "#1b9e77", EPI2ME = "#d95f02", "CZ ID" = "#7570b3", "MEGAN-LR" = "#e6ab02")

read_species <- function(dataset, file, samples = NULL) {
  purrr::map_dfr(names(TOOLS), \(t) {
    dirs <- list.dirs(file.path(MET[[t]], dataset, t), recursive = FALSE)
    if (!is.null(samples)) dirs <- dirs[basename(dirs) %in% samples]
    purrr::map_dfr(dirs, \(d) read_tsv(file.path(d, file), show_col_types = FALSE) |>
                     mutate(sample = basename(d)))  |>
      mutate(tool = TOOLS[[t]])
  }) |>
    filter(rank == "species") |>
    mutate(tool = factor(tool, levels = TOOLS))
}
pooled <- function(m) {
  m |> group_by(tool) |>
    summarise(TP = sum(TP), FP = sum(FP), FN = sum(FN), .groups = "drop") |>
    mutate(Precision = TP / (TP + FP), Recall = TP / (TP + FN),
           F1 = 2 * Precision * Recall / (Precision + Recall))
}
theme_fig <- theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), plot.tag = element_text(face = "bold", size = 13),
        plot.title = element_text(size = 10.5, face = "bold"))

# ---- (a) Petersen: per-library F1 + pooled Pr / Rc / F1 ----------------------
pm <- read_species("peterson_2026", "metrics.tsv")
pa <- bind_rows(
  pm |> transmute(tool, column = unname(PETERSEN_VIRUS[sample]), value = F1, group = "F1 per library"),
  pooled(pm) |> tidyr::pivot_longer(c(Precision, Recall, F1), names_to = "column") |>
    transmute(tool, column, value, group = "Pooled")
) |>
  mutate(column = factor(column, levels = c("Zika", "WNV", "SARS-CoV-2", "DENV",
                                            "Precision", "Recall", "F1")),
         group  = factor(group, levels = c("F1 per library", "Pooled")),
         tool   = factor(tool, levels = rev(TOOLS)))

score_tiles <- function(d, tag, title) {
  ggplot(d, aes(column, tool, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.2f", value),
                  colour = if_else(value > 0.6, "dark", "light")), size = 3.2) +
    scale_colour_manual(values = c(dark = "white", light = "grey15"), guide = "none") +
    scale_fill_gradient(low = "#f2f0f7", high = "#1b7837", limits = c(0, 1), name = NULL) +
    facet_grid(~ group, scales = "free_x", space = "free_x") +
    labs(tag = tag, title = title, x = NULL, y = NULL) +
    theme_fig + theme(legend.position = "none", strip.text = element_text(size = 9.5))
}
p_a <- score_tiles(pa, "a", "Four single-virus field libraries (Petersen 2026)")

# ---- (b) Gauthier: pooled Pr / Rc / F1 and the false-positive matrix -------------
patients <- gauthier_patients()
gm <- read_species("gauthier_2021", "metrics.tsv", patients)
gp <- pooled(gm) |> tidyr::pivot_longer(c(Precision, Recall, F1), names_to = "column") |>
  transmute(tool = factor(tool, levels = rev(TOOLS)),
            column = factor(column, levels = c("Precision", "Recall", "F1")),
            value, group = "Pooled")
p_b1 <- score_tiles(gp, "b", "16 SARS-CoV-2 clinical samples (Gauthier 2021)") +
  theme(plot.title.position = "plot",   # full title fits; margin clears the "b" tag
        plot.title = element_text(size = 10.5, face = "bold", margin = margin(l = 14, b = 5.5)))

fp <- read_species("gauthier_2021", "per_taxon.tsv", patients) |>
  filter(status == "FP") |>
  count(tool, taxon, name = "samples") |>
  mutate(taxon = if_else(taxon == "unknown", "\"unknown\" (above species)", taxon))
fp_order <- fp |> group_by(taxon) |>
  summarise(first = min(as.integer(tool)), n = sum(samples), .groups = "drop") |>
  arrange(first, desc(n)) |> pull(taxon)
# Pipelines as rows (same order as the score tiles), species as 45-degree columns.
fp <- fp |> mutate(taxon = factor(taxon, levels = fp_order),
                   tool = factor(tool, levels = rev(TOOLS)),
                   group = "Number of samples reported in")
p_b2 <- ggplot(fp, aes(taxon, tool, fill = tool)) +
  geom_tile(colour = "white", linewidth = 0.6, alpha = 0.85) +
  geom_text(aes(label = samples), size = 2.9, colour = "white") +
  scale_fill_manual(values = COLS, guide = "none") +
  scale_y_discrete(drop = FALSE) +
  facet_grid(~ group) +   # a strip row like the score tiles', so the two align
  labs(title = "False-positive species (Gauthier 2021)", x = NULL, y = NULL) +
  theme_fig + theme(axis.text.x = element_text(face = "italic", size = 8,
                                               angle = 45, hjust = 1, vjust = 1),
                    axis.text.y = element_text(size = 8.5),
                    strip.text = element_text(size = 9.5))

# ---- layout: (a) on top; (b) score tiles left, false-positive matrix right ---------
out <- "figures/frameworkA_detection.png"
png(out, width = 11, height = 7.2, units = "in", res = 150)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(c(0.42, 0.58), "npc"),
                                           widths = unit(c(0.4, 0.6), "npc"))))
print(p_a,  vp = viewport(layout.pos.row = 1, layout.pos.col = 1:2))
# Same row heights in both (b) panels, so their pipeline rows line up.
g_b1 <- ggplotGrob(p_b1); g_b2 <- ggplotGrob(p_b2)
stopifnot(length(g_b1$heights) == length(g_b2$heights))
g_b1$heights <- g_b2$heights <- grid::unit.pmax(g_b1$heights, g_b2$heights)
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1)); grid.draw(g_b1); popViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2)); grid.draw(g_b2); popViewport()
invisible(dev.off())

# pooled values, to check against Fig. 2
print(pooled(pm) |> mutate(across(c(Precision, Recall, F1), \(x) round(x, 2))))
print(pooled(gm) |> mutate(across(c(Precision, Recall, F1), \(x) round(x, 2))))
cat("wrote", out, "\n")
