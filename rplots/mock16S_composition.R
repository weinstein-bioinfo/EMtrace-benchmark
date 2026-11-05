#!/usr/bin/env Rscript
# Framework 16S — ZymoBIOMICS mock community composition
#
# The headline manuscript figure: one stacked bar per pipeline configuration,
# `Expected` leftmost as the reference, species/genera as stacked fills, and
# everything outside the truth set pooled into `Other`. Reproduces the layout
# Martyna Krupinska used for the LazypipeNP comparison.
#
# Two panels are written, one per rank. They are separate files rather than
# facets because Spaghetti reports at genus level by design and has no species
# bar, so a shared facet would carry a hole.
#
# Input (repo root):
#   metrics/framework_16S/composition_long.tsv   pipeline, rank, taxon, rel_abundance
# Outputs:
#   figures/mock16S_composition_species.png
#   figures/mock16S_composition_genus.png
#
# Run from the repo root:
#   Rscript rplots/mock16S_composition.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(forcats); library(ggplot2)
})

in_file <- "metrics/framework_16S/composition_long.tsv"
out_dir <- "figures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read_tsv(in_file, show_col_types = FALSE)

# composition_long.tsv now carries every library. The figure is one library at a
# time, so select it here rather than keeping a file per library.
library_id <- Sys.getenv("LIBRARY", "barcode79")
if ("library" %in% names(dat)) {
  dat <- dat %>% filter(library == library_id)
  if (nrow(dat) == 0) stop("no rows for library '", library_id, "' in ", in_file)
}

# The figure's bars, in display order: composition_long.tsv also carries every
# sweep cell, which belongs in the metrics tables, not in the figure. Names are
# manifest labels; values are the axis labels. Labels missing from the data are
# dropped (Spaghetti on the species panel), so the list can stay fixed.
#
# Run v3 cells (all hits, 17-18 Sep 2026). Cells are named by reference only: every
# database reports through its own taxonomy, so there is no SILVA-reference +
# NCBI-taxonomy bar. chain1 is omitted because it equals rrnDB on this mock, and
# the rrnDB.dedup `--nospnorep` bar because it equals rrnDB.dedup to 0.4 points.
# rrnDB.dedup -N 25 is the best NCBI-taxonomy cell; there is no tuned bar yet.
# Bars are EM-adjusted read counts (`readn_em`) where the table has them, else `readn`
# (Emu, Spaghetti); the scorer picks the weight.
bars <- c(
  "Expected"                 = "Expected",
  "Emu"                      = "Emu",
  "Spaghetti"                = "Spaghetti",
  "EMtrace.silva"               = "EMtrace SILVA",
  "EMtrace.rrnDB"               = "EMtrace rrnDB",
  "EMtrace.rrnDB.dedup.N5"      = "EMtrace rrnDB.dedup\n-N 5",
  "EMtrace.rrnDB.dedup"         = "EMtrace rrnDB.dedup",
  "EMtrace.rrnDB.dedup.N25"     = "EMtrace rrnDB.dedup\n-N 25"
)
dat <- dat |>
  filter(pipeline %in% names(bars)) |>
  mutate(pipeline = factor(bars[pipeline], levels = unname(bars)))

# Expected is always the reference bar, Other is always the last fill. Everything
# between keeps the order it appears in, so the eight mock organisms sit in a
# stable position across every bar.
taxon_levels <- function(d) {
  expected <- d |> filter(pipeline == "Expected") |> arrange(desc(rel_abundance)) |> pull(taxon)
  rest <- setdiff(unique(d$taxon), c(expected, "Other"))
  c(expected, sort(rest), "Other")
}

# Colour-blind-safe qualitative set; Other is deliberately neutral grey so the
# eye reads it as "unassigned mass" rather than as another organism.
fill_palette <- function(levels) {
  base <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
            "#937860", "#DA8BC3", "#17A398", "#CCB974", "#64B5CD",
            "#B07AA1", "#76B7B2", "#F28E2B", "#59A14F", "#EDC948")
  named <- setNames(rep(base, length.out = length(levels)), levels)
  if ("Other" %in% levels) named["Other"] <- "#BDBDBD"
  named
}

plot_rank <- function(rank_name, out_file, width, height) {
  d <- dat |> filter(rank == rank_name)
  if (nrow(d) == 0) {
    message("No rows for rank '", rank_name, "' — skipping ", out_file)
    return(invisible(NULL))
  }
  levels_t <- taxon_levels(d)
  d <- d |>
    mutate(
      taxon    = factor(taxon, levels = levels_t),
      pipeline = droplevels(pipeline)
    )

  # Label only the slices with room for text; below ~3 % the label collides with
  # its neighbours and the reader gets the value from the table instead.
  labels <- d |> filter(rel_abundance >= 3)

  p <- ggplot(d, aes(x = pipeline, y = rel_abundance, fill = taxon)) +
    geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
    geom_text(
      data = labels,
      aes(label = sprintf("%.1f", rel_abundance)),
      position = position_stack(vjust = 0.5),
      size = 2.7, colour = "white", fontface = "bold"
    ) +
    scale_fill_manual(values = fill_palette(levels_t), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = paste0("ZymoBIOMICS mock community — ", rank_name, "-level composition"),
      subtitle = "Relative abundance (%). Expected is the 16S-only basis; taxa outside the truth set are pooled into Other.",
      x = NULL, y = "Relative abundance (%)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x     = element_text(angle = 30, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.text     = element_text(face = "italic"),
      legend.key.size = unit(0.8, "lines"),
      plot.subtitle   = element_text(size = 8, colour = "grey35")
    )

  ggsave(file.path(out_dir, out_file), p, width = width, height = height, dpi = 300)
  message("Wrote ", file.path(out_dir, out_file))
}

n_bars <- length(unique(dat$pipeline))
# barcode79 keeps the established filenames, so existing figure references
# stay valid; any other library is suffixed.
stem <- if (library_id == "barcode79") "mock16S_composition" else
        paste0("mock16S_composition_", library_id)
plot_rank("species", paste0(stem, "_species.png"), width = max(7, 1.1 * n_bars + 4), height = 5.5)
plot_rank("genus",   paste0(stem, "_genus.png"),   width = max(7, 1.1 * n_bars + 4), height = 5.5)
