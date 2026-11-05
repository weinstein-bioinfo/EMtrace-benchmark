#!/usr/bin/env Rscript
# Framework C — real-time replay, manuscript Fig. 5 (one column, four panels)
#
#   a  viral: target reads retrieved against reads processed (log), per library and
#      annotation chain (frameworkC_realtime_combined.R, panel a1)
#   b  viral: genome coverage of the reference-based consensus against reads processed;
#      legend gives reads aligned in the final assembly (frameworkC_realtime_combined.R)
#   c  16S: L1 error at genus (solid) and species (dashed), barcode79 and the
#      gut mock (mean +/- SD over four replicates)
#   d  16S: species precision (solid) and recall (dashed), same libraries
# The viral panels use reads processed on a log axis (the tick pools have no sequencing
# time; libraries span 44,503-464,000 reads); the 16S panels use cumulative reads as a
# percentage of the library. Viral libraries: the project notes.
# 16S scoring: near-synonyms collapsed; a taxon counts if it has >= 10 EM-expected
# reads and >= 0.01 % of them (`*_floor` columns of the ranks tables).
#
# Inputs (repo root):
#   inputs of rplots/frameworkC_realtime_combined.R (sourced for panels a, b)
#   metrics/framework_C/<library>.ranks.tsv     framework_C.barcode79.ranks.py --library
# Outputs:
#   figures/frameworkC_realtime_{ab,c,d}_*.png  one PNG per panel (a, b together)
#   figures/frameworkC_realtime.png            the assembled figure
#
# Run from the repo root:  Rscript rplots/frameworkC_realtime_figure.R

invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))   # keep the plus-minus, en dash and >= signs
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2); library(grid)
})

MET <- "metrics/framework_C"
FIG <- "figures"
W   <- 7.2        # inches, every panel
DPI <- 200

S16_COLS <- c(barcode79 = "#2166ac", gutmock = "black")  # blue and black, 16S panels
GUT_REPS <- paste0("gutmock_sup_r", 1:4)

fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)
pct_axis <- scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                               labels = \(x) paste0(x, "%"), expand = expansion(mult = 0.01))
theme_panel <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.tag = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 7.5), legend.title = element_text(size = 8),
        legend.key.width = unit(1.6, "lines"), legend.key.height = unit(0.75, "lines"),
        legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
        legend.box = "horizontal", legend.box.just = "top",
        plot.margin = margin(4, 12, 4, 4))
# Legend inside the panel at (x, y) in npc, anchored by its corner `just`.
inside <- function(x, y, just = c(1, 0)) theme(legend.position = "inside",
  legend.position.inside = c(x, y), legend.justification.inside = just)
X_LAB <- "Cumulative reads (% of library)"

# ---- viral panels a, b ------------------------------------------------------------
# Built by frameworkC_realtime_combined.R (also writes its own frameworkC_viral_* files).
# Here a and b sit side by side and share one legend below them: the same five lines,
# labelled with the reads aligned in b's final assembly, plus a's annotation chain.
viral <- new.env()
sys.source("rplots/frameworkC_realtime_combined.R", envir = viral)
no_legend <- theme(legend.position = "none", plot.margin = margin(4, 8, 2, 4))
p_a <- viral$p_a1 + labs(title = "Target reads retrieved") + no_legend
p_b <- viral$p_b + labs(title = "Genome coverage") + no_legend
p_leg <- suppressMessages(
  viral$p_a1 + scale_colour_manual(values = viral$LINES, labels = viral$cov_labels,
                                   name = "Reads aligned in the final assembly (% of library):") +
    guides(colour = guide_legend(ncol = 2, order = 1, title.position = "top"),
           linetype = guide_legend(order = 2)) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
          legend.key.height = unit(0.9, "lines"), legend.spacing.y = unit(2, "pt"),
          legend.margin = margin(0, 0, 0, 0)))
viral_legend <- gtable::gtable_filter(ggplotGrob(p_leg), "guide-box-bottom")
# Its row is as tall as the legend itself (+ a small gap), in units of H.
legend_rel <- (convertHeight(sum(viral_legend$heights), "in", valueOnly = TRUE) + 0.15) / 2.15

# ---- 16S ------------------------------------------------------------------------
read_ranks <- function(lib) read_tsv(file.path(MET, paste0(lib, ".ranks.tsv")),
                                     show_col_types = FALSE) |>
  mutate(pct = 100 * reads_cum / max(reads_cum), total = max(reads_cum))
COLS16 <- c(L1_genus = "L1_genus_coll", L1_species = "L1_species_coll",
            Pr = "Pr_species_floor", Rc = "Rc_species_floor")

b79 <- read_ranks("barcode79")
gut <- purrr::map_dfr(GUT_REPS, \(l) read_ranks(l) |> mutate(rep = l))
gut_sizes <- gut |> distinct(rep, total)

# Replicates have different batch boundaries in % terms: interpolate each onto one
# grid (from the latest first batch to 100 %), then take mean and SD.
grid_pct <- seq(ceiling(max(tapply(gut$pct, gut$rep, min))), 100, by = 0.5)
gut_grid <- gut |> group_by(rep) |>
  reframe(across(all_of(unname(COLS16)), \(v) approx(pct, v, xout = grid_pct)$y)) |>
  mutate(pct = grid_pct, .by = rep)
gut_sum <- gut_grid |> group_by(pct) |>
  summarise(across(all_of(unname(COLS16)), list(mean = mean, sd = sd)), .groups = "drop")

lab16 <- c(barcode79 = sprintf("barcode79 (%s reads)", fmt_n(max(b79$reads_cum))),
           gutmock   = sprintf("gut mock, mean ± SD of 4 replicates (%s–%s reads)",
                               fmt_n(min(gut_sizes$total)), fmt_n(max(gut_sizes$total))))
cols16 <- setNames(S16_COLS, lab16)

long16 <- function(metrics, names) {
  b <- b79 |> select(pct, all_of(unname(metrics))) |>
    pivot_longer(-pct, names_to = "col", values_to = "mean") |> mutate(sd = NA_real_,
                                                                       set = lab16[["barcode79"]])
  g <- gut_sum |> select(pct, ends_with("_mean"), ends_with("_sd")) |>
    pivot_longer(-pct, names_to = c("col", ".value"), names_pattern = "(.*)_(mean|sd)$") |>
    filter(col %in% metrics) |> mutate(set = lab16[["gutmock"]])
  bind_rows(b, g) |> mutate(series = factor(names[match(col, metrics)], levels = names),
                            set = factor(set, levels = lab16))
}

# L1 as two stacked sub-panels with their own y ranges, so each set's batch-to-batch
# movement is visible (on a shared 0-1 axis barcode79's early peak, +0.017, vanishes).
# Points are batches. Gut mock: replicates averaged by batch number (x = mean read %,
# band = +/- SD) over the batches all four have (1-18; they have 19-22), plus each
# replicate's final batch at 100 %, so every point is a four-replicate mean.
L1_RANKS <- c(L1_genus_coll = "Genus", L1_species_coll = "Species")
b79_L1 <- b79 |> select(batch, pct, all_of(names(L1_RANKS))) |>
  pivot_longer(all_of(names(L1_RANKS)), names_to = "col", values_to = "mean") |>
  mutate(sd = NA_real_, n = 1L)
n_common <- min(tapply(gut$batch, gut$rep, max)) - 1
gut_L1 <- gut |> group_by(rep) |>
  filter(batch <= n_common | batch == max(batch)) |>
  mutate(batch = if_else(batch == max(batch), n_common + 1L, batch)) |> ungroup() |>
  select(rep, batch, pct, all_of(names(L1_RANKS))) |>
  pivot_longer(all_of(names(L1_RANKS)), names_to = "col", values_to = "v") |>
  group_by(batch, col) |>
  summarise(pct = mean(pct), mean = mean(v), sd = sd(v), n = n(), .groups = "drop")
# Events, as in frameworkC_barcode79_stabilisation.R:
#   first correct call  first batch at which every truth genus and species is recovered
#                       and stays recovered (never, on the gut mock: four genera sit below
#                       the read threshold by design, so recall stays at 0.73-0.80)
#   stabilised          first batch after which genus and species L1 both stay within
#                       +/-L1_TOL of their final values; gut mock: the latest replicate
L1_TOL <- 0.01
stabilised_at <- function(v, tol = L1_TOL) {
  ok <- abs(v - v[length(v)]) <= tol
  min(which(rev(cumprod(rev(ok))) == 1))
}
first_call_at <- function(d) {
  ok <- d$Rc_genus_coll == 1 & d$Rc_species_coll == 1
  if (!any(ok)) return(NA_integer_)
  min(which(ok & rev(cumprod(rev(ok))) == 1))
}
stab_of <- function(d) max(stabilised_at(d$L1_genus_coll), stabilised_at(d$L1_species_coll))
ev_text <- function(what, b) sprintf("%s: batch %d\n%s reads, %.1f h seq., %.0f min compute",
                                     what, b$batch, fmt_n(b$reads_cum), b$seq_time_s / 3600,
                                     b$wall_s / 60)
b79_first <- first_call_at(b79); b79_stab <- stab_of(b79)
ev_b79 <- tibble(pct = b79$pct[c(b79_first, b79_stab)],
                 label = c(ev_text("First correct call", b79[b79_first, ]),
                           ev_text(sprintf("Stabilised (L1 ±%.2f)", L1_TOL), b79[b79_stab, ])))
gut_stab <- vapply(split(gut, gut$rep), stab_of, 1L)
gut_at <- gut |> filter(batch == max(gut_stab))
ev_gut <- tibble(
  pct = mean(gut_at$pct),
  label = sprintf(paste0("Stabilised (L1 ±%.2f) in all 4 replicates by batch %d ",
                         "(per replicate: %s)\n%s reads, ≤%.2f h seq., ≤%.0f min compute"),
                  L1_TOL, max(gut_stab), paste(gut_stab, collapse = ", "),
                  fmt_n(unique(gut_at$reads_cum)), max(gut_at$seq_time_s) / 3600,
                  max(gut_at$wall_s) / 60))

# One L1 sub-panel. Headroom above the data holds the library label (row 1) and one
# event label per row below it, each starting just right of its dashed line.
L1_sub <- function(d, colour, label, top, events) {
  d <- d |> mutate(series = factor(L1_RANKS[col], levels = L1_RANKS))
  lo <- min(d$mean - coalesce(d$sd, 0)); hi <- max(d$mean + coalesce(d$sd, 0))
  span <- hi - lo; rows <- 1 + nrow(events)
  y_top <- hi + span * (0.15 + 0.5 * rows)
  row_y <- \(r) y_top - span * (0.25 + 0.5 * (r - 1))
  events <- events |> mutate(y = row_y(row_number() + 1))
  ggplot(d, aes(pct, mean, linetype = series, shape = series)) +
    geom_vline(data = events, aes(xintercept = pct), linetype = "dashed",
               colour = "grey45", linewidth = 0.4) +
    geom_ribbon(data = filter(d, !is.na(sd)), aes(ymin = mean - sd, ymax = mean + sd),
                fill = colour, alpha = 0.18, colour = NA, show.legend = FALSE) +
    geom_line(colour = colour, linewidth = 0.6) +
    geom_point(colour = colour, size = 1.3) +
    annotate("label", x = 1, y = row_y(1), label = label, hjust = 0, vjust = 0.5,
             colour = colour, size = 2.7, fill = "white", linewidth = 0,
             label.padding = unit(0.08, "lines")) +
    geom_label(data = events, aes(x = pct + 0.6, y = y, label = label), inherit.aes = FALSE,
               hjust = 0, vjust = 0.5, size = 2.3, colour = "grey25", lineheight = 0.95,
               fill = "white", linewidth = 0, label.padding = unit(0.08, "lines")) +
    scale_linetype_manual(values = c(Genus = "solid", Species = "22"), name = NULL) +
    scale_shape_manual(values = c(Genus = 16, Species = 1), name = NULL) +
    scale_y_continuous(limits = c(lo - 0.1 * span, y_top), expand = c(0, 0)) +
    pct_axis + labs(x = X_LAB, y = "L1 error",
                    title = if (top) "16S mock communities: L1 error of the abundance profile") +
    theme_panel +
    (if (top) inside(0.99, 0.98, c(1, 1)) +
       theme(axis.title.x = element_blank(), axis.text.x = element_blank())
     else theme(legend.position = "none"))
}
p_d1 <- L1_sub(b79_L1, S16_COLS[["barcode79"]], lab16[["barcode79"]], top = TRUE, ev_b79)
p_d2 <- L1_sub(gut_L1, S16_COLS[["gutmock"]],
               sub("mean ± SD of 4 replicates", "mean ± SD of 4 replicates by batch",
                   lab16[["gutmock"]]), top = FALSE, ev_gut)

e16 <- long16(COLS16[c("Pr", "Rc")], c("Precision", "Recall"))
p_e <- ggplot(e16, aes(pct, mean, colour = set, linetype = series)) +
  geom_ribbon(data = filter(e16, !is.na(sd)), aes(ymin = mean - sd, ymax = mean + sd, fill = set),
              colour = NA, alpha = 0.18, show.legend = FALSE) +
  geom_vline(data = bind_rows(mutate(ev_b79, set = lab16[["barcode79"]]),
                              mutate(ev_gut, set = lab16[["gutmock"]])) |>
               mutate(set = factor(set, levels = lab16)),
             aes(xintercept = pct, colour = set), linetype = "dashed", linewidth = 0.35,
             show.legend = FALSE) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = cols16, name = NULL,
                      labels = c("barcode79", "gut mock (mean \u00b1 SD, n = 4)")) +
  scale_fill_manual(values = cols16, guide = "none") +
  scale_linetype_manual(values = c(Precision = "solid", Recall = "22"), name = "Species rank") +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  pct_axis + scale_y_continuous(oob = scales::oob_squish, limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
  labs(tag = "e", title = "16S mock communities: species precision and recall (≥10 reads, ≥0.01%)",
       x = X_LAB, y = "Precision / recall") + theme_panel + inside(0.01, 0.04, c(0, 0)) +
  theme(legend.title = element_blank())

# ---- write one PNG per panel, then the assembled figure ---------------------------
# Row 1: a | b, equal heights; row 2: their shared legend; then the 16S plots c (two
# stacked L1 plots) and d, left-aligned with each other. Heights in units of H.
H <- 2.15
grob_tag <- \(p, t) ggplotGrob(p + labs(tag = t))
ga <- grob_tag(p_a, "a"); gb <- grob_tag(p_b, "b")
maxh <- unit.pmax(ga$heights, gb$heights)
ga$heights <- maxh; gb$heights <- maxh
viral_row <- gtable::gtable_row("viral", list(ga, gb), height = unit(1, "null"),
                                widths = unit(c(1, 1), "null"))
g16 <- list(grob_tag(p_d1, "c"), ggplotGrob(p_d2), grob_tag(p_e, "d"))
maxw <- do.call(unit.pmax, lapply(g16, \(g) g$widths))
g16 <- lapply(g16, \(g) { g$widths <- maxw; g })

draw_stack <- function(gs, hs, file) {
  ragg::agg_png(file, width = W, height = H * sum(hs), units = "in", res = DPI,
                background = "white")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(length(gs), 1, heights = unit(hs, "null"))))
  for (i in seq_along(gs)) {
    pushViewport(viewport(layout.pos.row = i, layout.pos.col = 1))
    grid.draw(gs[[i]]); popViewport()
  }
  invisible(dev.off())
}
# Stale panel files of this script only (other frameworkC_realtime_* figures are kept).
old <- Sys.glob(file.path(FIG, paste0("frameworkC_realtime_",
                                      c("?_viral_*", "?_16S_*", "ab_*"), ".png")))
invisible(file.remove(old))
pieces <- list(ab_viral = list(list(viral_row, viral_legend), c(1.2, legend_rel)),
               c_16S_L1 = list(g16[1:2], c(0.95, 0.9)),
               d_16S_precision_recall = list(g16[3], 1))
for (k in names(pieces))
  draw_stack(pieces[[k]][[1]], pieces[[k]][[2]],
             file.path(FIG, sprintf("frameworkC_realtime_%s.png", k)))
out <- file.path(FIG, "frameworkC_realtime.png")
draw_stack(unlist(lapply(pieces, `[[`, 1), recursive = FALSE),
           unlist(lapply(pieces, `[[`, 2)), out)

# Values for the caption.
print(gut_sizes)
cat("wrote", out, "and", length(pieces), "panel PNGs\n")
