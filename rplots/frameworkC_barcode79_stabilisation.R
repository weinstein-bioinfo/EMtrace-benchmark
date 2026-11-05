#!/usr/bin/env Rscript
# Framework C — barcode79 (ZymoBIOMICS 16S) stabilisation curves
#
# the project notes, required figure 2, for the profiling library: how the answer
# evolves batch by batch, with the two events that define "when the answer
# arrived" marked on the x-axis.
#
# Two panels over one shared x (cumulative reads; sequencing time on the top
# axis):
#   top    profile error — L1 against the pct_16S_only truth, genus and species
#   bottom taxon set     — precision and recall, genus and species, plus
#                          species precision on the 16S benchmark's 1 % tail-
#                          trimmed cell
#
# They are stacked rather than overlaid because the point of metric 4 is that
# the two do NOT move together: L1 settles at the first batch while precision
# decays all run as the false-positive tail accumulates. A shared x makes that
# readable; a shared y would not, the scales being unrelated.
#
# THREE SERIES COINCIDE AT 1.00 in the bottom panel — recall at both ranks and
# the tail-trimmed species precision — and that coincidence is the finding, so
# they are drawn at their true values rather than nudged apart. The tail-trim
# line is drawn first and the two recall lines dashed over it, so all three stay
# visible; the annotation names them.
#
# The tail-trimmed series answers a question the other two cannot: the
# Framework C precision (0.40 at end of run) and the 16S benchmark's
# precision (1.00) are the same run scored two ways. They differ because the
# whole false-positive tail lives inside the smallest 1 % of reads. Its scoring
# path also differs — the 16S trim drops abundance-table ROWS and collapses to a
# rank afterwards, where Module C sums to the rank first and thresholds after
# (see framework_C.barcode79.ranks.py) — so it is a third colour, not a fourth
# rank.
#
# Both events are drawn as vertical rules:
#   first correct call   all 8 truth genera and species recovered (recall 1.0)
#   stabilisation        first batch after which L1 stays within +/-0.01 of its
#                        final value. The roadmap's +/-0.02 band is already met
#                        at batch 1, which would put both rules on the same
#                        line and say nothing, so the tighter band is drawn and
#                        the roadmap one is stated in the caption.
#
# Species are near-synonym collapsed (data/truth/zymobiomics_16S_collapse.tsv),
# the basis the 16S benchmark's headline numbers use; strict species L1 sits at
# 0.58 throughout because the B. subtilis group splits across four labels, and
# the flat offset would only compress the panel.
#
# Input (repo root):
#   metrics/framework_C/barcode79.ranks.tsv
#     written by jobs/framework_C/framework_C.barcode79.ranks.py
# Output:
#   figures/frameworkC_barcode79_stabilisation.png
#
# Run from the repo root:
#   Rscript rplots/frameworkC_barcode79_stabilisation.R

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(ggplot2)
})

in_file  <- "metrics/framework_C/barcode79.ranks.tsv"
out_dir  <- "figures"
out_file <- "frameworkC_barcode79_stabilisation.png"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

L1_TOL <- 0.01   # stabilisation band drawn (roadmap band is 0.02)

# One entry per line in the figure: display name, colour, linetype. Drawing
# order is this order, so the tail-trim line goes down before the dashed recall
# lines that sit on top of it.
SERIES <- tribble(
  ~series,                              ~colour,   ~linetype,
  "Genus (L1, precision)",           "#1b9e77", "solid",
  "Species (L1, precision)",         "#7570b3", "solid",
  "Species precision, 1 % tail trim",   "#d95f02", "solid",
  "Genus recall",                       "#1b9e77", "42",
  "Species recall",                     "#7570b3", "42",
)

dat <- read_tsv(in_file, show_col_types = FALSE)

# --- the two events ----------------------------------------------------------
# First correct call: the roadmap's rank-1 rule never fires on this library (the
# top genus is Bacillus at every batch, the truth's most abundant is
# Limosilactobacillus, 17.4 % against 18.4 %), so the profiling equivalent is
# used — the first batch at which the truth set is complete and stays complete.
first_call_batch <- with(dat, {
  ok <- Rc_genus_coll == 1 & Rc_species_coll == 1
  min(which(ok & rev(cumprod(rev(ok))) == 1))
})

# First batch after which a curve stays within tol of its final value
# (ngs_benchmark.realtime.metrics.stabilised_at).
stabilised_at <- function(v, tol) {
  final <- v[length(v)]
  min(which(rev(cumprod(rev(abs(v - final) <= tol))) == 1))
}
stab_batch <- max(stabilised_at(dat$L1_genus_coll,   L1_TOL),
                  stabilised_at(dat$L1_species_coll, L1_TOL))

events <- tibble(
  batch = c(first_call_batch, stab_batch),
  label = c("First correct call", sprintf("Stabilised (L1 ±%.2f)", L1_TOL))
) |>
  left_join(dat |> select(batch, reads_cum, seq_time_s, wall_s), by = "batch") |>
  mutate(label = sprintf("%s — batch %d\n%s reads, %.1f h seq., %.0f min compute",
                         label, batch, format(reads_cum, big.mark = ","),
                         seq_time_s / 3600, wall_s / 60))

# --- long form, one row per (batch, panel, series) ---------------------------
P_TOP <- "Profile error — L1 vs truth"
P_BOT <- "Precision / Recall (≥ 10 reads)"

long <- bind_rows(
  dat |>
    select(batch, reads_cum,
           `Genus (L1, precision)`   = L1_genus_coll,
           `Species (L1, precision)` = L1_species_coll) |>
    pivot_longer(-c(batch, reads_cum), names_to = "series", values_to = "value") |>
    mutate(panel = P_TOP),
  dat |>
    select(batch, reads_cum,
           `Genus (L1, precision)`      = Pr_genus_coll,
           `Species (L1, precision)`    = Pr_species_coll,
           `Species precision, 1 % tail trim` = Pr_species_tail,
           `Genus recall`                     = Rc_genus_coll,
           `Species recall`                   = Rc_species_coll) |>
    pivot_longer(-c(batch, reads_cum), names_to = "series", values_to = "value") |>
    mutate(panel = P_BOT)
) |>
  mutate(panel  = factor(panel, levels = c(P_TOP, P_BOT)),
         series = factor(series, levels = SERIES$series))

pal <- setNames(SERIES$colour,   SERIES$series)
lty <- setNames(SERIES$linetype, SERIES$series)

# Right-edge labels, so a reader need not walk back to the legend. y is the
# series' own final value except where three of them coincide at 1.00.
ends <- long |>
  filter(batch == max(batch)) |>
  mutate(text = case_when(
    panel == P_TOP & series == "Genus (L1, precision)"   ~ "Genus L1  0.261",
    panel == P_TOP & series == "Species (L1, precision)" ~ "Species L1  0.257",
    series == "Genus (L1, precision)"                    ~ "Genus precision  0.44",
    series == "Species (L1, precision)"                  ~ "Species precision  0.40",
    series == "Species precision, 1 % tail trim" ~
      "Recall, both ranks  1.00\nSpecies precision, 1 % tail trim  1.00",
    TRUE ~ NA_character_
  ), y = if_else(series == "Species precision, 1 % tail trim", 1.0, value)) |>
  filter(!is.na(text))

read_breaks <- seq(0, 100000, by = 20000)
seq_at <- approx(dat$reads_cum, dat$seq_time_s / 3600, xout = read_breaks, rule = 2)$y

p <- ggplot(long, aes(reads_cum, value, colour = series, linetype = series)) +
  geom_vline(data = events, aes(xintercept = reads_cum),
             linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.1) +
  geom_text(data = events |> mutate(panel = factor(P_BOT, levels = levels(long$panel))),
            aes(x = reads_cum, y = 0.44, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = -0.3, angle = 90, size = 2.5,
            colour = "grey25", lineheight = 0.95) +
  geom_text(data = ends, aes(x = reads_cum, y = y, label = text, colour = series),
            inherit.aes = FALSE, hjust = 0, nudge_x = 2200, size = 2.6,
            lineheight = 1.05, show.legend = FALSE) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_linetype_manual(values = lty, name = NULL) +
  scale_x_continuous(
    breaks = read_breaks, labels = scales::label_comma(),
    limits = c(0, 132000),
    sec.axis = sec_axis(~ ., breaks = read_breaks,
                        labels = sprintf("%.0f h", seq_at),
                        name = "Sequencing time (from read start_time headers)")
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(
    title = "barcode79 (ZymoBIOMICS 16S) — what the real-time answer looks like, batch by batch",
    subtitle = paste(
      "26 MinKNOW batches, 100,128 reads, 71.8 h of sequencing replayed in 123 min of compute. Species near-synonym collapsed.",
      "Recall at both ranks and tail-trimmed species precision all sit at 1.00 every batch, so those three lines coincide.",
      sep = "\n"),
    x = "Cumulative reads", y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    strip.placement    = "outside",
    strip.text.y.left  = element_text(angle = 90, face = "bold", size = 9),
    panel.grid.minor   = element_blank(),
    plot.subtitle      = element_text(size = 8, colour = "grey35"),
    axis.title.x.top   = element_text(size = 8, colour = "grey35", margin = margin(b = 4)),
    legend.position    = "bottom",
    legend.key.width   = unit(2.0, "lines"),
    legend.text        = element_text(size = 8),
    plot.margin        = margin(6, 12, 6, 6)
  )

ggsave(file.path(out_dir, out_file), p, width = 10.5, height = 6.6, dpi = 300)
message("Wrote ", file.path(out_dir, out_file),
        "  (first call batch ", first_call_batch, ", stabilised batch ", stab_batch, ")")
