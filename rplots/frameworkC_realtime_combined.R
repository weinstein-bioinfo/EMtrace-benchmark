#!/usr/bin/env Rscript
# Framework C — manuscript figures across the four reported libraries
# (design: the project notes)
#
#   a1  main target reads retrieved against reads processed, one line per library and
#       annotation chain: Orbivirus_kevo (genus Orbivirus), TBEV_FBC92293 (species TBEV),
#       ticks_58-37 and ticks_35-33 (genus Orthonairovirus, CCHFV)
#   a2  all panel taxa of the two tick pools, faceted by pool
#   b   genome coverage against reads processed: KEVOOV, TBEV, CCHFV (both pools),
#       JMTV (58-37); one line per virus, from the RefSeq-only replay's reads
#
# Reads are each batch's cumulative `readn` at the taxon's rank (the count the panel
# scoring uses). x is cumulative reads processed, log scale: the tick pools have no
# sequencing time, and the libraries differ ten-fold in size.
#
# Inputs (repo root):
#   metrics/framework_C/realtime_reads.curve.tsv   framework_C.curve.py --reads
#   metrics/framework_C/coverage_over_time.tsv     framework_C.coverage_over_time.bash
#   metrics/framework_C/coverage_mapped_reads.tsv  framework_C.coverage_mapped_reads.bash
# Outputs:
#   figures/frameworkC_viral_a_reads.png
#   figures/frameworkC_viral_a2_tick_panels.png
#   figures/frameworkC_viral_b_coverage.png
#
# Run from the repo root:  Rscript rplots/frameworkC_realtime_combined.R

invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
suppressPackageStartupMessages({ library(readr); library(dplyr); library(ggplot2) })

MET <- "metrics/framework_C"
FIG <- "figures"
W   <- 7.2
DPI <- 200
THRESHOLD <- 10

# Figure labels (the project notes, Labels): virus, then the sample in the reader's
# terms. Tick pools keep the pool IDs of Ergunay 2023; the SRA runs are in the caption.
# One colour per line of a1 and b; JMTV appears in b only.
LAB_KEVO  <- "KEVOOV, isolate fi663 (cell culture)"
LAB_TBEV  <- "TBEV, 4 pooled isolates"
LAB_JMTV  <- "JMTV, tick pool 58-37"
LINES <- setNames(c("#1b9e77", "#7570b3", "#66a61e", "#d95f02", "#e7298a"),
                  c(LAB_KEVO, LAB_TBEV, LAB_JMTV, "CCHFV, tick pool 58-37",
                    "CCHFV, tick pool 35-33"))
CHAINS <- c("+ DIAMOND" = "solid", "RefSeq only" = "22")
TARGET <- c(Orbivirus_kevo = LAB_KEVO, TBEV_FBC92293 = LAB_TBEV,
            `ticks_58-37` = "CCHFV, tick pool 58-37", `ticks_35-33` = "CCHFV, tick pool 35-33")

# Reads aligned to each figure-b reference in the final assembly, and their share of the
# library's raw reads (framework_C.coverage_mapped_reads.bash). They label b and set the
# legend order of every panel: descending, independent of the annotation chain.
mapped <- read_tsv(file.path(MET, "coverage_mapped_reads.tsv"), show_col_types = FALSE) |>
  mutate(lib = sub("\\..*", "", library),
         line = if_else(startsWith(segment, "MW5611"), LAB_JMTV, TARGET[lib])) |>
  group_by(line) |>
  summarise(mapped = sum(mapped_reads), total = first(library_reads)) |>
  arrange(desc(mapped)) |>
  mutate(label = sprintf("%s: %s (%s%%)", line, scales::comma(mapped),
                         trimws(formatC(100 * mapped / total, format = "fg", digits = 2))))
LINES <- LINES[mapped$line]
cov_labels <- with(mapped, setNames(label, line))

reads_axis <- scale_x_log10(labels = scales::label_comma(),
                            breaks = c(4e3, 1e4, 3e4, 1e5, 3e5), name = "Reads processed")
theme_fig <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.tag = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 7.5), legend.title = element_text(size = 8),
        legend.key.width = unit(1.8, "lines"),
        strip.text = element_text(size = 9, face = "bold"),
        plot.margin = margin(4, 12, 4, 4))
threshold_line <- geom_hline(yintercept = THRESHOLD, linetype = "dotted", colour = "grey40")
save <- function(p, file, h) {
  ragg::agg_png(file.path(FIG, file), width = W, height = h, units = "in", res = DPI,
                background = "white")
  print(p); invisible(dev.off()); cat("wrote", file.path(FIG, file), "\n")
}

rc <- read_tsv(file.path(MET, "realtime_reads.curve.tsv"), show_col_types = FALSE) |>
  mutate(chain = factor(chain, names(CHAINS)))

# ---- a1: main target of each library, plus JMTV (the second PCR target of 58-37) --------
a1 <- rc |>
  filter(main | (library == "ticks_58-37" & taxon == "Jingmenvirus")) |>
  mutate(line = factor(if_else(taxon == "Jingmenvirus", LAB_JMTV, TARGET[library]),
                       names(LINES))) |>
  filter(readn > 0)
p_a1 <- ggplot(a1, aes(reads_cum, readn, colour = line, linetype = chain)) +
  threshold_line +
  annotate("text", x = 4.6e5, y = THRESHOLD * 1.25, label = "10-read threshold", hjust = 1,
           vjust = 0, size = 2.5, colour = "grey30") +
  geom_line(linewidth = 0.6) + geom_point(data = filter(a1, batch == 1), size = 1.6) +
  reads_axis +
  scale_y_log10(labels = scales::label_comma(), name = "Target reads retrieved") +
  scale_colour_manual(values = LINES, name = NULL, drop = TRUE) +
  scale_linetype_manual(values = CHAINS, name = "Annotation:") +
  labs(tag = "a", title = "Target reads against reads processed") + theme_fig

# ---- a2: panel taxa of the tick pools --------------------------------------------
# The same taxon can be a PCR target in one pool and a PCR negative in the other, so the
# evidence is written at the end of each line, not in the colour legend.
EVIDENCE <- c(pcr_pos = "PCR+", pcr_neg = "PCR−", reported = "Nanopore+")
a2 <- rc |> filter(startsWith(library, "ticks_")) |>
  mutate(pool = sub("ticks_", "Tick pool ", library),
         taxon_lab = sprintf("%s (%s)", taxon, sub(",.*", "", source)))
ends <- a2 |> group_by(pool, taxon) |> filter(batch == max(batch)) |>
  summarise(reads_cum = first(reads_cum), readn = max(readn),
            evidence = EVIDENCE[first(evidence)], .groups = "drop")
taxon_labs <- with(distinct(a2, taxon, taxon_lab), setNames(taxon_lab, taxon))
p_a2 <- ggplot(a2, aes(reads_cum, readn, colour = taxon)) +
  threshold_line +
  geom_line(aes(linetype = chain), linewidth = 0.6) +
  geom_text(data = ends, aes(label = evidence), hjust = -0.15, size = 2.5, show.legend = FALSE) +
  facet_wrap(~ pool, nrow = 1, scales = "free_x") +
  scale_x_log10(labels = scales::label_comma(), name = "Reads processed",
                expand = expansion(mult = c(0.03, 0.28))) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1, base = 10),
                     breaks = c(0, 1, 10, 100, 1000, 1e4), labels = scales::label_comma(),
                     name = "Reads retrieved") +
  scale_colour_brewer(palette = "Set1", name = "Panel taxon:", labels = taxon_labs,
                      guide = guide_legend(nrow = 2)) +
  scale_linetype_manual(values = CHAINS, name = "Annotation:") +
  labs(tag = "a2", title = "Tick pools: every panel taxon (Ergunay 2023)") + theme_fig +
  theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left",
        legend.title.position = "left", legend.title = element_text(size = 8, vjust = 1), legend.margin = margin(0, 0, 0, 0))

save(p_a1, "frameworkC_viral_a_reads.png", 3.0)
save(p_a2, "frameworkC_viral_a2_tick_panels.png", 3.6)

# ---- b: genome coverage ------------------------------------------------------------
COV_LIBS <- c("Orbivirus_kevo.refseq_abv1g_k11w5", "TBEV_FBC92293.refseq_abv1g",
              "ticks_35-33.refseq_abv1g_k11w5", "ticks_58-37.refseq_abv1g_k11w5")
cov <- read_tsv(file.path(MET, "coverage_over_time.tsv"), show_col_types = FALSE) |>
  filter(library %in% COV_LIBS) |>
  mutate(lib = sub("\\..*", "", library),
         virus = if_else(startsWith(segment, "MW5611"), "JMTV", NA_character_),
         line = case_when(!is.na(virus) ~ LAB_JMTV, .default = TARGET[lib])) |>
  group_by(line, batch, polished, reads_cum) |>
  summarise(segments = n(), identity = weighted.mean(identity, covered_bp, na.rm = TRUE),
            ref_len = sum(ref_len), covered_bp = sum(covered_bp), .groups = "drop") |>
  mutate(coverage = 100 * covered_bp / ref_len, line = factor(line, names(LINES)))
missing <- setdiff(names(LINES), unique(as.character(cov$line)))
if (length(missing)) cat("no coverage rows for:", paste(missing, collapse = ", "), "\n")

print(mapped)

curve <- filter(cov, polished == 0)
p_b <- ggplot(curve, aes(reads_cum, coverage, colour = line)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.1) +
  geom_point(data = filter(cov, polished == 1), shape = 18, size = 3) +
  reads_axis +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     labels = \(x) paste0(x, "%"), name = "Genome coverage") +
  scale_colour_manual(values = LINES, labels = cov_labels, drop = TRUE,
                      name = "Final assembly: reads aligned (% of library)") +
  labs(tag = "b", title = "Genome coverage against reads processed") + theme_fig
save(p_b, "frameworkC_viral_b_coverage.png", 3.0)

# ---- multipanel figure: a (a1), b (b) stacked, panel areas aligned ---------------------
# a2 (tick-pool panel taxa) is left out: its CCHFV and JMTV lines repeat a1, and the panel
# outcomes are reported in the text; it stays a standalone (supplementary) figure.
panels <- list(p_a1 + labs(tag = "a"), p_b + labs(tag = "b"))
heights <- c(3.0, 3.0)
grobs <- lapply(panels, ggplotGrob)
maxw <- do.call(grid::unit.pmax, lapply(grobs, \(g) g$widths))
grobs <- lapply(grobs, \(g) { g$widths <- maxw; g })
draw_stack <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    length(grobs), 1, heights = grid::unit(heights, "null"))))
  for (i in seq_along(grobs)) {
    grid::pushViewport(grid::viewport(layout.pos.row = i, layout.pos.col = 1))
    grid::grid.draw(grobs[[i]]); grid::popViewport()
  }
}
OUT_ALL <- file.path(FIG, "frameworkC_viral")
ragg::agg_png(paste0(OUT_ALL, ".png"), width = W, height = sum(heights), units = "in",
              res = DPI, background = "white")
draw_stack(); invisible(dev.off())
grDevices::cairo_pdf(paste0(OUT_ALL, ".pdf"), width = W, height = sum(heights))
draw_stack(); invisible(dev.off())
cat("wrote", paste0(OUT_ALL, c(".png", ".pdf")), "\n")

# Values for the captions.
print(rc |> group_by(run, taxon) |> slice_max(batch, n = 1) |>
        select(run, taxon, batches = batch, reads_cum, readn), n = Inf)
print(cov |> group_by(line, polished) |> filter(batch %in% range(batch)) |>
        select(line, batch, polished, reads_cum, segments, coverage, identity), n = Inf)
