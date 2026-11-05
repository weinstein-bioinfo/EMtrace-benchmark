# Generate frozen R golden output for the simmet R-parity integration test.
#
# Reimplements bmark_stats() logic directly using the "division" field,
# bypassing the local NGSlib.R which has a bug in the output table loop
# (line 205: uses stale `taxgroup` variable instead of `taxgroups[i]`).
#
# Run once from the repo root (EMtrace-benchmark/):
#   Rscript ngs-benchmark/tests/scripts/generate_r_golden_simmet.R
#
# Commit the resulting tests/data/golden_r_simmet.tsv.
#
# Dependencies (install once in R):
#   install.packages("reshape")

library(reshape)

pred_file <- "Legacy/data/LazypipeX.simmet.vi_nt2/abund_table.tsv"
true_file <- "Legacy/data/2025_10_14.simmet.abund_table.tsv"
out_tsv   <- "ngs-benchmark/tests/data/golden_r_simmet.tsv"

taxranks  <- c("species", "genus", "family", "division")
taxgroups <- c("Viruses", "Bacteria")
tail_vi   <- 0
tail_ba   <- 0

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

pred <- read.table(pred_file, sep = "\t", header = TRUE,
                   comment.char = "", stringsAsFactors = FALSE, quote = "")
true <- read.table(true_file, sep = "\t", header = TRUE,
                   comment.char = "", stringsAsFactors = FALSE, quote = "")

pred <- pred[order(pred$readn, decreasing = TRUE), ]
true <- true[order(true$readn, decreasing = TRUE), ]

# ---------------------------------------------------------------------------
# NA / empty taxon labels -> "unknown"
# ---------------------------------------------------------------------------

for (rank in taxranks) {
  pred[grep("^$", pred[[rank]]), rank] <- "unknown"
  pred[is.na(pred[[rank]]), rank]      <- "unknown"
  true[grep("^$", true[[rank]]), rank] <- "unknown"
  true[is.na(true[[rank]]), rank]      <- "unknown"
}

# ---------------------------------------------------------------------------
# Filter phages
# ---------------------------------------------------------------------------

if ("bphage" %in% colnames(pred)) pred <- pred[pred$bphage != "yes", ]
if ("bphage" %in% colnames(true)) true <- true[true$bphage != "yes", ]

# ---------------------------------------------------------------------------
# Select only needed columns, then group-by + sum readn
# ---------------------------------------------------------------------------

pred <- pred[, c("readn", taxranks)]
true <- true[, c("readn", taxranks)]

formula_str <- sprintf("%s~variable", paste(taxranks, collapse = "+"))

pred <- cast(melt(pred, id = taxranks), formula_str, sum)
pred <- pred[order(pred$readn, decreasing = TRUE), ]

true <- cast(melt(true, id = taxranks), formula_str, sum)
true <- true[order(true$readn, decreasing = TRUE), ]

# ---------------------------------------------------------------------------
# Cumulative share on pred (global, as in bmark_stats)
# ---------------------------------------------------------------------------

pred$csum <- (cumsum(pred$readn) - pred$readn) / sum(pred$readn)

# ---------------------------------------------------------------------------
# Compute stats per taxgroup x rank
# ---------------------------------------------------------------------------

rank_list <- taxranks[seq_len(length(taxranks) - 1)]   # species, genus, family

rows <- list()
for (tg in taxgroups) {
  tmp_tail <- if (tg == "Viruses") tail_vi else if (tg == "Bacteria") tail_ba else 0
  threshold <- (100.0 - tmp_tail) / 100.0

  true_tg <- true[true$division == tg, ]
  pred_tg <- pred[pred$division == tg & pred$csum < threshold, ]

  for (rank in rank_list) {
    true_taxa <- unique(true_tg[[rank]])
    pred_taxa <- unique(pred_tg[[rank]])

    tp <- length(intersect(pred_taxa, true_taxa))
    fp <- length(setdiff(pred_taxa,  true_taxa))
    fn <- length(setdiff(true_taxa,  pred_taxa))
    pr <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    rc <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    f1 <- if ((pr + rc) > 0) 2 * pr * rc / (pr + rc) else 0

    rows[[length(rows) + 1]] <- data.frame(
      taxgroup = tg,
      rank     = rank,
      pred     = length(pred_taxa),
      true     = length(true_taxa),
      TP       = tp,
      FP       = fp,
      FN       = fn,
      Pr       = pr,
      Rc       = rc,
      F1       = f1,
      tail     = tmp_tail,
      stringsAsFactors = FALSE
    )
  }
}

stats <- do.call(rbind, rows)

write.table(stats, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Written:", out_tsv, "\n")
print(stats)
