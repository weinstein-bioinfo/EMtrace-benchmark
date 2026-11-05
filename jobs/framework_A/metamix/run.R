#!/usr/bin/env Rscript
# =============================================================================
# Framework A — metaMix: full 4-step Bayesian mixture model pipeline
# =============================================================================
# Runs metaMix's full pipeline on one sample's DIAMOND output:
#   generative.prob -> reduce.space -> parallel.temper -> bayes.model.aver
#
# Uses metaMix 0.3 (CRAN Archive)
#
# Input is framework_A/metamix/preprocess/diamond/*_diamond.sbatch's DIAMOND blastx
# output in the CUSTOM 10-column format generative.prob(blast.default=FALSE) expects:
#   qseqid qlen sseqid slen stitle bitscore length pident evalue staxids
# staxids comes straight from DIAMOND (its DB was built with --taxonmap), so unlike
# blast.default=TRUE mode, this needs neither a separate read-length file (qlen IS the
# read length, column 2) nor NCBI's prot.accession2taxid mapping at all.
#
# Usage: Rscript run.R <ACC> <RESULTS_DIR> <TAXONOMY_NAMES_DMP>
# <RESULTS_DIR>/<ACC>_diamond.tab must already exist. <TAXONOMY_NAMES_DMP> is NCBI's
# names.dmp, matching the abvh DIAMOND DB's taxonomy source at
# $databases/2026/NCBI/refseq.abvh.1g.234/taxdump/names.dmp.
#
# Writes <RESULTS_DIR>/<ACC>_metamix_species.tsv (taxid, species, Reads, prob,
# log10BF — the columns ngs_benchmark's metamix adapter looks for) plus the full
# step1..step4 .RData bundle for debugging/reanalysis.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript run.R <ACC> <RESULTS_DIR> <TAXONOMY_NAMES_DMP>")
}
acc         <- args[1]
results_dir <- args[2]
taxon_names <- args[3]

diamond_file <- file.path(results_dir, paste0(acc, "_diamond.tab"))
if (!file.exists(diamond_file)) {
  stop(sprintf("Diamond output not found: %s", diamond_file))
}
if (!file.exists(taxon_names)) {
  stop(sprintf("Taxonomy names.dmp not found: %s", taxon_names))
}

library(metaMix)

output_rdata <- file.path(results_dir, paste0(acc, "_metamix_results.RData"))
output_tsv   <- file.path(results_dir, paste0(acc, "_metamix_species.tsv"))

# outDir must be real, not NULL: reduce.space() (via step1$outDir) and step3/step4 all
# auto-save their own *.RData/*.tsv/*.pdf into it.
message(sprintf("[metaMix] %s — step1: generative.prob", acc))
step1 <- generative.prob(
  blast.output.file  = diamond_file,
  blast.default      = FALSE,   # custom 10-col outfmt (qseqid qlen sseqid slen stitle
                                 # bitscore length pident evalue staxids) — avoids both
                                 # of metaMix 0.3's blast.default=TRUE problems.
  contig.weight.file = 1,       # unassembled ONT reads, weight=1 for all reads
  outDir             = results_dir
)

message(sprintf("[metaMix] %s — step2: reduce.space", acc))
step2 <- reduce.space(step1 = step1)

message(sprintf("[metaMix] %s — step3: parallel.temper", acc))
step3 <- parallel.temper(step2 = step2)

message(sprintf("[metaMix] %s — step4: bayes.model.aver", acc))
# poster.prob.thr is left at metaMix 0.3's own published default (0.9) rather than
# overridden. stop()s if no species besides "unknown" clears that threshold — expected for
# a genuinely clean negative control, not a bug; framework_A/metamix/nc_filter.py
# handles the missing output.
step4 <- bayes.model.aver(step2 = step2, step3 = step3, taxon.name.map = taxon_names)

save(step1, step2, step3, step4, file = output_rdata)

species <- step4$presentSpecies.allInfo
out <- data.frame(
  taxid   = species$taxonID,
  species = species$scientName,
  Reads   = species$finalAssignments,
  prob    = species$poster.prob,
  log10BF = species$log10BF,
  stringsAsFactors = FALSE
)
write.table(out, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

message(sprintf("[metaMix] %s — done: %s (%d species incl. unknown)", acc, output_tsv, nrow(out)))
