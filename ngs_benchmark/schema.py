"""Canonical AbundanceDF schema constants.

Every adapter normalises its input to this schema so the downstream metric code
can be schema-agnostic.
"""

from __future__ import annotations

# Ordered low → high; matches Lazypipe/EMtrace `abund_table.tsv`.
TAXRANKS = ("species", "genus", "family", "order", "class", "division")

# Required columns that every loader must produce.
REQUIRED_COLS = ("readn", "taxid", *TAXRANKS)

# Optional columns — passed through when present, absent otherwise.
OPTIONAL_COLS = (
    "contign",   # legacy
    "csum",      # legacy cumulative-share; superseded by `prob` in newer output
    "prob",      # EMtrace EM posterior
    "pvalue",    # EMtrace EM likelihood-ratio P-value
    "bphage",    # phage flag ("yes"/"no"); used by §3.1 step 3
    "sample_id", # batch / dilution-series mode
)

# Top-level taxonomic groups recognised at the `division` rank.
# Match the values that appear in Lazypipe/EMtrace output.
TAXGROUPS_DEFAULT = ("Viruses", "Bacteria")
# Set of values treated as "microbial" for the Buddle PMR denominator (§3.5).
MICROBIAL_DIVISIONS = ("Viruses", "Bacteria", "Archaea")

# Short labels used in output filenames / tab names (§2.3).
TAXGROUP_SHORT = {
    "Viruses": "vi",
    "Bacteria": "ba",
    "Eukaryota": "euk",
    "Archaea": "ar",
}

UNKNOWN = "unknown"
