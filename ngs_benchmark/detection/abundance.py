"""Abundance-accuracy metrics for composition benchmarks.

Detection metrics (`metrics.set_metrics`) score *which* taxa were reported.
These score *how much* of each — the half that matters for a mock community
where every method finds all 8 organisms and the question is proportion.

Two families of metric are produced, because the two literatures differ:

``L1``
    Sum of absolute differences between predicted and true proportions over
    the **union** of predicted and true taxa. This is the metric Emu is scored
    on (Curry et al. 2022, Nat Methods 19:845-853), so it is the one that makes
    our numbers comparable to a published bar. Range 0-2; lower is better.

``MeanRelErr``
    Sum of absolute differences over the **truth taxa only** — the legacy
    LazypipeNP metric (``Bmark.16S.mockcomm.R``). Despite the name it is not a
    mean and not relative; it is L1 restricted to expected taxa. Kept for
    continuity with the legacy tables, and because it is blind to false
    positives it should never be reported without ``other_pct`` beside it.

``other_pct``
    Percentage of predicted abundance mass sitting on taxa absent from the
    truth table. This is the "Other" bar of the mock-composition plot, and it
    is the single number the legacy comparison turned on (Emu 0.7 %, legacy
    LazypipeNP 8.8-47.8 %).

Taxonomic near-synonyms (*Shigella* for *Escherichia*, the *B. subtilis*
group) are database artifacts rather than genuine errors. `apply_collapse`
folds them onto their parent using an explicit, version-controlled map so the
choice is visible rather than buried in the scoring code. Out-of-date names in
a comparator's reference (*Lactobacillus fermentum*) are fixed the same way,
not in the comparator's converter.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from ngs_benchmark.schema import TAXRANKS


# ------------------------------------------------------------------
# Collapse map
# ------------------------------------------------------------------

COLLAPSE_COLS = ("rank", "from", "to")


def load_collapse_map(path: str | Path) -> pd.DataFrame:
    """Load a TSV of ``rank, from, to [, if_species] [, reason]`` near-synonym rules.

    `if_species`, when set, scopes a rule to rows whose *original* species
    label equals it. That is how a rename is applied to one species' genus
    (``genus Lactobacillus -> Limosilactobacillus if_species Lactobacillus
    fermentum``) without folding every other member of the old genus.
    """
    cmap = pd.read_csv(path, sep="\t", header=0, dtype="string")
    missing = [c for c in COLLAPSE_COLS if c not in cmap.columns]
    if missing:
        raise ValueError(
            f"{path}: collapse map missing column(s) {missing}. "
            f"Found: {list(cmap.columns)}"
        )
    bad = sorted(set(cmap["rank"].dropna()) - set(TAXRANKS))
    if bad:
        raise ValueError(f"{path}: unknown rank(s) in collapse map: {bad}")
    if "if_species" not in cmap.columns:
        cmap["if_species"] = pd.NA
    cmap["if_species"] = cmap["if_species"].replace("", pd.NA)
    return cmap


def apply_collapse(df: pd.DataFrame, cmap: pd.DataFrame) -> pd.DataFrame:
    """Rewrite near-synonym labels onto their parent, per rank.

    Returns a new DataFrame; rows are relabelled, never dropped or merged, so
    downstream grouping decides how the folded rows combine.

    Scoped rules (`if_species` set) are matched against the species labels as
    they were *before* any species rule ran, so the order of rules in the map
    does not matter. A table without a `species` column cannot satisfy a scope,
    and scoped rules are skipped for it.
    """
    out = df.copy()
    if "if_species" not in cmap.columns:
        cmap = cmap.assign(if_species=pd.NA)
    orig_species = df["species"] if "species" in df.columns else None
    for rank, grp in cmap.groupby("rank", dropna=True):
        if rank not in out.columns:
            continue
        before = out[rank].copy()
        plain = grp[grp["if_species"].isna()]
        lookup = dict(zip(plain["from"], plain["to"]))
        out[rank] = before.map(lambda v: lookup.get(v, v))
        if orig_species is None:
            continue
        scoped = grp[grp["if_species"].notna()]
        for src, dst, sp in zip(scoped["from"], scoped["to"], scoped["if_species"]):
            out.loc[(before == src) & (orig_species == sp), rank] = dst
    return out


# ------------------------------------------------------------------
# SILVA species normalisation
# ------------------------------------------------------------------

# A binomial, optionally with a strain designation after it. SILVA has no
# species rank, so the EMtrace SILVA-taxonomy path fills `species` from the
# free-text `organism_name` column. Only ~9.5 % of those are plain binomials;
# the rest carry strain suffixes ("Bacillus subtilis BEST7003"), which score as
# distinct taxa and inflate `Other` for reasons that have nothing to do with
# the pipeline.
_BINOMIAL = re.compile(r"^([A-Z][a-z]+) ([a-z][a-z-]+)(?:\s|$)")

# Labels that carry no species information at all. Left as-is rather than
# normalised, so they stay visible as unassigned mass.
_UNINFORMATIVE = re.compile(
    r"^(uncultured|unidentified|unclassified|metagenome|bacterium|archaeon)\b",
    re.IGNORECASE,
)


def normalise_species(name: str | None) -> str | None:
    """Reduce a free-text organism label to its binomial, where it has one.

    ``"Bacillus subtilis BEST7003"`` -> ``"Bacillus subtilis"``. Labels with no
    binomial (``"uncultured bacterium"``, ``"Escherichia sp. K-12"``) are
    returned unchanged so they remain countable as unassigned mass rather than
    being silently folded into a real species.

    Apply this only to tables produced under the SILVA taxonomy backbone. NCBI
    output already carries a curated species rank and needs no normalisation.
    """
    if not isinstance(name, str) or not name:
        return name
    if _UNINFORMATIVE.match(name):
        return name
    m = _BINOMIAL.match(name)
    if not m:
        return name
    genus, epithet = m.groups()
    # "Escherichia sp." and friends name a genus, not a species.
    if epithet in ("sp", "spp", "cf", "aff"):
        return name
    return f"{genus} {epithet}"


# A species label counts as NAMED only if it is a real binomial. NCBI's
# bracketed provisional genera ("[Clostridium] innocuum", "[Ruminococcus]
# gnavus") ARE named — the brackets flag a misclassified genus, not a missing
# name — so the genus word may be wrapped in them.
_NAMED_SPECIES = re.compile(r"^\[?([A-Z][a-z]+)\]? ([a-z][a-z-]+)(?:\s|$)")


def is_named_species(name: object) -> bool:
    """True when `name` resolves to a genus + epithet binomial.

    False for ``"Escherichia sp."``, ``"uncultured bacterium"``, empty labels
    and anything non-string. Mass on those entries is reported separately by
    ``abundance_metrics`` as ``nospname_pct``: moving reads onto an unnamed bin
    lowers species L1 without classifying anything better, so the two numbers
    have to be read together (the project notes).
    """
    if not isinstance(name, str) or not name.strip():
        return False
    if _UNINFORMATIVE.match(name):
        return False
    m = _NAMED_SPECIES.match(name)
    if not m:
        return False
    return m.group(2) not in ("sp", "spp", "cf", "aff")


# SILVA marks members of a taxon that are not resolved to a named child with a
# "--other" suffix, at every rank: "Salmonella--other" sits beside "Salmonella"
# (22,868 rows against 90,790 in SILVA 144). 601 of 7,086 genera carry it, plus
# 106 families and 38 orders. Left alone it splits a single organism across two
# bins — the Goal 3 SILVA smoke run reported Salmonella at 4.1 % and
# Salmonella--other at 3.2 % instead of one 7.3 % bin.
_SILVA_OTHER = re.compile(r"--other$")


def strip_silva_other(name: str | None) -> str | None:
    """Fold SILVA's ``X--other`` bin onto its parent taxon ``X``."""
    if not isinstance(name, str):
        return name
    return _SILVA_OTHER.sub("", name)


def apply_species_normalisation(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise SILVA labels. Apply to SILVA-taxonomy tables only.

    Two mechanical fixes, both labelling artifacts rather than pipeline
    behaviour: `species` is reduced to a binomial, and ``--other`` bins are
    folded onto their parent at every rank.

    NCBI output carries a curated species rank and no ``--other`` convention,
    so it must not be passed through here.
    """
    out = df.copy()
    for rank in TAXRANKS:
        if rank in out.columns:
            out[rank] = out[rank].map(strip_silva_other)
    if "species" in out.columns:
        out["species"] = out["species"].map(normalise_species)
    return out


# ------------------------------------------------------------------
# Proportions
# ------------------------------------------------------------------

def to_proportions(
    df: pd.DataFrame,
    rank: str,
    weight_col: str = "readn",
) -> pd.Series:
    """Collapse a table to ``taxon -> proportion`` at one rank.

    `weight_col` selects the abundance basis — benchmark axis 3. Use ``readn``
    for raw read counts and ``prob`` for the EM posterior. Rows whose label at
    this rank is missing are dropped, so proportions are over labelled mass
    only.
    """
    if rank not in df.columns:
        raise ValueError(f"Rank column '{rank}' missing. Found: {list(df.columns)}")
    if weight_col not in df.columns:
        raise ValueError(
            f"Weight column '{weight_col}' missing. Found: {list(df.columns)}"
        )
    sub = df[[rank, weight_col]].copy()
    sub[rank] = sub[rank].astype("string")
    sub = sub[sub[rank].notna() & (sub[rank].str.len() > 0)]
    sub[weight_col] = pd.to_numeric(sub[weight_col], errors="coerce").fillna(0.0)
    grouped = sub.groupby(rank, dropna=True)[weight_col].sum()
    total = grouped.sum()
    if total <= 0:
        return grouped.astype(float)
    return (grouped / total).astype(float)


# ------------------------------------------------------------------
# Metrics
# ------------------------------------------------------------------

def abundance_metrics(
    pred: pd.DataFrame,
    truth: pd.DataFrame,
    ranks: tuple[str, ...] = ("species", "genus", "family"),
    pred_weight: str = "readn",
    truth_weight: str = "pct_16S_only",
) -> pd.DataFrame:
    """Compute L1, MeanRelErr and other_pct per rank.

    `truth_weight` names the abundance basis of the truth table. For the
    ZymoBIOMICS mock use ``pct_16S_only``: a 16S read count measures rRNA
    operon copies, not cells or genomes, so the 16S-only column is the only
    defensible basis. The two yeasts are 0 % on that basis and therefore drop
    out of the truth proportions automatically, while any yeast call in `pred`
    still lands in `other_pct` as a specificity failure.
    """
    rows: list[dict] = []
    for rank in ranks:
        p = to_proportions(pred, rank, pred_weight)
        q = to_proportions(truth, rank, truth_weight)
        q = q[q > 0]

        union = p.index.union(q.index)
        p_u = p.reindex(union, fill_value=0.0)
        q_u = q.reindex(union, fill_value=0.0)

        l1 = float((p_u - q_u).abs().sum())
        mre = float((p.reindex(q.index, fill_value=0.0) - q).abs().sum())
        other = float(p[~p.index.isin(q.index)].sum())

        missed = sorted(q.index[~q.index.isin(p.index)])
        # Predicted mass sitting on entries with no species name. Only defined
        # at species rank; genus/family labels are not binomials.
        if rank == "species" and len(p):
            unnamed = float(p[[not is_named_species(v) for v in p.index]].sum())
            nospname = 100.0 * unnamed
        else:
            nospname = float("nan")
        rows.append(
            {
                "rank": rank,
                "n_pred": int(len(p)),
                "n_true": int(len(q)),
                "L1": l1,
                "MeanRelErr": mre,
                "other_pct": 100.0 * other,
                "nospname_pct": nospname,
                "n_missing": len(missed),
                "missing": ";".join(missed),
            }
        )
    return pd.DataFrame(rows)


def composition_long(
    pred: pd.DataFrame,
    rank: str,
    pipeline: str,
    truth_taxa: set[str] | None = None,
    pred_weight: str = "readn",
    other_label: str = "Other",
) -> pd.DataFrame:
    """Long-format rows for the stacked mock-composition plot.

    Every taxon outside `truth_taxa` is pooled into a single `Other` row, which
    is what makes the bars readable and what the legacy figure compared on. Pass
    `truth_taxa=None` to keep every taxon separate.

    Returns columns ``pipeline, rank, taxon, rel_abundance`` (percent), the
    single shape that EMtrace configurations, the Expected truth bar and
    externally supplied Emu/Spaghetti tables all enter the plot through.
    """
    p = to_proportions(pred, rank, pred_weight) * 100.0
    if truth_taxa is not None:
        keep = p[p.index.isin(truth_taxa)]
        other = float(p[~p.index.isin(truth_taxa)].sum())
        p = keep
        if other > 0:
            p = pd.concat([p, pd.Series({other_label: other})])
    return pd.DataFrame(
        {
            "pipeline": pipeline,
            "rank": rank,
            "taxon": p.index.astype(str),
            "rel_abundance": p.to_numpy(dtype=float),
        }
    )
