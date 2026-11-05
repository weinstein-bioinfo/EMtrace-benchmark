#!/usr/bin/env python3
"""Framework 16S — score mock-community runs against their mock-community truth.

Scores every library in one pass. Each manifest row names its library, and each
library brings its own truth table; the abundance basis is always the truth's
`pct_16S_only` column. --truth/--collapse override the registry for a one-off
library that is not in it.

Reads a manifest of runs, scores each one, and writes the two tables every
downstream artifact is built from:

    metrics/framework_16S/metrics.tsv           one row per run x rank x filter cell
    metrics/framework_16S/composition_long.tsv  pipeline, rank, taxon, rel_abundance

The composition table is the single path into the mock-composition plot, so
EMtrace configurations, the Expected truth bar and externally supplied
Emu/Spaghetti tables all arrive in the same shape.

Filters (readn_min, then the tail trim) and the composition bars use the
EM-expected read count `readn_em` where a table has it, like Emu's own 10-read
threshold; other tables fall back to `readn` (`--filter-on readn` forces it).
The count used is recorded per row in `filter_basis`.

Proportions (L1, MeanRelErr, Other %) come from each table's EM estimate where
it has one (EMtrace `prob`), so EMtrace and Emu are both scored on their
EM output; tables without one use `readn`. The basis is recorded in `basis`.

Every run is scored twice on the near-synonym question (the project notes): `collapsed` folds *Shigella* onto *Escherichia* and the
*B. subtilis* group onto *B. subtilis* using the version-controlled map, and
`strict` counts them as distinct taxa. Collapsed numbers are the headline;
strict numbers ship beside them so the artifact penalty stays visible.

Manifest is a TSV of:

    label  [key_settings]  library  panel  format  taxonomy  path  [notes]

`library` selects the truth table (barcode79 -> ZymoBIOMICS, synthgut -> the
synthetic gut mock), so one manifest covers both. `panel` is A, B.1, B.2 or
diag, and is carried into metrics.tsv so a headline table can exclude the
diagnostic rows without parsing free text. Both columns are echoed into
metrics.tsv and composition_long.tsv; a row is identified by label + library,
because the same label legitimately appears once per library.

A `taxonomy` column may carry `silva` for SILVA-reference cells. Omit the
column, or leave it `ncbi`, for everything else.

SILVA has no native species rank: its taxid identifies a lineage down to genus.
Species are scored on SILVA cells only when the table was built with a
species-aware SILVA taxonomy DB (EMtrace `clean_silva_species`), recognisable
by species keys of the form `<taxid>.<k>`. Those species names are DERIVED from
the ENA `organism_name` by a documented procedure — binomial whose genus word
belongs to the SILVA genus, no `sp.`/`uncultured`/strain labels — and are not
SILVA taxonomy entries. Every metrics row records this in `species_source`.
Older SILVA tables, whose species column is an arbitrary organism name, are
still scored at genus and family only.

Usage:
    python3 jobs/framework_16S/framework_16S.score.py --manifest jobs/framework_16S/framework_16S.manifest.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]   # jobs/<framework>/<script>
sys.path.insert(0, str(REPO / "ngs-benchmark"))

from ngs_benchmark.tailtrim import resolve_count_col, tail_trim  # noqa: E402
from ngs_benchmark.detection.abundance import (  # noqa: E402
    abundance_metrics,
    apply_collapse,
    apply_species_normalisation,
    composition_long,
    load_collapse_map,
)
from ngs_benchmark.detection.metrics import set_metrics  # noqa: E402

TRUTH = REPO / "data/truth/zymobiomics_16S_truth.tsv"
COLLAPSE = REPO / "data/truth/zymobiomics_16S_collapse.tsv"
OUTDIR = REPO / "metrics/framework_16S"

# Each library is scored against its own truth. One manifest spans both, so the
# truth cannot be a single global file: scoring synthgut rows against the
# ZymoBIOMICS truth would produce plausible, wrong numbers and no error.
TRUTHS = {
    "barcode79": (TRUTH, COLLAPSE),
    "synthgut": (
        REPO / "data/truth/synthgut_16S_truth.tsv",
        REPO / "data/truth/synthgut_16S_collapse.tsv",
    ),
    # ZymoBIOMICS Gut Microbiome Standard (D6331), R10.4.1 (the project notes).
    # One truth for every replicate library: gutmock_<model>_r<N>.
    "gutmock": (
        REPO / "data/truth/gutmock_16S_truth.tsv",
        REPO / "data/truth/gutmock_16S_collapse.tsv",
    ),
}


def truth_key(library: str) -> str | None:
    """Registry key for a library: exact name, else its prefix before the first
    underscore, so replicate libraries (gutmock_sup_r1) share one truth."""
    if library in TRUTHS:
        return library
    prefix = library.split("_", 1)[0]
    return prefix if prefix in TRUTHS else None

RANKS = ("species", "genus", "family")

# Benchmark axis 4 — these are variables, not settings. readn_min 10 is both the
# legacy default and Emu's own minimum-abundance threshold, so it is the cell
# where our numbers and Emu's are directly comparable.
READN_MIN_SWEEP = (0, 10)
TAIL_SWEEP = (0.0, 0.01)
# Emu's published floor (Curry et al. 2022, Methods): every method's ZymoBIOMICS
# result was "trimmed to include only taxa with abundance >= 0.01%" before the
# metrics, after ZymoBIOMICS' guaranteed contamination level. Unlike the tail
# (a cumulative share), this is a per-taxon relative abundance, applied after
# readn_min; min_rel > 0 is combined with tail 0 only.
MIN_REL_SWEEP = (0.0, 0.0001)
# Benchmark axis 3 — proportions come from the EM estimate (`prob`) where a
# table has one (EMtrace), so both EM tools are scored on their EM output;
# tables without it (Emu's estimated counts, Kraken 2/Bracken, Spaghetti) use
# `readn`. One basis per table; the resolved one is recorded in `basis`.
ABUNDANCE_BASES = ("prob", "readn")   # order of preference

# Count the axis-4 filters (readn_min, then tail) act on, and the weight of the
# composition bars. Default `readn_em`: EMtrace's EM-expected reads
# (prob x reads in the EM input), the same kind of estimate Emu thresholds with
# its 10-read minimum. Tables without it (Emu, Spaghetti, EMtrace runs
# without EM) fall back to `readn`; every metrics row records which count was
# used in `filter_basis`. `--filter-on readn` restores top-hit read counts.
#
# The tail trim is benchmark scoring only, not a pipeline step. It is the
# simplest robust false-positive filter, and it can remove genuinely rare taxa:
# anything among the smallest taxa that together hold `tail` of the reads.
FILTER_ON_CHOICES = ("readn_em", "readn")


# ------------------------------------------------------------------
# Loaders — one per comparator output format
# ------------------------------------------------------------------

def load_emtrace(path: Path) -> pd.DataFrame:
    """EMtrace / LazypipeNP `abund_table.tsv`."""
    # taxid as text: SILVA species keys ("3514.2") must not be parsed as floats.
    df = pd.read_csv(path, sep="\t", dtype={"taxid": str})
    if "readn" not in df.columns:
        raise ValueError(f"{path}: no 'readn' column")
    return df


def load_emu(path: Path) -> pd.DataFrame:
    """Emu `*_rel-abundance.tsv`.

    Emu reports a relative `abundance` and, with --keep-counts, an
    `estimated counts` column. Use counts when present so the readn_min sweep
    means the same thing as it does for EMtrace; otherwise scale the
    relative abundance by the run's read total.

    Emu names the domain column `superkingdom`; it is renamed to `division`,
    which set_metrics filters on — without it every Emu taxon scores as FN.
    """
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"tax_id": "taxid", "superkingdom": "division"})
    # Read-accounting rows (unmapped, mapped_filtered, mapped_unclassified)
    # carry counts but no taxon; they would distort the tail cumsum.
    df = df[pd.to_numeric(df["taxid"], errors="coerce").notna()].copy()
    if "estimated counts" in df.columns:
        df["readn"] = pd.to_numeric(df["estimated counts"], errors="coerce")
    elif "abundance" in df.columns:
        df["readn"] = pd.to_numeric(df["abundance"], errors="coerce")
    else:
        raise ValueError(f"{path}: no 'estimated counts' or 'abundance' column")
    return df


def load_table(path: Path) -> pd.DataFrame:
    """Per-lineage table: one row per taxon, rank columns plus `readn`.

    For comparators whose output we convert ourselves and that carry real read
    counts (Spaghetti, via framework_16S.spaghetti2table.py). Label fixes are
    not the converter's job; they come from the collapse map.
    """
    df = pd.read_csv(path, sep="\t")
    if "readn" not in df.columns:
        raise ValueError(f"{path}: no 'readn' column")
    return df


def load_long(path: Path) -> pd.DataFrame:
    """Pre-scored long table: taxon proportions supplied by a collaborator.

    Expects `rank, taxon, rel_abundance`, optionally `readn` and `division`.
    Used for comparator bars delivered as a per-rank table (Spaghetti, or
    collaborator tables where only the final percentages are available).

    Pivoted so each rank becomes its own column, one row per taxon, which is
    the shape the scorer expects. Only ranks present in the table get a column,
    so e.g. a genus-only table is scored at genus only. Without `readn` the
    proportions stand in for counts, and readn_min filtering is meaningless.
    """
    df = pd.read_csv(path, sep="\t")
    missing = {"rank", "taxon", "rel_abundance"} - set(df.columns)
    if missing:
        raise ValueError(f"{path}: long table missing column(s) {sorted(missing)}")
    if "readn" not in df.columns:
        print(f"[note] {path}: no readn column — readn_min cells use proportions",
              file=sys.stderr)
        df["readn"] = df["rel_abundance"]
    if "division" not in df.columns:
        df["division"] = "Bacteria"
    return pd.concat(
        [
            grp[["taxon", "division", "readn"]].rename(columns={"taxon": rank})
            for rank, grp in df.groupby("rank")
        ],
        ignore_index=True,
    )


LOADERS = {
    "emtrace": load_emtrace,
    "emu": load_emu,
    "table": load_table,
    "long": load_long,
}


# ------------------------------------------------------------------
# Scoring
# ------------------------------------------------------------------

# Provenance of the species labels scored, written to every metrics row.
SPECIES_NATIVE = "native"                         # NCBI taxonomy / comparator's own
SPECIES_SILVA_DERIVED = "silva_organism_name"     # derived, not a SILVA rank
SPECIES_NOT_SCORED = "not_scored"


def has_silva_species_keys(df: pd.DataFrame) -> bool:
    """True if a SILVA abundance table was grouped by derived species keys."""
    if "taxid" not in df.columns:
        return False
    return df["taxid"].astype(str).str.fullmatch(r"\d+\.\d+").any()


def filter_cell(
    df: pd.DataFrame, readn_min: float, tail: float, count_col: str = "readn",
    min_rel: float = 0.0,
) -> pd.DataFrame:
    """Axis 4 — drop taxa below `readn_min`, then below `min_rel` relative
    abundance, then trim the smallest `tail` share.

    All steps act on `count_col`, in that order: the relative abundance and the
    tail share are computed over the taxa that passed `readn_min`.

    The tail itself is `ngs_benchmark.tailtrim`, which is the same code the
    pipeline's `report.tail` runs — this function is where those semantics were
    defined, and it now calls them rather than restating them.
    `preserve_order=False` reproduces what this function returned when the trim
    was inline: sorted by `count_col` descending, with that column coerced to
    numeric.
    """
    counts = pd.to_numeric(df[count_col], errors="coerce").fillna(0)
    keep = counts >= readn_min
    out = df[keep].copy()
    if min_rel > 0:
        kept = counts[keep]
        out = out[(kept / kept.sum()) >= min_rel].copy()
    return tail_trim(out, tail, count_col, preserve_order=False)


def filter_basis(df: pd.DataFrame, filter_on: str) -> str:
    """The count column the filters use for this table: `filter_on`, else readn."""
    return resolve_count_col(df, prefer=filter_on, fallback="readn")


def score_run(
    label: str,
    df: pd.DataFrame,
    truth: pd.DataFrame,
    cmap: pd.DataFrame,
    ranks: tuple[str, ...] = RANKS,
    filter_on: str = "readn_em",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score one run across the synonym x filter grid. Returns (metrics, composition)."""
    truth_taxa = {
        rank: set(truth.loc[truth["pct_16S_only"] > 0, rank]) for rank in ranks
    }
    # Rows at 0 % are expected-ABSENT taxa (Zymo's two yeasts; the synthetic gut
    # mock's B. dentium, which the ONT primer does not amplify). abundance_metrics
    # already drops them with `q = q[q > 0]`, so reporting one lands in other_pct
    # as a false positive, which is the intent. set_metrics knows nothing about
    # abundance, so it must be given the same filtered truth — otherwise a 0 %
    # taxon scores as a true positive when reported and a false negative when
    # not. The Zymo truth hid this because its 0 % rows are Eukaryota and the
    # taxgroups=("Bacteria",) filter removed them anyway.
    truth_present = truth[truth["pct_16S_only"] > 0]
    metric_rows: list[pd.DataFrame] = []
    comp_rows: list[pd.DataFrame] = []
    # Axis 3 — abundance basis: the EM estimate where the table has one, else
    # read counts. Filtering (axis 4) acts on `count_col`.
    basis = next(b for b in ABUNDANCE_BASES if b in df.columns)
    count_col = filter_basis(df, filter_on)

    for synonyms, table in (("strict", df), ("collapsed", apply_collapse(df, cmap))):
        for readn_min in READN_MIN_SWEEP:
            for tail, min_rel in [(t, 0.0) for t in TAIL_SWEEP] + \
                                 [(0.0, r) for r in MIN_REL_SWEEP if r > 0]:
                cell = filter_cell(table, readn_min, tail, count_col, min_rel)
                if cell.empty:
                    continue
                try:
                    d = set_metrics(cell, truth_present, ranks=ranks, taxgroups=("Bacteria",))
                    d = d[["rank", "TP", "FP", "FN", "Pr", "Rc", "F1"]]
                except ValueError:
                    d = None  # comparator table without a division column
                m = abundance_metrics(cell, truth, ranks=ranks, pred_weight=basis)
                if d is not None:
                    m = m.merge(d, on="rank", how="left")
                m.insert(0, "basis", basis)
                m.insert(0, "filter_basis", count_col)
                m.insert(0, "min_rel", min_rel)
                m.insert(0, "tail", tail)
                m.insert(0, "readn_min", readn_min)
                m.insert(0, "synonyms", synonyms)
                m.insert(0, "pipeline", label)
                metric_rows.append(m)

                # One composition bar per run: the headline cell only, weighted
                # by the same count the filters used (EM-adjusted where present).
                if synonyms == "collapsed" and readn_min == 10 and tail == 0.0 and min_rel == 0.0:
                    for rank in ("species", "genus"):
                        if rank not in ranks:
                            continue
                        comp_rows.append(
                            composition_long(
                                cell, rank, label, truth_taxa=truth_taxa[rank],
                                pred_weight=count_col,
                            )
                        )

    metrics = pd.concat(metric_rows, ignore_index=True) if metric_rows else pd.DataFrame()
    comp = pd.concat(comp_rows, ignore_index=True) if comp_rows else pd.DataFrame()
    return metrics, comp


def expected_bar(truth: pd.DataFrame) -> pd.DataFrame:
    """The `Expected` reference bar, on the 16S-only basis."""
    rows = []
    for rank in ("species", "genus"):
        t = truth[truth["pct_16S_only"] > 0]
        g = t.groupby(rank)["pct_16S_only"].sum()
        rows.append(
            pd.DataFrame(
                {
                    "pipeline": "Expected",
                    "rank": rank,
                    "taxon": g.index.astype(str),
                    "rel_abundance": g.to_numpy(dtype=float),
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True,
                    help="TSV of label, library, panel, format, taxonomy, path")
    ap.add_argument("--outdir", default=str(OUTDIR))
    ap.add_argument("--library", default=None,
                    help="score only this library (default: every library in "
                         "the manifest)")
    ap.add_argument("--truth", default=None,
                    help="override the truth table for every scored row. Only "
                         "for a library absent from the registry; with a "
                         "multi-library manifest it would score rows against "
                         "the wrong truth, so it requires --library too.")
    ap.add_argument("--collapse", default=None,
                    help="near-synonym collapse map matching --truth")
    ap.add_argument("--filter-on", choices=FILTER_ON_CHOICES, default="readn_em",
                    help="count the readn_min and tail filters (and the composition "
                         "bars) use; readn_em falls back to readn where a table "
                         "has no EM estimate (default: readn_em)")
    args = ap.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", comment="#")
    if "library" not in manifest.columns:
        # A pre-merge manifest names no library; treat the whole file as one.
        manifest["library"] = args.library or "barcode79"
    if "panel" not in manifest.columns:
        manifest["panel"] = ""
    manifest["library"] = manifest["library"].astype(str)

    if args.truth and not args.library:
        print("--truth requires --library: with a multi-library manifest it "
              "would score every row against one truth.", file=sys.stderr)
        return 2
    if args.library:
        manifest = manifest[manifest["library"] == args.library]
        if manifest.empty:
            print(f"No rows for library '{args.library}'.", file=sys.stderr)
            return 1

    all_metrics, all_comp = [], []
    for library, rows in manifest.groupby("library", sort=True):
        if args.truth:
            truth_p, collapse_p = Path(args.truth), Path(args.collapse or COLLAPSE)
        elif truth_key(library):
            truth_p, collapse_p = TRUTHS[truth_key(library)]
        else:
            print(f"[skip] library '{library}': no truth registered. Add it to "
                  "TRUTHS, or pass --library/--truth.", file=sys.stderr)
            continue
        if not truth_p.exists():
            print(f"[skip] library '{library}': truth not found: {truth_p}",
                  file=sys.stderr)
            continue

        truth = pd.read_csv(truth_p, sep="\t")
        cmap = load_collapse_map(collapse_p)
        print(f"\n== {library}  (truth: {truth_p.name}) ==", file=sys.stderr)

        exp = expected_bar(truth)
        exp.insert(0, "library", library)
        all_comp.append(exp)
        _score_library(library, rows, truth, cmap, args, all_metrics, all_comp)

    return _write(all_metrics, all_comp, Path(args.outdir))


def _score_library(library, manifest, truth, cmap, args, all_metrics, all_comp):
    """Score every manifest row belonging to one library."""
    for row in manifest.itertuples():
        path = REPO / row.path if not Path(row.path).is_absolute() else Path(row.path)
        if not path.exists():
            print(f"[skip] {row.label}: {path} not found", file=sys.stderr)
            continue
        loader = LOADERS.get(row.format)
        if loader is None:
            print(f"[skip] {row.label}: unknown format '{row.format}'", file=sys.stderr)
            continue
        df = loader(path)
        ranks = tuple(r for r in RANKS if r in df.columns)
        species_source = SPECIES_NATIVE if "species" in ranks else SPECIES_NOT_SCORED
        if str(getattr(row, "taxonomy", "")).lower() == "silva":
            # Fold SILVA's "X--other" bins onto their parent at every rank, and
            # reduce species labels to binomials.
            df = apply_species_normalisation(df)
            if "species" not in ranks:
                pass  # genus-only table (Spaghetti): nothing to decide
            elif has_silva_species_keys(df):
                # Species derived from organism_name by the documented rule in
                # EMtrace clean_silva_species — scored, but labelled as such.
                species_source = SPECIES_SILVA_DERIVED
            else:
                # Pre-species DB: each lineage taxid carries one arbitrary name
                # out of thousands, so species is meaningless. See the project notes.
                ranks = tuple(r for r in ranks if r != "species")
                species_source = SPECIES_NOT_SCORED
                print(f"[note] {row.label}: SILVA table without species keys — "
                      "species rank not scored", file=sys.stderr)
        m, c = score_run(row.label, df, truth, cmap, ranks=ranks,
                         filter_on=args.filter_on)
        if not m.empty:
            m.insert(1, "species_source", species_source)
            # label alone is not unique across libraries -- the same cell name
            # is run on both -- so both columns travel with every row.
            m.insert(1, "panel", str(getattr(row, "panel", "") or ""))
            m.insert(1, "library", library)
            # The settings that define the row, from the manifest, for readers of
            # metrics.tsv who should not need the run scripts to decode a label.
            m.insert(1, "key_settings", str(getattr(row, "key_settings", "") or ""))
        if not c.empty:
            c.insert(0, "library", library)
        all_metrics.append(m)
        all_comp.append(c)
        head_rank = "species" if "species" in ranks else "genus"
        head = m[
            (m.synonyms == "collapsed")
            & (m.readn_min == 10)
            & (m["tail"] == 0.0)
            & (m.min_rel == 0.0)
            & (m["rank"] == head_rank)
        ]
        if len(head):
            r = head.iloc[0]
            print(
                f"{row.label:<28} {head_rank:<8} [{species_source}] L1={r.L1:.3f}  "
                f"MeanRelErr={r.MeanRelErr:.3f}  Other={r.other_pct:.1f}%"
            )


def _write(all_metrics, all_comp, outdir: Path) -> int:
    """Write the two combined tables, one row set per library."""
    if not all_metrics:
        print("No runs scored.", file=sys.stderr)
        return 1

    outdir.mkdir(parents=True, exist_ok=True)
    metrics = pd.concat(all_metrics, ignore_index=True)
    comp = pd.concat(all_comp, ignore_index=True)
    metrics.to_csv(outdir / "metrics.tsv", sep="\t", index=False)
    comp.to_csv(outdir / "composition_long.tsv", sep="\t", index=False)

    libs = ", ".join(sorted(metrics["library"].unique()))
    print(f"\nWrote {outdir}/metrics.tsv and {outdir}/composition_long.tsv "
          f"({libs})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
