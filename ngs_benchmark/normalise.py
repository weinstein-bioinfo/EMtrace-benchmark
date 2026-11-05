"""§3.1 Normalisation pipeline.

Port of `bmark_stats()` (R/NGSlib.R:122). The pipeline is deterministic and
side-effect free — every downstream metric consumes its output.

Order of operations (must not be reordered):

1.  Load via adapter (caller responsibility).
2.  NA/empty taxon labels → literal "unknown".
3.  Drop phages (`bphage == "yes"`) unless --keep-phages.
4.  Drop retroviruses (Retroviridae family) when --exclude-retro.
5.  Drop taxa with readn < --read-min.
6.  Group-by + sum at the top rank (collapses duplicate rows).
7.  csum is computed (descending readn order) for downstream reporting.
8.  Tail trim: drop the smallest taxa whose cumulative share of the library
    exceeds 1 - tail, where tail is `tail_vi` for Viruses and `tail_ba` for
    everything else. Both are FRACTIONS.

The trim itself is `ngs_benchmark.tailtrim`, shared verbatim with the
pipeline's `report.tail`, so "1 % tail trim" means one thing in the pipeline,
in the scorers and in the published numbers. Two consequences worth knowing:

  * The share is computed over EVERY row, not within each division. Grouping
    sets a row's threshold, not its denominator -- so `tail_vi=0` protects
    viral taxa from removal while their reads still count toward everyone
    else's share. On a virus-dominated library that makes the trim harsher on
    the background, not gentler.
  * Step 7's `csum` keeps its own, different convention: PRECEDING cumulative
    share on `readn`. It is a reported column, not the trim basis, and
    changing it would move published tables for no gain.
"""

from __future__ import annotations

from dataclasses import dataclass
import pandas as pd

from ngs_benchmark.schema import TAXRANKS, UNKNOWN
from ngs_benchmark.tailtrim import resolve_count_col, tail_trim


@dataclass(frozen=True)
class NormaliseOpts:
    keep_phages: bool = False
    exclude_retro: bool = False
    # `read_min` is applied to `filter_col`, never to `tail_col` or the score
    # column: a score column may be `prob`, on an entirely different scale, and
    # "--read-min 10" must not silently become "prob >= 10". `filter_col` is a
    # read count (readn, or readn_em falling back to readn).
    read_min: float = 0.0
    # Per-taxon relative-abundance floor, as a FRACTION (0.0001 = 0.01%) of the
    # `filter_col` total that survives `read_min` (Emu's published floor).
    min_rel: float = 0.0
    # Divisions left out of the floor's denominator (e.g. "Primates" for host
    # reads a pipeline reports instead of filtering), so the floor means the
    # same for pipelines that remove host reads and pipelines that keep them.
    min_rel_exclude: tuple[str, ...] = ()
    filter_col: str = "readn"  # count column read_min and min_rel act on
    tail_vi: float = 0.0    # FRACTION of the library to trim from viral taxa
    tail_ba: float = 0.0    # FRACTION to trim from everything else
    tail_col: str = "readn"  # count column the tail trim ranks on


def normalise(df: pd.DataFrame, opts: NormaliseOpts | None = None) -> pd.DataFrame:
    """Apply §3.1 to a canonical AbundanceDF; return a new DataFrame."""
    if opts is None:
        opts = NormaliseOpts()

    out = df.copy()

    # Step 2 — NA / empty → "unknown" at every rank.
    for rank in TAXRANKS:
        if rank in out.columns:
            out[rank] = out[rank].astype("string")
            out[rank] = out[rank].fillna(UNKNOWN)
            out.loc[out[rank].str.fullmatch(r"\s*", na=False), rank] = UNKNOWN

    # Step 3 — drop phages.
    if not opts.keep_phages and "bphage" in out.columns:
        out = out[out["bphage"].astype("string").str.lower() != "yes"]

    # Step 4 — drop retroviruses.
    if opts.exclude_retro and "family" in out.columns:
        out = out[out["family"].astype("string").str.lower() != "retroviridae"]

    # Step 5 — read_min filter, then the relative-abundance floor, both on
    # `filter_col`. Unparseable counts (readn_em "NA" rows) count as 0.
    if opts.read_min > 0 or opts.min_rel > 0:
        fcol = resolve_count_col(out, prefer=opts.filter_col, fallback="readn")
        counts = pd.to_numeric(out[fcol], errors="coerce").fillna(0.0)
        keep = counts >= opts.read_min
        if opts.min_rel > 0:
            in_total = keep
            if opts.min_rel_exclude and "division" in out.columns:
                in_total = keep & ~out["division"].astype("string").isin(opts.min_rel_exclude)
            total = counts[in_total].sum()
            if total > 0:
                keep &= counts / total >= opts.min_rel
        out = out[keep]

    # Step 6 — group by every rank and sum readn. Other numeric columns are summed too;
    # carry first value for string columns.
    if len(out) > 0:
        rank_cols = [r for r in TAXRANKS if r in out.columns]
        agg: dict[str, str] = {"readn": "sum"}
        for col in out.columns:
            if col in rank_cols or col == "readn":
                continue
            if pd.api.types.is_numeric_dtype(out[col]):
                agg[col] = "sum"
            else:
                agg[col] = "first"
        out = out.groupby(rank_cols, dropna=False, as_index=False).agg(agg)

    # Step 7 — csum (preceding-sum cumulative share). Sorted by readn desc.
    out = out.sort_values("readn", ascending=False, kind="mergesort").reset_index(drop=True)
    readn_sum = out["readn"].sum() if len(out) else 0.0
    if readn_sum > 0:
        cs = out["readn"].cumsum() - out["readn"]
        out["csum"] = cs / readn_sum
    else:
        out["csum"] = 0.0

    # Step 8 — tail trim. `tail_ba` is the default threshold and `tail_vi`
    # overrides it for Viruses; the share behind both is global (see the
    # module docstring).
    if opts.tail_vi > 0 or opts.tail_ba > 0:
        out = tail_trim(
            out,
            tail=opts.tail_ba,
            count_col=resolve_count_col(out, prefer=opts.tail_col, fallback="readn"),
            group_col="division",
            group_tails={"Viruses": opts.tail_vi},
            preserve_order=False,
        ).reset_index(drop=True)

    return out
