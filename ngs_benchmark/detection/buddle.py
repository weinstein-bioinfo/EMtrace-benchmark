"""§3.5 Buddle 2024 PMR/RPMR threshold mode.

Filters a prediction DataFrame using the rule from Buddle et al. 2024
(Genome Medicine 16:111):

    For viruses:        PMR > 0.01% AND (RPMR > 5 OR readn_neg_ctrl == 0)
    For other organisms: PMR > 1%    AND  RPMR > 10

Where:
    PMR  = readn_taxon / readn_total_microbial * 100
    RPMR = RPM_sample / RPM_neg_ctrl
    RPM  = readn / total_reads * 1e6

The total-reads counts come from the sample's library QC; in the Buddle
truth-file flavour they are columns of the truth file or are supplied via
explicit kwargs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
import pandas as pd

from ngs_benchmark.schema import MICROBIAL_DIVISIONS


@dataclass(frozen=True)
class BuddleOpts:
    """Knobs for the §3.5 rule. Defaults match Buddle et al. 2024."""

    pmr_virus: float = 0.01     # %
    pmr_other: float = 1.0      # %
    rpmr_virus: float = 5.0
    rpmr_other: float = 10.0
    microbial_divisions: tuple[str, ...] = field(default_factory=lambda: tuple(MICROBIAL_DIVISIONS))


def _rpm(readn: pd.Series, total_reads: float) -> pd.Series:
    if total_reads <= 0:
        return pd.Series(np.zeros(len(readn)), index=readn.index)
    return readn / total_reads * 1e6


def annotate(
    sample: pd.DataFrame,
    neg_ctrl: pd.DataFrame,
    total_reads_sample: float,
    total_reads_nc: float,
    opts: BuddleOpts | None = None,
    rank: str = "species",
) -> pd.DataFrame:
    """Return `sample` with `PMR`, `RPMR`, `buddle_pass` columns added.

    `neg_ctrl` is the parallel prediction DataFrame for the matched negative
    control. Both should be §3.1-normalised. Taxa not present in `neg_ctrl`
    are assigned readn_neg_ctrl = 0 ⇒ RPMR = +inf (rule passes).
    """
    if opts is None:
        opts = BuddleOpts()

    out = sample.copy()

    # PMR — denominator is the sum of microbial reads in the sample.
    if "division" in out.columns:
        microbial_mask = out["division"].isin(opts.microbial_divisions)
        total_microbial = float(out.loc[microbial_mask, "readn"].sum())
    else:
        total_microbial = float(out["readn"].sum())
    if total_microbial > 0:
        out["PMR"] = out["readn"] / total_microbial * 100.0
    else:
        out["PMR"] = 0.0

    # Per-taxon RPMR — match on `rank` (the manuscript-relevant taxonomic level).
    rpm_sample = _rpm(out["readn"], total_reads_sample)
    nc_lookup = (
        neg_ctrl.groupby(rank, dropna=False)["readn"].sum() if len(neg_ctrl) else pd.Series(dtype=float)
    )
    readn_nc = out[rank].map(nc_lookup).fillna(0.0).astype(float)
    rpm_nc = _rpm(readn_nc, total_reads_nc)
    with np.errstate(divide="ignore", invalid="ignore"):
        rpmr = np.where(rpm_nc > 0, rpm_sample / rpm_nc, np.inf)
    out["readn_neg_ctrl"] = readn_nc.to_numpy()
    out["RPMR"] = rpmr

    # Apply rule.
    is_virus = (out.get("division") == "Viruses") if "division" in out.columns else False
    if isinstance(is_virus, bool):
        # No division column — treat every row as "other".
        out["buddle_pass"] = (out["PMR"] > opts.pmr_other) & (out["RPMR"] > opts.rpmr_other)
    else:
        virus_pass = (out["PMR"] > opts.pmr_virus) & (
            (out["RPMR"] > opts.rpmr_virus) | (out["readn_neg_ctrl"] == 0)
        )
        other_pass = (out["PMR"] > opts.pmr_other) & (out["RPMR"] > opts.rpmr_other)
        out["buddle_pass"] = np.where(is_virus, virus_pass, other_pass)

    return out


def filter_by_rule(annotated: pd.DataFrame) -> pd.DataFrame:
    """Drop rows where buddle_pass is False; return the survivor DataFrame."""
    if "buddle_pass" not in annotated.columns:
        raise ValueError("Expected `buddle_pass` column; call `annotate()` first.")
    return annotated[annotated["buddle_pass"]].drop(columns=["buddle_pass"]).reset_index(drop=True)
