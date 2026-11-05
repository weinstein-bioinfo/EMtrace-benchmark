"""Truth-set loading for the background-filtering benchmark (Module D).

The simmet simulated metagenome ships per-read truth tables keyed by ``readid``.
Rather than load the full 20.7M-read master, we load only the small subset tables
and derive every category by set arithmetic (the simmet contains only human +
bacteria + viruses):

    Homo sapiens               = Total − nohuman
    Bacteria                   = nohuman − vi          (vi ⊆ nohuman)
    Viruses                    = vi
    Viruses excl. Retroviridae = (vi − retro)          = vi_nophix_noretro ∪ phix
    HERV / Retroviridae        = vi − vi_nophix_noretro − phix

``phix`` is identified within the viral table by ``species_id == 10847``
(*Sinsheimervirus phiX174*) and is **retained** in the "Viruses excl.
Retroviridae" category (only Retroviridae/HERV are excluded there).
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

#: Sinsheimervirus phiX174 ``species_id`` in the simmet truth tables.
PHIX_SPECIES_ID = 10847


def _open(path: str | Path):
    """Open a plain or gzipped text file transparently."""
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def _read_id_column(path: str | Path, id_col: str = "readid") -> set[str]:
    """Return the set of values in the ``id_col`` column of a TSV (header required)."""
    ids: set[str] = set()
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = header.index(id_col)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > idx:
                ids.add(parts[idx])
    return ids


def _read_vi_and_phix(
    path: str | Path,
    species_id: int = PHIX_SPECIES_ID,
    id_col: str = "readid",
    sp_col: str = "species_id",
) -> tuple[set[str], set[str]]:
    """Single pass over the viral truth table → ``(all_vi_ids, phix_ids)``."""
    vi: set[str] = set()
    phix: set[str] = set()
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id = header.index(id_col)
        i_sp = header.index(sp_col)
        sp_str = str(species_id)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_id, i_sp):
                continue
            rid = parts[i_id]
            vi.add(rid)
            if parts[i_sp] == sp_str:
                phix.add(rid)
    return vi, phix


def count_data_lines(path: str | Path) -> int:
    """Count non-header data lines in a TSV (used to derive ``total`` from the master)."""
    with _open(path) as fh:
        return max(sum(1 for _ in fh) - 1, 0)


@dataclass
class BgfilterTruth:
    """Truth id-sets plus the read total, with per-category denominators.

    Invariants relied on by the maths: ``vi ⊆ nohuman``,
    ``vi_nophix_noretro ⊆ vi``, ``phix ⊆ vi``, and
    ``vi_nophix_noretro ∩ phix = ∅``.
    """

    total: int
    nohuman: set[str]
    vi: set[str]
    vi_nophix_noretro: set[str]
    phix: set[str]

    @property
    def denominators(self) -> dict[str, int]:
        """Read count per true-origin category (the fraction-removed denominators)."""
        n_nohuman = len(self.nohuman)
        n_vi = len(self.vi)
        n_vinpnr = len(self.vi_nophix_noretro)
        n_phix = len(self.phix)
        return {
            "total": self.total,
            "homo_sapiens": self.total - n_nohuman,
            "bacteria": n_nohuman - n_vi,
            "viruses": n_vi,
            "viruses_excl_retro": n_vinpnr + n_phix,
            "herv": n_vi - n_vinpnr - n_phix,
        }


@dataclass
class DivisionTruth:
    """Per-read division labels for a single-table ground truth.

    Used by the simulated-Nanopore datasets (``combined_<rate>.truth.tsv``), where
    each read carries a ``division`` (Viruses / Bacteria / Eukaryota / …). Unlike
    :class:`BgfilterTruth` (set arithmetic over simmet subset files), the category
    of each read is read directly, so the denominators are just the per-division
    read counts.
    """

    #: readid → division label.
    div_of: dict[str, str]

    @property
    def total(self) -> int:
        return len(self.div_of)

    @property
    def divisions(self) -> list[str]:
        """Distinct division labels, preferred biological order then alphabetical."""
        preferred = ["Viruses", "Bacteria", "Archaea", "Eukaryota"]
        seen = set(self.div_of.values())
        ordered = [d for d in preferred if d in seen]
        ordered += sorted(d for d in seen if d not in preferred)
        return ordered

    @property
    def denominators(self) -> dict[str, int]:
        """``{"total": N, <division>: count, ...}`` — fraction-removed denominators."""
        counts: dict[str, int] = {"total": self.total}
        for div in self.div_of.values():
            counts[div] = counts.get(div, 0) + 1
        return counts


def load_division_truth(
    path: str | Path, id_col: str = "readid", div_col: str = "division"
) -> DivisionTruth:
    """Load a per-read division truth table (``readid``…``division``… TSV)."""
    div_of: dict[str, str] = {}
    with _open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_id = header.index(id_col)
        i_div = header.index(div_col)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(i_id, i_div):
                div_of[f[i_id]] = f[i_div]
    return DivisionTruth(div_of=div_of)


def load_truth(
    *,
    total: int,
    nohuman: str | Path,
    vi: str | Path,
    vi_nophix_noretro: str | Path,
    phix_species_id: int = PHIX_SPECIES_ID,
) -> BgfilterTruth:
    """Load the simmet truth subset tables into a :class:`BgfilterTruth`.

    ``total`` is the full simmet read count (e.g. ``count_data_lines`` of the
    master ``readid_taxpath`` table). Only the three small subset tables are read
    into memory; the human set is never materialised (handled by subtraction).
    """
    vi_ids, phix_ids = _read_vi_and_phix(vi, species_id=phix_species_id)
    return BgfilterTruth(
        total=total,
        nohuman=_read_id_column(nohuman),
        vi=vi_ids,
        vi_nophix_noretro=_read_id_column(vi_nophix_noretro),
        phix=phix_ids,
    )
