"""Tests for the comparator adapters (CZ ID, EPI2ME, Kraken/Bracken, MEGAN-LR, metaMix).

Where possible the tests are hermetic (synthetic fixtures, no NCBI DB). Tests
that need the NCBI taxonomy DB are wrapped in `_needs_ncbi_db` and skipped
when the DB cannot be loaded.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from ngs_benchmark.io.adapters import czid, epi2me, kraken_bracken, megan_lr, metamix


def _ncbi_db_available() -> bool:
    try:
        from ngs_benchmark.io.taxonomy import get_taxonomy
        get_taxonomy()
        return True
    except (ImportError, FileNotFoundError):
        return False


_needs_ncbi_db = pytest.mark.skipif(
    not _ncbi_db_available(),
    reason="NCBI taxonomy DB not available; set NGS_BMARK_NCBI_TAXONOMY_DB to enable.",
)


# ---------------------------------------------------------------------------
# CZ ID — no NCBI lookup required; file already carries the full lineage.
# ---------------------------------------------------------------------------


def test_czid_uses_nt_count_when_present(tmp_path: Path) -> None:
    p = tmp_path / "czid.tsv"
    pd.DataFrame(
        {
            "taxid": [562],
            "name": ["E. coli"],
            "is_phage": ["FALSE"],
            "nt_count": [100],
            "nr_count": [50],
            "species": ["Escherichia coli"],
            "genus": ["Escherichia"],
            "family": ["Enterobacteriaceae"],
            "order": ["Enterobacterales"],
            "class": ["Gammaproteobacteria"],
            "division": ["Bacteria"],
        }
    ).to_csv(p, sep="\t", index=False)
    df = czid.load(p)
    assert df["readn"].iloc[0] == 100.0


def test_czid_falls_back_to_nr_count_when_nt_missing(tmp_path: Path) -> None:
    p = tmp_path / "czid.tsv"
    pd.DataFrame(
        {
            "taxid": [562],
            "name": ["E. coli"],
            "is_phage": ["FALSE"],
            "nt_count": [None],   # NA
            "nr_count": [50],
            "species": ["Escherichia coli"],
            "genus": ["Escherichia"],
            "family": ["Enterobacteriaceae"],
            "order": ["Enterobacterales"],
            "class": ["Gammaproteobacteria"],
            "division": ["Bacteria"],
        }
    ).to_csv(p, sep="\t", index=False)
    df = czid.load(p)
    assert df["readn"].iloc[0] == 50.0


def test_czid_takes_larger_count_when_nr_exceeds_nt(tmp_path: Path) -> None:
    """Buddle DNA CMV case: nt 1 / nr 75 must read as 75, not 1."""
    p = tmp_path / "czid.tsv"
    pd.DataFrame(
        {
            "taxid": [10359],
            "name": ["Human betaherpesvirus 5"],
            "is_phage": ["FALSE"],
            "nt_count": [1],
            "nr_count": [75],
            "species": ["Cytomegalovirus humanbeta5"],
            "genus": ["Cytomegalovirus"],
            "family": ["Orthoherpesviridae"],
            "order": ["Herpesvirales"],
            "class": ["Herviviricetes"],
            "division": ["Viruses"],
        }
    ).to_csv(p, sep="\t", index=False)
    df = czid.load(p)
    assert df["readn"].iloc[0] == 75.0


def test_czid_drops_is_phage_rows(tmp_path: Path) -> None:
    p = tmp_path / "czid.tsv"
    pd.DataFrame(
        {
            "taxid": [1, 2],
            "name": ["v1", "phage1"],
            "is_phage": ["FALSE", "TRUE"],
            "nt_count": [100, 100],
            "nr_count": [0, 0],
            "species": ["v1", "phage1"],
            "genus": ["g1", "p1"],
            "family": ["f1", "pf1"],
            "order": ["o1", "po1"],
            "class": ["c1", "pc1"],
            "division": ["Viruses", "Viruses"],
        }
    ).to_csv(p, sep="\t", index=False)
    df = czid.load(p)
    assert list(df["species"]) == ["v1"]


# ---------------------------------------------------------------------------
# EPI2ME / Bracken
# ---------------------------------------------------------------------------


@_needs_ncbi_db
def test_epi2me_uses_new_est_reads(tmp_path: Path) -> None:
    p = tmp_path / "sample.bracken.tsv"
    pd.DataFrame(
        {
            "name": ["Escherichia coli"],
            "taxonomy_id": [562],
            "taxonomy_lvl": ["S"],
            "kraken_assigned_reads": [80],
            "added_reads": [20],
            "new_est_reads": [100],
            "fraction_total_reads": [0.5],
        }
    ).to_csv(p, sep="\t", index=False)
    df = epi2me.load(p)
    assert df["readn"].iloc[0] == 100.0
    assert df["taxid"].iloc[0] == 562
    assert df["genus"].iloc[0] == "Escherichia"
    assert df["division"].iloc[0] == "Bacteria"


@_needs_ncbi_db
def test_epi2me_filters_to_species_rank(tmp_path: Path) -> None:
    p = tmp_path / "sample.bracken.tsv"
    pd.DataFrame(
        {
            "name": ["Escherichia coli", "Escherichia"],
            "taxonomy_id": [562, 561],
            "taxonomy_lvl": ["S", "G"],     # only S rows should survive
            "new_est_reads": [100, 200],
        }
    ).to_csv(p, sep="\t", index=False)
    df = epi2me.load(p)
    assert len(df) == 1
    assert df["taxid"].iloc[0] == 562


# ---------------------------------------------------------------------------
# Kraken2 / Bracken
# ---------------------------------------------------------------------------


@_needs_ncbi_db
def test_kraken_kreport_species_only(tmp_path: Path) -> None:
    p = tmp_path / "kreport.txt"
    p.write_text(
        # %      reads_clade  reads_taxon  rank  taxid  name
        "100.0\t1000\t0\tR\t1\troot\n"
        "50.0\t500\t10\tD\t2\t  Bacteria\n"
        "30.0\t300\t5\tG\t561\t    Escherichia\n"
        "20.0\t200\t200\tS\t562\t      Escherichia coli\n"
        "10.0\t100\t100\tS\t1351\t    Enterococcus faecalis\n"
    )
    df = kraken_bracken.load(p)
    assert set(df["taxid"]) == {562, 1351}
    assert int(df.loc[df["taxid"] == 562, "readn"].iloc[0]) == 200
    assert int(df.loc[df["taxid"] == 1351, "readn"].iloc[0]) == 100
    # NCBI fill-in works.
    assert df.loc[df["taxid"] == 562, "species"].iloc[0] == "Escherichia coli"


@_needs_ncbi_db
def test_kraken_bracken_delegates_when_header_matches(tmp_path: Path) -> None:
    p = tmp_path / "merged.tsv"
    pd.DataFrame(
        {
            "name": ["Escherichia coli"],
            "taxonomy_id": [562],
            "taxonomy_lvl": ["S"],
            "kraken_assigned_reads": [80],
            "added_reads": [20],
            "new_est_reads": [100],
            "fraction_total_reads": [0.5],
        }
    ).to_csv(p, sep="\t", index=False)
    df = kraken_bracken.load(p)
    assert df["readn"].iloc[0] == 100.0
    assert df["taxid"].iloc[0] == 562


# ---------------------------------------------------------------------------
# MEGAN-LR
# ---------------------------------------------------------------------------


@_needs_ncbi_db
def test_megan_lr_taxid_first(tmp_path: Path) -> None:
    p = tmp_path / "megan.tsv"
    pd.DataFrame({"TaxId": [562, 1351], "Count": [100, 50]}).to_csv(p, sep="\t", index=False)
    df = megan_lr.load(p)
    assert set(df["taxid"]) == {562, 1351}
    assert df.loc[df["taxid"] == 562, "species"].iloc[0] == "Escherichia coli"


@_needs_ncbi_db
def test_megan_lr_name_first(tmp_path: Path) -> None:
    p = tmp_path / "megan.tsv"
    pd.DataFrame(
        {"Taxonomy": ["Escherichia coli"], "Count": [100]}
    ).to_csv(p, sep="\t", index=False)
    df = megan_lr.load(p)
    assert df["readn"].iloc[0] == 100.0
    assert df["species"].iloc[0] == "Escherichia coli"


# ---------------------------------------------------------------------------
# metaMix
# ---------------------------------------------------------------------------


@_needs_ncbi_db
def test_metamix_with_taxid_column(tmp_path: Path) -> None:
    p = tmp_path / "metamix.tsv"
    pd.DataFrame(
        {
            "Species": ["Escherichia coli"],
            "taxid": [562],
            "Reads": [123],
            "MarginalLogLike": [-45.6],
        }
    ).to_csv(p, sep="\t", index=False)
    df = metamix.load(p)
    assert df["readn"].iloc[0] == 123.0
    assert df["taxid"].iloc[0] == 562
    assert df["division"].iloc[0] == "Bacteria"
    assert df["prob"].iloc[0] == pytest.approx(-45.6)


@_needs_ncbi_db
def test_metamix_name_only(tmp_path: Path) -> None:
    p = tmp_path / "metamix.tsv"
    pd.DataFrame(
        {"Species": ["Escherichia coli"], "Reads": [50]}
    ).to_csv(p, sep="\t", index=False)
    df = metamix.load(p)
    assert df["readn"].iloc[0] == 50.0
    assert df["species"].iloc[0] == "Escherichia coli"


def test_metamix_missing_required_column_errors(tmp_path: Path) -> None:
    p = tmp_path / "metamix.tsv"
    pd.DataFrame({"Species": ["x"]}).to_csv(p, sep="\t", index=False)
    with pytest.raises(ValueError, match="read-count column"):
        metamix.load(p)
