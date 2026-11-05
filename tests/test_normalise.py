"""Unit tests for §3.1 normalisation."""

from __future__ import annotations

import pandas as pd
import pytest

from ngs_benchmark.normalise import NormaliseOpts, normalise


def _df(**kw) -> pd.DataFrame:
    base = dict(
        readn=[10.0],
        taxid=[1],
        species=["s"],
        genus=["g"],
        family=["f"],
        order=["o"],
        **{"class": ["c"]},
        division=["Viruses"],
    )
    base.update(kw)
    return pd.DataFrame(base)


def test_na_taxon_label_becomes_unknown() -> None:
    df = _df(readn=[10.0, 20.0], taxid=[1, 2], species=[None, ""], genus=["g", "g"],
             family=["f", "f"], order=["o", "o"], **{"class": ["c", "c"]},
             division=["Viruses", "Viruses"])
    out = normalise(df)
    assert (out["species"] == "unknown").all()


def test_phage_drop() -> None:
    df = pd.DataFrame({
        "readn": [10.0, 20.0],
        "taxid": [1, 2],
        "species": ["s1", "s2"],
        "genus": ["g", "g"],
        "family": ["f", "f"],
        "order": ["o", "o"],
        "class": ["c", "c"],
        "division": ["Viruses", "Viruses"],
        "bphage": ["yes", "no"],
    })
    out = normalise(df)
    assert list(out["species"]) == ["s2"]


def test_phage_kept_when_flag_set() -> None:
    df = pd.DataFrame({
        "readn": [10.0, 20.0],
        "taxid": [1, 2],
        "species": ["s1", "s2"],
        "genus": ["g", "g"],
        "family": ["f", "f"],
        "order": ["o", "o"],
        "class": ["c", "c"],
        "division": ["Viruses", "Viruses"],
        "bphage": ["yes", "no"],
    })
    out = normalise(df, NormaliseOpts(keep_phages=True))
    assert set(out["species"]) == {"s1", "s2"}


def test_retro_drop() -> None:
    df = pd.DataFrame({
        "readn": [10.0, 20.0],
        "taxid": [1, 2],
        "species": ["s1", "s2"],
        "genus": ["g", "g"],
        "family": ["Retroviridae", "f"],
        "order": ["o", "o"],
        "class": ["c", "c"],
        "division": ["Viruses", "Viruses"],
    })
    out = normalise(df, NormaliseOpts(exclude_retro=True))
    assert list(out["species"]) == ["s2"]


def test_read_min_filter() -> None:
    df = pd.DataFrame({
        "readn": [5.0, 50.0, 500.0],
        "taxid": [1, 2, 3],
        "species": ["a", "b", "c"],
        "genus": ["g", "g", "g"],
        "family": ["f", "f", "f"],
        "order": ["o", "o", "o"],
        "class": ["c", "c", "c"],
        "division": ["Viruses", "Viruses", "Viruses"],
    })
    out = normalise(df, NormaliseOpts(read_min=10))
    assert set(out["species"]) == {"b", "c"}


def _em_table() -> pd.DataFrame:
    return pd.DataFrame({
        "readn": [12.0, 50.0, 99_938.0],
        "readn_em": [8.0, 50.0, 99_942.0],
        "taxid": [1, 2, 3],
        "species": ["a", "b", "c"],
        "genus": ["g", "g", "g"],
        "family": ["f", "f", "f"],
        "order": ["o", "o", "o"],
        "class": ["c", "c", "c"],
        "division": ["Viruses", "Viruses", "Viruses"],
    })


def test_read_min_on_filter_col() -> None:
    out = normalise(_em_table(), NormaliseOpts(read_min=10, filter_col="readn_em"))
    assert set(out["species"]) == {"b", "c"}


def test_min_rel_floor() -> None:
    # a is 0.012% and b 0.05% of the readn total: both pass a 0.01% floor,
    # neither passes 0.1%.
    out = normalise(_em_table(), NormaliseOpts(min_rel=0.0001))
    assert set(out["species"]) == {"a", "b", "c"}
    out = normalise(_em_table(), NormaliseOpts(min_rel=0.001))
    assert set(out["species"]) == {"c"}


def test_filter_col_falls_back_to_readn() -> None:
    df = _em_table().drop(columns="readn_em")
    out = normalise(df, NormaliseOpts(read_min=10, filter_col="readn_em"))
    assert set(out["species"]) == {"a", "b", "c"}


def test_groupby_sums_duplicate_taxa() -> None:
    df = pd.DataFrame({
        "readn": [10.0, 20.0],
        "taxid": [1, 1],
        "species": ["s1", "s1"],
        "genus": ["g", "g"],
        "family": ["f", "f"],
        "order": ["o", "o"],
        "class": ["c", "c"],
        "division": ["Viruses", "Viruses"],
    })
    out = normalise(df)
    assert len(out) == 1
    assert out["readn"].iloc[0] == 30.0


def test_tail_trim_viruses_drops_smallest_taxa() -> None:
    # 4 viruses: 1000, 100, 10, 1 reads, total 1111. tail_vi is now a
    # FRACTION, so 5 % is 0.05 and the threshold is 0.95.
    #
    # The share is INCLUSIVE (a row's own reads count toward it):
    #   1000/1111 = 0.900, 1100/1111 = 0.990, 1110/1111 = 0.999, 1.000
    # kept: share <= 0.95 -> "a" alone.
    #
    # This expectation CHANGED with the move to the shared tailtrim: the old
    # per-division, percent-unit, preceding-cumsum, strict-< rule kept {a, b}.
    # No completed benchmark result used a non-zero tail, so nothing published
    # moves -- see BMARK/docs.
    df = pd.DataFrame({
        "readn": [1000.0, 100.0, 10.0, 1.0],
        "taxid": [1, 2, 3, 4],
        "species": ["a", "b", "c", "d"],
        "genus": ["g"] * 4,
        "family": ["f"] * 4,
        "order": ["o"] * 4,
        "class": ["c"] * 4,
        "division": ["Viruses"] * 4,
    })
    out = normalise(df, NormaliseOpts(tail_vi=0.05, tail_ba=0.0))
    assert set(out["species"]) == {"a"}


def test_tail_trim_zero_keeps_everything() -> None:
    df = pd.DataFrame({
        "readn": [100.0, 1.0],
        "taxid": [1, 2],
        "species": ["a", "b"],
        "genus": ["g"] * 2,
        "family": ["f"] * 2,
        "order": ["o"] * 2,
        "class": ["c"] * 2,
        "division": ["Viruses"] * 2,
    })
    out = normalise(df, NormaliseOpts(tail_vi=0.0, tail_ba=0.0))
    assert set(out["species"]) == {"a", "b"}


def test_csum_is_computed_and_monotonic() -> None:
    df = pd.DataFrame({
        "readn": [100.0, 50.0, 10.0],
        "taxid": [1, 2, 3],
        "species": ["a", "b", "c"],
        "genus": ["g"] * 3,
        "family": ["f"] * 3,
        "order": ["o"] * 3,
        "class": ["c"] * 3,
        "division": ["Viruses"] * 3,
    })
    out = normalise(df)
    assert "csum" in out.columns
    # Sorted by readn desc; csum starts at 0 (preceding-sum) and grows.
    assert out["csum"].iloc[0] == pytest.approx(0.0)
    assert out["csum"].is_monotonic_increasing


def test_min_rel_exclude_leaves_host_out_of_denominator() -> None:
    df = _em_table()
    df.loc[2, "division"] = "Primates"   # c: 99,938 host reads
    out = normalise(df, NormaliseOpts(min_rel=0.001))
    assert set(out["species"]) == {"c"}
    out = normalise(df, NormaliseOpts(min_rel=0.001, min_rel_exclude=("Primates",)))
    assert set(out["species"]) == {"a", "b", "c"}
