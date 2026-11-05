"""`--tail` / `--tail-col`, and the deprecated percent spellings.

`--tail-vi` / `--tail-ba` predate the unification and mean PERCENT. They keep
that meaning so job scripts written against them do not silently ask for a
100x larger trim; `--tail` is the fraction everything else uses.
"""

from __future__ import annotations

import warnings

import pandas as pd
import pytest

from ngs_benchmark.cli import _resolve_tails
from ngs_benchmark.detection.runner import DetectionConfig, _make_norm_opts
from ngs_benchmark.normalise import NormaliseOpts, normalise
from ngs_benchmark.realtime.runner import RealtimeConfig, score_snapshot


# ------------------------------------------------------------------ #
# _resolve_tails
# ------------------------------------------------------------------ #

def test_global_tail_applies_to_both_divisions():
    with warnings.catch_warnings():
        warnings.simplefilter("error")          # must not warn
        assert _resolve_tails(0.02, None, None) == (0.02, 0.02)


def test_deprecated_flags_convert_from_percent_and_warn():
    with pytest.warns(DeprecationWarning, match="PERCENT"):
        vi, ba = _resolve_tails(0.0, 5.0, None)
    assert vi == 0.05
    assert ba == 0.0


def test_deprecated_flags_override_the_global_tail():
    with pytest.warns(DeprecationWarning):
        assert _resolve_tails(0.02, 5.0, None) == (0.05, 0.02)


def test_the_spelling_every_job_script_uses_still_disables():
    # `--tail-vi 0 --tail-ba 0` is what all 84 completed detection runs passed.
    with pytest.warns(DeprecationWarning):
        assert _resolve_tails(0.0, 0.0, 0.0) == (0.0, 0.0)


# ------------------------------------------------------------------ #
# Defaults must not trim
# ------------------------------------------------------------------ #

def _mixed():
    return pd.DataFrame({
        "readn": [1000.0, 100.0, 10.0, 1.0],
        "taxid": [1, 2, 3, 4],
        "species": ["a", "b", "c", "d"],
        "genus": ["g"] * 4, "family": ["f"] * 4,
        "order": ["o"] * 4, "class": ["c"] * 4,
        "division": ["Viruses", "Viruses", "Bacteria", "Bacteria"],
    })


def test_default_opts_do_not_trim():
    # tail_ba used to default to 1.0 (i.e. 1 %); it is now 0.0, off.
    assert NormaliseOpts().tail_ba == 0.0
    assert NormaliseOpts().tail_vi == 0.0
    assert len(normalise(_mixed())) == 4


def test_default_detection_config_does_not_trim():
    cfg = DetectionConfig(pred="p", truth="t", outdir="o")
    assert (cfg.tail_vi, cfg.tail_ba, cfg.tail_col) == (0.0, 0.0, "readn")
    assert _make_norm_opts(cfg).tail_ba == 0.0


# ------------------------------------------------------------------ #
# Grouping and tail_col in normalise
# ------------------------------------------------------------------ #

def test_viruses_can_be_protected_while_the_rest_is_trimmed():
    out = normalise(_mixed(), NormaliseOpts(tail_vi=0.0, tail_ba=0.01))
    assert set(out[out["division"] == "Viruses"]["species"]) == {"a", "b"}


def test_tail_col_falls_back_when_readn_em_is_absent():
    # No readn_em column: must rank on readn rather than raise.
    out = normalise(_mixed(), NormaliseOpts(tail_ba=0.01, tail_col="readn_em"))
    assert len(out) < 4


def test_csum_keeps_its_own_preceding_convention():
    # The trim is inclusive, `csum` is preceding. They are different columns
    # answering different questions; unifying them would move published tables.
    out = normalise(_mixed())
    assert out["csum"].iloc[0] == pytest.approx(0.0)
    assert out["csum"].is_monotonic_increasing


# ------------------------------------------------------------------ #
# realtime: additive, off by default
# ------------------------------------------------------------------ #

def _abund():
    return pd.DataFrame({
        "species": ["Target virus", "noise1", "noise2"],
        "readn": [1000, 12, 11],
        "prob": [0.99, 0.005, 0.005],
    })


def _cfg(**kw):
    return RealtimeConfig(fastq_dir=".", library="L",
                          truth_taxon="Target virus", **kw)


def test_realtime_tail_defaults_to_off():
    assert _cfg().tail == 0.0
    # Same answer as before the flag existed: two noise taxa clear read_min=10.
    assert score_snapshot(_abund(), _cfg(read_min=0), truth=None)["FP"] == 2


def test_realtime_tail_removes_the_noise_rows():
    # readn 1000, 12, 11 over 1023: inclusive shares 0.9775, 0.9893, 1.0.
    # tail 0.02 -> threshold 0.98 -> only the true taxon survives.
    got = score_snapshot(_abund(), _cfg(read_min=0, tail=0.02), truth=None)
    assert got["FP"] == 0
    assert got["TP"] == 1


def test_realtime_tail_and_read_min_answer_different_questions():
    # Both noise taxa clear read_min=10, so that filter alone leaves them in.
    assert score_snapshot(_abund(), _cfg(read_min=10), truth=None)["FP"] == 2
    # The tail catches them because they are a rounding error of the library.
    assert score_snapshot(_abund(), _cfg(read_min=10, tail=0.02), truth=None)["FP"] == 0
