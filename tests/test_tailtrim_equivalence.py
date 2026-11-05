"""The tail trim is defined once and implemented twice — prove they agree.

`ngs_benchmark` and `emtrace` are separate packages and cannot import each
other, so `tailtrim.py` is duplicated verbatim in both. This file is the only
mechanical link between the copies: it fails the moment they diverge, whether
by an edit to one that was not copied to the other, or by a behavioural drift
that a byte comparison would miss.

Skipped when `emtrace` is not importable, so the benchmark suite still runs
on its own.
"""

from __future__ import annotations

import inspect

import pandas as pd
import pytest

from ngs_benchmark import tailtrim as bench

emtrace = pytest.importorskip(
    "emtrace.analysis.tailtrim",
    reason="emtrace not installed; cross-repo equivalence not checkable here",
)


# ------------------------------------------------------------------ #
# The source is the contract
# ------------------------------------------------------------------ #

def test_modules_are_byte_identical():
    """The two files must be copies, not merely similar.

    If this fails, copy whichever you edited over the other:
        EMtrace/src/emtrace/analysis/tailtrim.py
        EMtrace-benchmark/ngs-benchmark/ngs_benchmark/tailtrim.py
    """
    assert inspect.getsource(bench) == inspect.getsource(emtrace)


@pytest.mark.parametrize("name", ["tail_trim", "resolve_count_col"])
def test_docstrings_match(name):
    assert getattr(bench, name).__doc__ == getattr(emtrace, name).__doc__


@pytest.mark.parametrize("name", ["tail_trim", "resolve_count_col"])
def test_signatures_match(name):
    assert (inspect.signature(getattr(bench, name))
            == inspect.signature(getattr(emtrace, name)))


# ------------------------------------------------------------------ #
# ...and so is the behaviour
# ------------------------------------------------------------------ #

VI, BA = "Viruses", "Bacteria"


def _frame(counts, divisions=None, readn_em=None):
    data = {"taxid": [f"t{i}" for i in range(len(counts))], "readn": counts}
    if divisions is not None:
        data["division"] = divisions
    if readn_em is not None:
        data["readn_em"] = readn_em
    return pd.DataFrame(data)


# (label, frame, kwargs) -- every edge the two copies must agree on.
CASES = [
    ("long tail",
     _frame([1000, 500, 100, 50, 10, 5, 2, 1]), {"tail": 0.01}),
    ("ties at the boundary",
     _frame([100, 10, 10, 10, 10]), {"tail": 0.05}),
    ("tail zero",
     _frame([1000, 100, 10, 1]), {"tail": 0.0}),
    ("tail one",
     _frame([1000, 100, 10, 1]), {"tail": 1.0}),
    ("single row",
     _frame([42]), {"tail": 0.01}),
    ("empty frame",
     _frame([]), {"tail": 0.01}),
    ("all-zero counts",
     _frame([0, 0, 0]), {"tail": 0.01}),
    ("readn_em all NA",
     _frame([1000, 100, 1], readn_em=["NA", "NA", "NA"]), {"tail": 0.01}),
    ("readn_em partly numeric",
     _frame([1000, 100, 1], readn_em=["NA", 900.0, 2.0]), {"tail": 0.01}),
    ("readn_em ranks differently",
     _frame([1, 500, 1000], readn_em=[1000.0, 500.0, 1.0]), {"tail": 0.01}),
    ("protected viral division",
     _frame([1000, 100, 2, 1], [BA, BA, VI, VI]),
     {"tail": 0.01, "group_col": "division", "group_tails": {VI: 0.0}}),
    ("virus-dominated, viruses protected",
     _frame([10000, 100, 10, 5, 1], [VI, BA, BA, BA, BA]),
     {"tail": 0.01, "group_col": "division", "group_tails": {VI: 0.0}}),
    ("group column absent",
     _frame([1000, 100, 10, 1]),
     {"tail": 0.01, "group_col": "division", "group_tails": {VI: 0.0}}),
]
IDS = [c[0] for c in CASES]


@pytest.mark.parametrize("_, df, kwargs", CASES, ids=IDS)
@pytest.mark.parametrize("preserve_order", [True, False])
def test_tail_trim_agrees(_, df, kwargs, preserve_order):
    a = bench.tail_trim(df, preserve_order=preserve_order, **kwargs)
    b = emtrace.tail_trim(df, preserve_order=preserve_order, **kwargs)
    pd.testing.assert_frame_equal(a, b)


@pytest.mark.parametrize("_, df, kwargs", CASES, ids=IDS)
def test_resolve_count_col_agrees(_, df, kwargs):
    assert bench.resolve_count_col(df) == emtrace.resolve_count_col(df)
