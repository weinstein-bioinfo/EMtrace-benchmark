"""Module C metrics — when the answer arrives, and whether the pipeline keeps up.

Three outputs, all computed from the per-batch snapshots collected during a
replay (see :mod:`ngs_benchmark.realtime.replay`):

``timeseries.tsv``
    one row per batch: both clocks, cumulative reads, the current top taxon and
    its EM probability, and — for a profiling run with a multi-taxon truth —
    the L1 distance to the expected composition.

``first_call.tsv``
    the first batch whose top taxon is the truth taxon **and is never displaced
    afterwards**. A transient correct call at batch 1 that later flips is not a
    call; requiring persistence is what makes the number mean "from here on, a
    clinician reading the screen would have been right".

``first_detection.tsv``
    the first batch from which the truth taxon stays at or above the read
    threshold (``read_min``) in every later batch, whatever its rank. On a field
    library the virus is a minor component beside the host microbiome, so it is
    never top-ranked and ``first_call`` never fires; detection is the question
    there, not dominance.

``panel.tsv``, ``panel_summary.tsv``
    field libraries only (``--panel``). A field library has no complete truth —
    the host microbiome and any unlisted virus are unknown, not false — so
    precision, recall and FP are not computed (left empty). Instead each taxon
    the source study reports is scored as detected yes/no (at or above
    ``read_min`` at the end of the run, with its persistent first detection),
    by the kind of evidence behind it: ``pcr_pos`` (independently confirmed;
    summarised as "PCR targets detected k/n"), ``pcr_neg`` (a specificity check:
    the target should stay below ``read_min``) and ``reported`` (found by the
    study's own sequence analysis, e.g. BLAST; agreement, not truth).

``latency.tsv``
    median and p95 per-batch latency, normalised per 1k reads so batches of
    different sizes compare, plus the throughput ratio — sequencing time covered
    divided by compute time. Above 1 the pipeline is faster than the sequencer,
    so an accelerated replay's sequencing-time results would also hold live;
    below 1 a live run would fall behind by a growing backlog.

``accuracy.tsv``
    set metrics at the three moments a reader cares about — the first correct
    call, the stabilisation point, and the end of the run. "When" is only half
    the question: a call that arrives in 90 seconds beside four false positives
    is not the same result as a clean one, and none of the three tables above
    would show the difference.

Per-batch TP/FP/FN/Pr/Rc/F1 also ride along in ``timeseries.tsv``, so precision
and recall can be plotted against the stabilisation curve. On a single-target
viral library these are coarse by construction — with one true taxon, recall is
0 or 1 and precision is 1/(taxa reported) — and they earn their keep on the
multi-taxon mock. ``n_true`` is reported beside them so a reader can see which
regime a row is in.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import pandas as pd

TIMESERIES_COLUMNS = (
    "library", "batch", "seq_time_s", "wall_s", "reads_cum",
    "top_taxon", "em_prob", "L1", "n_truth_taxa", "latency_s",
    "n_pred", "n_true", "TP", "FP", "FN", "Pr", "Rc", "F1",
)
FIRST_CALL_COLUMNS = ("library", "seq_time_s", "wall_s", "reads", "em_prob", "batch")
LATENCY_COLUMNS = (
    "library", "batches", "median_s", "p95_s", "s_per_1k_reads",
    "seq_time_covered_s", "compute_time_s", "ratio",
)
PANEL_COLUMNS = (
    "library", "taxon", "rank", "evidence", "source", "detected", "outcome",
    "reads_final", "first_detection_batch", "first_detection_seq_time_s",
    "first_detection_wall_s", "first_detection_reads_cum",
)
PANEL_SUMMARY_COLUMNS = (
    "library", "pcr_pos_detected", "pcr_pos_n", "pcr_neg_clear", "pcr_neg_n",
    "reported_detected", "reported_n",
)
#: Panel evidence classes, and the outcome label for (not detected, detected).
PANEL_OUTCOMES = {
    "pcr_pos": ("missed", "detected"),
    "pcr_neg": ("clear", "false_detection"),
    "reported": ("not_found", "found"),
}
ACCURACY_COLUMNS = (
    "library", "event", "batch", "seq_time_s", "wall_s", "reads_cum",
    "top_taxon", "em_prob", "L1", "n_pred", "n_true", "TP", "FP", "FN",
    "Pr", "Rc", "F1",
)


@dataclass
class Snapshot:
    """What one batch produced, as seen by the scorer."""

    batch: int
    seq_time_s: float | None
    wall_s: float
    reads_cum: int
    latency_s: float
    #: Highest-abundance taxon in the cumulative table after this batch.
    top_taxon: str | None = None
    #: Its EM posterior, where the pipeline reports one.
    em_prob: float | None = None
    #: L1 against the expected composition (profiling runs only).
    L1: float | None = None
    #: Truth taxa recovered so far (profiling runs only).
    n_truth_taxa: int | None = None
    #: Taxa reported above the abundance threshold after this batch.
    n_pred: int | None = None
    #: Taxa in the truth set.
    n_true: int | None = None
    TP: int | None = None
    FP: int | None = None
    FN: int | None = None
    Pr: float | None = None
    Rc: float | None = None
    F1: float | None = None
    #: Reads per panel taxon after this batch (field libraries only).
    panel_reads: dict = field(default_factory=dict)
    extra: dict = field(default_factory=dict)


def _safe_div(num: float, den: float) -> float:
    return num / den if den else 0.0


def set_scores(pred_taxa: set[str], true_taxa: set[str]) -> dict:
    """TP/FP/FN/Pr/Rc/F1 over two taxon sets.

    Deliberately the same definitions as
    :func:`ngs_benchmark.detection.metrics.set_metrics`, so a Framework C
    accuracy number means what the Framework A one means — a test pins the two
    together. It is reimplemented here rather than called because `set_metrics`
    requires a ``division`` column and taxgroups, which a single-target viral
    run (``--truth-taxon``, no truth table) does not have.
    """
    tp = len(pred_taxa & true_taxa)
    fp = len(pred_taxa - true_taxa)
    fn = len(true_taxa - pred_taxa)
    pr = _safe_div(tp, tp + fp)
    rc = _safe_div(tp, tp + fn)
    return {
        "n_pred": len(pred_taxa),
        "n_true": len(true_taxa),
        "TP": tp, "FP": fp, "FN": fn,
        "Pr": round(pr, 4), "Rc": round(rc, 4),
        "F1": round(_safe_div(2 * pr * rc, pr + rc), 4),
    }


def timeseries(snapshots: list[Snapshot], library: str) -> pd.DataFrame:
    """One row per batch, in replay order."""
    rows = []
    for s in snapshots:
        d = asdict(s)
        d.pop("extra", None)
        rows.append({"library": library, **d})
    return pd.DataFrame(rows, columns=list(TIMESERIES_COLUMNS))


def first_call(
    snapshots: list[Snapshot], truth_taxon: str, library: str
) -> pd.DataFrame:
    """The first correct call that is never displaced afterwards.

    Scans backwards for the last batch whose top taxon is *not* the truth taxon;
    the call is the batch after it. Batches with no top taxon (nothing mapped
    yet) count as incorrect, so an empty prefix delays the call rather than
    being skipped over.

    Returns a single row, or an empty frame when the run never settles on the
    truth taxon.
    """
    if not snapshots:
        return pd.DataFrame(columns=list(FIRST_CALL_COLUMNS))

    last_wrong = -1
    for i, s in enumerate(snapshots):
        if s.top_taxon != truth_taxon:
            last_wrong = i
    if last_wrong == len(snapshots) - 1:          # still wrong at the end
        return pd.DataFrame(columns=list(FIRST_CALL_COLUMNS))

    s = snapshots[last_wrong + 1]
    return pd.DataFrame(
        [{
            "library": library,
            "seq_time_s": s.seq_time_s,
            "wall_s": s.wall_s,
            "reads": s.reads_cum,
            "em_prob": s.em_prob,
            "batch": s.batch,
        }],
        columns=list(FIRST_CALL_COLUMNS),
    )


def first_detection(snapshots: list[Snapshot], library: str) -> pd.DataFrame:
    """The first batch from which the truth taxon is never lost again.

    A batch detects the truth taxon when it is in the predicted set, i.e. at or
    above ``read_min`` reads at the scored rank (``TP`` > 0); rank does not
    matter. As for :func:`first_call`, persistence is required: the scan finds
    the last batch that does *not* detect it, and detection is the batch after.
    Batches never scored (``TP`` is None) count as not detected.

    Returns a single row, or an empty frame when the last batch does not detect
    the truth taxon.
    """
    s = _persistent_from(snapshots, lambda s: bool(s.TP)) if snapshots else None
    if s is None:
        return pd.DataFrame(columns=list(FIRST_CALL_COLUMNS))

    return pd.DataFrame(
        [{
            "library": library,
            "seq_time_s": s.seq_time_s,
            "wall_s": s.wall_s,
            "reads": s.reads_cum,
            "em_prob": s.em_prob,
            "batch": s.batch,
        }],
        columns=list(FIRST_CALL_COLUMNS),
    )


def _persistent_from(snapshots: list[Snapshot], hit) -> Snapshot | None:
    """First snapshot from which ``hit(s)`` holds in every later batch."""
    last_miss = -1
    for i, s in enumerate(snapshots):
        if not hit(s):
            last_miss = i
    if last_miss == len(snapshots) - 1:
        return None
    return snapshots[last_miss + 1]


def panel(
    snapshots: list[Snapshot], targets: pd.DataFrame, library: str, read_min: float
) -> pd.DataFrame:
    """Detected yes/no per panel taxon, with its persistent first detection.

    *targets* has columns ``taxon``, ``rank``, ``evidence`` and ``source`` (see
    :func:`ngs_benchmark.realtime.runner.load_panel`). A taxon is detected when
    it has at least *read_min* reads in the last batch; its first detection is
    the batch from which that holds in every later batch, as for
    :func:`first_detection`. For a ``pcr_neg`` taxon a detection is the
    specificity failure, so it is reported the same way.
    """
    rows = []
    for t in targets.itertuples(index=False):
        def hit(s, taxon=t.taxon):
            return s.panel_reads.get(taxon, 0) >= read_min

        first = _persistent_from(snapshots, hit) if snapshots else None
        detected = first is not None
        rows.append({
            "library": library, "taxon": t.taxon, "rank": t.rank,
            "evidence": t.evidence, "source": t.source,
            "detected": "yes" if detected else "no",
            "outcome": PANEL_OUTCOMES[t.evidence][detected],
            "reads_final": snapshots[-1].panel_reads.get(t.taxon, 0) if snapshots else None,
            "first_detection_batch": first.batch if detected else None,
            "first_detection_seq_time_s": first.seq_time_s if detected else None,
            "first_detection_wall_s": first.wall_s if detected else None,
            "first_detection_reads_cum": first.reads_cum if detected else None,
        })
    df = pd.DataFrame(rows, columns=list(PANEL_COLUMNS))
    # Nullable ints: a taxon never detected must not turn batch 3 into 3.0.
    for col in ("first_detection_batch", "first_detection_reads_cum"):
        df[col] = df[col].astype("Int64")
    return df


def panel_summary(panel_table: pd.DataFrame, library: str) -> pd.DataFrame:
    """The panel in three counts: PCR targets detected, PCR negatives clear,
    study-reported taxa found — each as k and n."""
    def count(evidence: str, outcome: str) -> tuple[int, int]:
        sel = panel_table[panel_table["evidence"] == evidence]
        return int((sel["outcome"] == outcome).sum()), len(sel)

    pos_k, pos_n = count("pcr_pos", "detected")
    neg_k, neg_n = count("pcr_neg", "clear")
    rep_k, rep_n = count("reported", "found")
    return pd.DataFrame(
        [{
            "library": library,
            "pcr_pos_detected": pos_k, "pcr_pos_n": pos_n,
            "pcr_neg_clear": neg_k, "pcr_neg_n": neg_n,
            "reported_detected": rep_k, "reported_n": rep_n,
        }],
        columns=list(PANEL_SUMMARY_COLUMNS),
    )


def latency(snapshots: list[Snapshot], library: str) -> pd.DataFrame:
    """Per-batch latency summary and the throughput ratio.

    ``s_per_1k_reads`` uses total compute over total reads rather than a mean of
    per-batch rates, so a library of 13 fat batches and one of 174 thin batches
    are comparable. ``ratio`` is left empty when the library carries no
    sequencing times.
    """
    if not snapshots:
        return pd.DataFrame(columns=list(LATENCY_COLUMNS))

    lat = pd.Series([s.latency_s for s in snapshots], dtype=float)
    compute_s = float(lat.sum())
    reads = snapshots[-1].reads_cum
    seq_times = [s.seq_time_s for s in snapshots if s.seq_time_s is not None]
    seq_covered = max(seq_times) if seq_times else None

    return pd.DataFrame(
        [{
            "library": library,
            "batches": len(snapshots),
            "median_s": round(float(lat.median()), 3),
            "p95_s": round(float(lat.quantile(0.95)), 3),
            "s_per_1k_reads": round(compute_s / reads * 1000, 3) if reads else None,
            "seq_time_covered_s": None if seq_covered is None else round(seq_covered, 1),
            "compute_time_s": round(compute_s, 3),
            "ratio": (
                None if not seq_covered or compute_s == 0
                else round(seq_covered / compute_s, 2)
            ),
        }],
        columns=list(LATENCY_COLUMNS),
    )


def stabilised_at(
    snapshots: list[Snapshot], column: str = "L1", tol: float = 0.02
) -> Snapshot | None:
    """First batch after which *column* stays within *tol* of its final value.

    The stabilisation point the roadmap asks for: the read count beyond which
    further sequencing does not change the answer. Returns ``None`` when the
    curve never settles, or when the column was never populated.
    """
    values = [(s, getattr(s, column)) for s in snapshots]
    values = [(s, v) for s, v in values if v is not None]
    if not values:
        return None

    final = values[-1][1]
    stable_from = None
    for s, v in reversed(values):
        if abs(v - final) > tol:
            break
        stable_from = s
    return stable_from


def accuracy(
    snapshots: list[Snapshot],
    library: str,
    *,
    truth_taxon: str | None = None,
    stabilise_on: str = "L1",
    tol: float = 0.02,
) -> pd.DataFrame:
    """Set metrics at the moments worth quoting: first call, first detection,
    stabilisation, end.

    Answers "how good was the answer *when it arrived*" — a call that lands in
    90 s beside four false positives is not the result a headline number should
    imply. Rows are emitted only for events that occurred, always in run order,
    so a library that never settles simply has fewer rows.

    ``stabilise_on`` picks the curve that defines stabilisation: ``L1`` for a
    profiling run, ``F1`` or ``em_prob`` for a single-target one, where L1 is
    never computed.
    """
    events: list[tuple[str, Snapshot]] = []

    if truth_taxon:
        fc = first_call(snapshots, truth_taxon, library)
        if not fc.empty:
            idx = int(fc.iloc[0]["batch"])
            match = next((s for s in snapshots if s.batch == idx), None)
            if match is not None:
                events.append(("first_call", match))
        fd = first_detection(snapshots, library)
        if not fd.empty:
            idx = int(fd.iloc[0]["batch"])
            match = next((s for s in snapshots if s.batch == idx), None)
            if match is not None:
                events.append(("first_detection", match))

    stable = stabilised_at(snapshots, column=stabilise_on, tol=tol)
    if stable is not None:
        events.append(("stabilised", stable))

    if snapshots:
        events.append(("final", snapshots[-1]))

    rows = []
    for event, s in events:
        rows.append({
            "library": library, "event": event, "batch": s.batch,
            "seq_time_s": s.seq_time_s, "wall_s": s.wall_s, "reads_cum": s.reads_cum,
            "top_taxon": s.top_taxon, "em_prob": s.em_prob, "L1": s.L1,
            "n_pred": s.n_pred, "n_true": s.n_true,
            "TP": s.TP, "FP": s.FP, "FN": s.FN,
            "Pr": s.Pr, "Rc": s.Rc, "F1": s.F1,
        })
    return pd.DataFrame(rows, columns=list(ACCURACY_COLUMNS))
