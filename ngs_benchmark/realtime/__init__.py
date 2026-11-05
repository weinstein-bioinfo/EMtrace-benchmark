"""Module C — real-time performance (when the answer arrives, not just what it is).

Replays a MinKNOW ``fastq_pass`` directory through a real-time pipeline in its
original batch files and scores the per-batch snapshots: time to the first
correct call that is never displaced, per-batch latency against the sequencer's
own pace, and how the profile stabilises as reads accumulate.

Public entry points:

    from ngs_benchmark.realtime.replay import scan_library, replay
    from ngs_benchmark.realtime.runner import RealtimeConfig, EMtraceDriver, run, rescore
    from ngs_benchmark.realtime.metrics import first_call, first_detection, latency, timeseries
"""

from .metrics import (
    Snapshot,
    first_call,
    first_detection,
    latency,
    panel,
    panel_summary,
    stabilised_at,
    timeseries,
)
from .replay import Batch, Release, contiguity_gaps, replay, scan_library
from .runner import RealtimeConfig, EMtraceDriver, load_panel, rescore, run, score_snapshot

__all__ = [
    "Batch",
    "Release",
    "RealtimeConfig",
    "EMtraceDriver",
    "Snapshot",
    "contiguity_gaps",
    "first_call",
    "first_detection",
    "latency",
    "load_panel",
    "panel",
    "panel_summary",
    "replay",
    "rescore",
    "run",
    "scan_library",
    "score_snapshot",
    "stabilised_at",
    "timeseries",
]
