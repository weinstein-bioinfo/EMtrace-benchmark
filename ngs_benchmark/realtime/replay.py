"""Replay a MinKNOW ``fastq_pass`` directory batch file by batch file.

The batches a live pipeline would have seen already exist: MinKNOW wrote one
FASTQ per few thousand reads during the run. Replay therefore releases *those*
files, in their original order, rather than re-cutting the library into
synthetic time windows — batch sizes, arrival order and within-batch read mix
stay the ones a real run produced.

Two clocks are tracked per batch and both are reported downstream:

* **sequencing time** — from the ``start_time`` headers inside the file, so it
  comes from the data rather than from file ordering;
* **compute time** — wall clock from release to the pipeline reporting, which
  is what decides whether the pipeline keeps up with the sequencer.

Release is *accelerated*: the next file goes out as soon as the previous one has
been reported, not at the original pace. A faithful replay would take as long as
the original run (2 h to 72 h for the Framework C libraries).
"""

from __future__ import annotations

import gzip
import re
import shutil
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterator, Sequence

#: MinKNOW batch files end in ``_<n>.fastq[.gz]``; *n* is the ordering key.
_BATCH_SUFFIX = re.compile(r"_(\d+)\.fastq(?:\.gz)?$")
_START_TIME = re.compile(r"start_time=(\S+)")


@dataclass(frozen=True)
class Batch:
    """One MinKNOW batch file, with what replay needs to know about it."""

    #: Position in the replay, 1-based (batch 1 is the first file released).
    index: int
    #: The file MinKNOW wrote.
    path: Path
    #: Reads in this file.
    reads: int
    #: Latest read ``start_time`` in the file, as seconds since the run's first
    #: read. ``None`` when the file carries no timestamps.
    seq_time_s: float | None


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def batch_files(fastq_dir: str | Path) -> list[Path]:
    """Return a library's batch files in replay order (``_0``, ``_1``, ``_2``, …).

    Ordering is by the numeric suffix, not lexicographic: ``_10`` sorts after
    ``_9``, which a plain ``sorted()`` would get wrong.
    """
    paths = [p for p in Path(fastq_dir).iterdir() if _BATCH_SUFFIX.search(p.name)]
    if not paths:
        raise ValueError(
            f"No MinKNOW batch files (*_<n>.fastq[.gz]) in {fastq_dir}"
        )
    return sorted(paths, key=lambda p: int(_BATCH_SUFFIX.search(p.name).group(1)))


def scan_batch(path: Path) -> tuple[int, datetime | None, datetime | None]:
    """Count reads in *path* and return its earliest and latest ``start_time``.

    Both ends are needed: the earliest across the library fixes the run's t0,
    while each batch is timed by its *latest* read, since a batch is only
    complete once its last read has been sequenced — that is the moment a live
    pipeline could first have seen the whole file.
    """
    reads, earliest, latest = 0, None, None
    with _open_text(path) as fh:
        for i, line in enumerate(fh):
            if i % 4:
                continue
            reads += 1
            m = _START_TIME.search(line)
            if m:
                ts = datetime.fromisoformat(m.group(1))
                if earliest is None or ts < earliest:
                    earliest = ts
                if latest is None or ts > latest:
                    latest = ts
    return reads, earliest, latest


def scan_library(fastq_dir: str | Path, reads_only: bool = False) -> list[Batch]:
    """Scan every batch file once, returning them in replay order.

    Sequencing times are expressed relative to the library's earliest read, so
    batch 1 starts near 0 s regardless of when the flow cell ran.

    A library without ``start_time`` headers is refused unless *reads_only* is
    set: SRA deposits have them stripped, and such a library can only be
    replayed on read count. Every batch then has ``seq_time_s=None``.
    """
    paths = batch_files(fastq_dir)
    scanned = [(p, *scan_batch(p)) for p in paths]
    firsts = [ts for _, _, ts, _ in scanned if ts is not None]
    if reads_only:
        if firsts:
            raise ValueError(
                f"{fastq_dir} carries start_time headers — replay it on sequencing "
                "time, not with reads_only"
            )
        return [
            Batch(index=i, path=p, reads=reads, seq_time_s=None)
            for i, (p, reads, _first, _last) in enumerate(scanned, start=1)
        ]
    if not firsts:
        raise ValueError(
            f"No start_time headers in {fastq_dir} — this library cannot be "
            "replayed on sequencing time (SRA strips them); use reads_only"
        )
    t0 = min(firsts)
    return [
        Batch(
            index=i,
            path=p,
            reads=reads,
            seq_time_s=None if last is None else (last - t0).total_seconds(),
        )
        for i, (p, reads, _first, last) in enumerate(scanned, start=1)
    ]


def contiguity_gaps(batches: Sequence[Batch], max_gap_h: float = 24.0) -> list[tuple[int, float]]:
    """Return ``(batch index, gap in hours)`` for jumps larger than *max_gap_h*.

    A combined library — two flow-cell runs merged into one directory — replays
    as a single run with an impossible gap in the middle. Callers should refuse such a
    library rather than report its timings.
    """
    gaps: list[tuple[int, float]] = []
    prev = None
    for b in batches:
        if b.seq_time_s is None:
            continue
        if prev is not None:
            gap_h = (b.seq_time_s - prev) / 3600.0
            if gap_h > max_gap_h:
                gaps.append((b.index, gap_h))
        prev = b.seq_time_s
    return gaps


@dataclass
class Release:
    """One released batch and the wall-clock cost of processing it."""

    batch: Batch
    #: Cumulative reads released, including this batch.
    reads_cum: int
    #: Seconds from release to the pipeline reporting this batch.
    latency_s: float
    #: Seconds since the replay started, at the moment this batch was reported.
    wall_s: float


def replay(
    batches: Sequence[Batch],
    watch_dir: str | Path,
    await_batch: Callable[[int], None],
    *,
    copy: bool = True,
) -> Iterator[Release]:
    """Release batches one at a time, yielding a :class:`Release` for each.

    *await_batch* is called with the 1-based batch index after the file lands in
    *watch_dir*, and must block until the pipeline has finished reporting that
    batch; the time it blocks is the batch's latency. Passing ``copy=False``
    hard-links instead of copying, which avoids duplicating a 300 MB library —
    use it only when *watch_dir* is on the same filesystem as the source.

    The generator form matters: a caller can write each row as it arrives, so a
    run interrupted after six hours still leaves five hours of usable timings.
    """
    watch = Path(watch_dir)
    watch.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    reads_cum = 0

    for b in batches:
        dest = watch / b.path.name
        released = time.monotonic()
        if copy:
            shutil.copy2(b.path, dest)
        else:
            dest.unlink(missing_ok=True)
            dest.hardlink_to(b.path)
        await_batch(b.index)
        now = time.monotonic()
        reads_cum += b.reads
        yield Release(
            batch=b,
            reads_cum=reads_cum,
            latency_s=now - released,
            wall_s=now - started,
        )
