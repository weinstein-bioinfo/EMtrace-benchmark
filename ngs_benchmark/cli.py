"""ngs-benchmark CLI entry point.

Beta slice exposes:
    ngs-benchmark detection         — single pred × truth
    ngs-benchmark detection-batch   — many preds × one truth (manifest TSV)
    ngs-benchmark assembly          — single consensus × reference (Module B)
    ngs-benchmark assembly-batch    — many consensus/reference pairs (manifest TSV)
    ngs-benchmark bgfilter          — single removed-reads FASTQ (Module D)
    ngs-benchmark bgfilter-batch    — many filter settings (simmet truth files)
    ngs-benchmark bgfilter-divbatch — many filter settings (per-read division truth)
"""

from __future__ import annotations

import warnings

import click
import pandas as pd

from .assembly.runner import AssemblyConfig
from .assembly.runner import run as assembly_run
from .assembly.runner import run_batch as assembly_run_batch
from .bgfilter.runner import BgfilterConfig
from .bgfilter.runner import run as bgfilter_run
from .bgfilter.runner import run_batch as bgfilter_run_batch
from .bgfilter.runner import run_division_batch as bgfilter_run_division_batch
from .bgfilter.truth import count_data_lines, load_division_truth, load_truth
from .detection.runner import DetectionConfig, run, run_batch
from .realtime.runner import RealtimeConfig, EMtraceDriver
from .realtime.runner import rescore as realtime_rescore
from .realtime.runner import run as realtime_run


def _split_csv(value: str) -> tuple[str, ...]:
    return tuple(s.strip() for s in value.split(",") if s.strip())


def _resolve_tails(
    tail: float, tail_vi: float | None, tail_ba: float | None,
) -> tuple[float, float]:
    """Resolve the viral and non-viral tail thresholds, as FRACTIONS.

    `--tail` is the fraction applied to every division. The older
    `--tail-vi` / `--tail-ba` still work and still mean PERCENT, so scripts
    written against them keep their meaning instead of silently asking for a
    100x larger trim; each warns and is divided by 100.

    The conversion lives here rather than in NormaliseOpts so the config
    carries one unit system -- which is what `summary.json` provenance
    records -- and so the warning fires where the deprecated spelling was
    typed, not on every internal NormaliseOpts construction.
    """
    resolved = []
    for name, percent in (("--tail-vi", tail_vi), ("--tail-ba", tail_ba)):
        if percent is None:
            resolved.append(tail)
            continue
        fraction = percent / 100.0
        message = (
            f"{name} is deprecated and is read as PERCENT ({percent:g} -> "
            f"tail={fraction:g}). Use --tail, which takes a fraction."
        )
        warnings.warn(message, DeprecationWarning, stacklevel=2)
        click.echo(f"warning: {message}", err=True)
        resolved.append(fraction)
    return resolved[0], resolved[1]


@click.group()
@click.version_option()
def cli() -> None:
    """Benchmarking metrics for Nanopore mNGS pipelines."""


@cli.command()
@click.option("--pred", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--truth", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--outdir", required=True, type=click.Path(file_okay=False))
@click.option("--pred-fmt", default="emtrace", show_default=True,
              help="Adapter name. Beta slice supports: emtrace.")
@click.option("--ranks", default="species,genus,family", show_default=True,
              callback=lambda c, p, v: _split_csv(v))
@click.option("--taxgroup", "taxgroups", default="Viruses", show_default=True,
              callback=lambda c, p, v: _split_csv(v),
              help="Top-level group to restrict metrics to (comma-separated for several).")
@click.option("--score-col", default="readn", show_default=True,
              help="Column used for the PR threshold sweep. Use 'prob' or 'pvalue' for EMtrace EM output.")
@click.option("--read-min", type=float, default=0.0, show_default=True)
@click.option("--min-rel", type=float, default=0.0, show_default=True,
              help="Per-taxon relative-abundance floor, as a FRACTION "
                   "(0.0001 = 0.01%) of the --filter-col total left after "
                   "--read-min.")
@click.option("--min-rel-exclude", default="", show_default=True,
              callback=lambda c, p, v: _split_csv(v) if v else (),
              help="Divisions left out of the --min-rel denominator, comma-separated "
                   "(e.g. Primates, for host reads a pipeline reports rather than filters).")
@click.option("--filter-col", type=click.Choice(["readn", "readn_em"]),
              default="readn", show_default=True,
              help="Count column --read-min and --min-rel act on. readn_em "
                   "falls back to readn when absent or non-numeric.")
@click.option("--tail", type=float, default=0.0, show_default=True,
              help="Cumulative-share tail trim, as a FRACTION (0.01 = 1%). "
                   "Applies to every division unless overridden by "
                   "--tail-vi / --tail-ba.")
@click.option("--tail-col", type=click.Choice(["readn", "readn_em"]),
              default="readn", show_default=True,
              help="Count column the tail trim ranks on. readn_em falls back "
                   "to readn when absent or non-numeric (EM did not run). "
                   "Does not affect --read-min, which acts on --filter-col.")
@click.option("--tail-vi", type=float, default=None,
              help="DEPRECATED, PERCENT units (5 = 5%). Overrides --tail for "
                   "Viruses. Use --tail with a fraction instead.")
@click.option("--tail-ba", type=float, default=None,
              help="DEPRECATED, PERCENT units (1 = 1%). Overrides --tail for "
                   "non-Viruses. Use --tail with a fraction instead.")
@click.option("--exclude-retro/--no-exclude-retro", default=False)
@click.option("--drop-unknown/--keep-unknown", default=False,
              help='Do not count "unknown" (reads placed above a rank) as a taxon call at that rank.')
@click.option("--keep-phages/--drop-phages", default=False,
              help="Keep phages in both pred and truth. Default is to drop them.")
@click.option("--buddle/--no-buddle", default=False,
              help="Enable Buddle PMR/RPMR threshold (§3.5). Requires --nc-pred, --total-reads-sample, --total-reads-nc.")
@click.option("--nc-pred", default=None, type=click.Path(exists=True, dir_okay=False),
              help="Negative-control prediction file (same format as --pred).")
@click.option("--total-reads-sample", default=None, type=float)
@click.option("--total-reads-nc", default=None, type=float)
@click.option("--tn-universe", default=None, type=click.Path(exists=True, dir_okay=False),
              help="(Deferred from beta.) Curated TN-universe file enabling §3.4 ROC.")
def detection(
    pred: str,
    truth: str,
    outdir: str,
    pred_fmt: str,
    ranks: tuple[str, ...],
    taxgroups: tuple[str, ...],
    score_col: str,
    read_min: float,
    min_rel: float,
    min_rel_exclude: tuple[str, ...],
    filter_col: str,
    tail: float,
    tail_col: str,
    tail_vi: float | None,
    tail_ba: float | None,
    exclude_retro: bool,
    keep_phages: bool,
    drop_unknown: bool,
    buddle: bool,
    nc_pred: str | None,
    total_reads_sample: float | None,
    total_reads_nc: float | None,
    tn_universe: str | None,
) -> None:
    """Score one pipeline prediction against a truth table."""
    _tails = _resolve_tails(tail, tail_vi, tail_ba)
    cfg = DetectionConfig(
        pred=pred,
        truth=truth,
        outdir=outdir,
        pred_fmt=pred_fmt,
        ranks=ranks,
        taxgroups=taxgroups,
        score_col=score_col,
        read_min=read_min,
        min_rel=min_rel,
        min_rel_exclude=min_rel_exclude,
        filter_col=filter_col,
        tail_vi=_tails[0],
        tail_ba=_tails[1],
        tail_col=tail_col,
        exclude_retro=exclude_retro,
        keep_phages=keep_phages,
        drop_unknown=drop_unknown,
        buddle=buddle,
        nc_pred=nc_pred,
        total_reads_sample=total_reads_sample,
        total_reads_nc=total_reads_nc,
        tn_universe=tn_universe,
    )
    summary = run(cfg)
    click.echo(f"Wrote metrics + curves to {outdir}/. AUPR={summary['AUPR']:.4f} F_max={summary['F_max']:.4f}")


@cli.command(name="detection-batch")
@click.option("--truth", required=True, type=click.Path(exists=True, dir_okay=False))
@click.option("--pred-list", required=True, type=click.Path(exists=True, dir_okay=False),
              help="TSV manifest: cols pipeline, sample_id, path, fmt, [gc_per_ml], [matrix_id], ...")
@click.option("--outdir", required=True, type=click.Path(file_okay=False))
@click.option("--pred-fmt", default="emtrace", show_default=True)
@click.option("--ranks", default="species,genus,family", show_default=True,
              callback=lambda c, p, v: _split_csv(v))
@click.option("--taxgroup", "taxgroups", default="Viruses", show_default=True,
              callback=lambda c, p, v: _split_csv(v))
@click.option("--score-col", default="readn", show_default=True)
@click.option("--read-min", type=float, default=0.0, show_default=True)
@click.option("--min-rel", type=float, default=0.0, show_default=True,
              help="Per-taxon relative-abundance floor, as a FRACTION "
                   "(0.0001 = 0.01%) of the --filter-col total left after "
                   "--read-min.")
@click.option("--min-rel-exclude", default="", show_default=True,
              callback=lambda c, p, v: _split_csv(v) if v else (),
              help="Divisions left out of the --min-rel denominator, comma-separated "
                   "(e.g. Primates, for host reads a pipeline reports rather than filters).")
@click.option("--filter-col", type=click.Choice(["readn", "readn_em"]),
              default="readn", show_default=True,
              help="Count column --read-min and --min-rel act on. readn_em "
                   "falls back to readn when absent or non-numeric.")
@click.option("--tail", type=float, default=0.0, show_default=True,
              help="Cumulative-share tail trim, as a FRACTION (0.01 = 1%). "
                   "Applies to every division unless overridden by "
                   "--tail-vi / --tail-ba.")
@click.option("--tail-col", type=click.Choice(["readn", "readn_em"]),
              default="readn", show_default=True,
              help="Count column the tail trim ranks on. readn_em falls back "
                   "to readn when absent or non-numeric (EM did not run). "
                   "Does not affect --read-min, which acts on --filter-col.")
@click.option("--tail-vi", type=float, default=None,
              help="DEPRECATED, PERCENT units. Overrides --tail for Viruses.")
@click.option("--tail-ba", type=float, default=None,
              help="DEPRECATED, PERCENT units. Overrides --tail for non-Viruses.")
@click.option("--exclude-retro/--no-exclude-retro", default=False)
@click.option("--drop-unknown/--keep-unknown", default=False,
              help='Do not count "unknown" (reads placed above a rank) as a taxon call at that rank.')
@click.option("--keep-phages/--drop-phages", default=False)
@click.option("--dilution-series/--no-dilution-series", default=False,
              help="Emit sensitivity_vs_gcml.tsv (requires manifest with gc_per_ml column).")
def detection_batch(
    truth: str,
    pred_list: str,
    outdir: str,
    pred_fmt: str,
    ranks: tuple[str, ...],
    taxgroups: tuple[str, ...],
    score_col: str,
    read_min: float,
    min_rel: float,
    min_rel_exclude: tuple[str, ...],
    filter_col: str,
    tail: float,
    tail_col: str,
    tail_vi: float | None,
    tail_ba: float | None,
    exclude_retro: bool,
    keep_phages: bool,
    drop_unknown: bool,
    dilution_series: bool,
) -> None:
    """Score many predictions against one truth."""
    manifest = pd.read_csv(pred_list, sep="\t")
    _tails = _resolve_tails(tail, tail_vi, tail_ba)
    cfg = DetectionConfig(
        pred="",  # unused in batch mode
        truth=truth,
        outdir=outdir,
        pred_fmt=pred_fmt,
        ranks=ranks,
        taxgroups=taxgroups,
        score_col=score_col,
        read_min=read_min,
        min_rel=min_rel,
        min_rel_exclude=min_rel_exclude,
        filter_col=filter_col,
        tail_vi=_tails[0],
        tail_ba=_tails[1],
        tail_col=tail_col,
        exclude_retro=exclude_retro,
        keep_phages=keep_phages,
        drop_unknown=drop_unknown,
        dilution_series=dilution_series,
    )
    summary = run_batch(cfg, manifest, truth_path=truth, outdir=outdir)
    click.echo(f"Wrote {len(manifest)} per-sample outputs + combined_metrics.tsv to {outdir}/.")
    if dilution_series:
        click.echo(f"Dilution series total_unique_fps={summary.get('total_unique_fps', 'NA')}")


@cli.command()
@click.option("--consensus", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Consensus FASTA (may be multi-segment: one record per segment).")
@click.option("--reference", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Reference FASTA (one record per segment).")
@click.option("--outdir", required=True, type=click.Path(file_okay=False))
@click.option("--sample", default=None, help="Sample label added to summary.tsv.")
@click.option("--stats-div", default=None, type=click.Path(dir_okay=False),
              help="EMtrace results/reports/stats_div.csv. read_count for "
                   "division=Viruses → viral_reads; missing file → 'no hits'.")
def assembly(
    consensus: str,
    reference: str,
    outdir: str,
    sample: str | None,
    stats_div: str | None,
) -> None:
    """Score one consensus assembly against a reference (Module B)."""
    cfg = AssemblyConfig(
        consensus=consensus,
        reference=reference,
        outdir=outdir,
        sample=sample,
        stats_div=stats_div,
    )
    summary = assembly_run(cfg)
    cov = summary.get("weighted_coverage_pct")
    ide = summary.get("weighted_identity_pct")
    pos = summary.get("weighted_positions_called_pct")
    cov_s = "NA" if cov is None else f"{cov:.2f}"
    ide_s = "NA" if ide is None else f"{ide:.2f}"
    pos_s = "NA" if pos is None else f"{pos:.2f}"
    click.echo(f"Wrote coverage.tsv + summary.tsv to {outdir}/. "
               f"coverage={cov_s}% identity={ide_s}% positions_called={pos_s}%")


@cli.command(name="assembly-batch")
@click.option("--manifest", required=True, type=click.Path(exists=True, dir_okay=False),
              help="TSV manifest: cols sample_id, consensus, reference, [stats_div].")
@click.option("--outdir", required=True, type=click.Path(file_okay=False))
def assembly_batch(manifest: str, outdir: str) -> None:
    """Score many consensus/reference pairs against their references."""
    man = pd.read_csv(manifest, sep="\t")
    required = {"sample_id", "consensus", "reference"}
    missing = required - set(man.columns)
    if missing:
        raise click.ClickException(
            f"manifest is missing required column(s): {', '.join(sorted(missing))}"
        )
    combined = assembly_run_batch(man, outdir=outdir)
    click.echo(f"Scored {len(combined)} samples → {outdir}/summary.tsv + {outdir}/coverage.tsv")


def _truth_options(fn):
    """Shared truth-source options for the bgfilter commands."""
    fn = click.option("--truth-nohuman", required=True, type=click.Path(exists=True, dir_okay=False),
                      help="readid_taxpath TSV of all non-human reads.")(fn)
    fn = click.option("--truth-vi", required=True, type=click.Path(exists=True, dir_okay=False),
                      help="readid_taxpath TSV of viral reads (with species_id column).")(fn)
    fn = click.option("--truth-vi-nophix-noretro", "truth_vinpnr", required=True,
                      type=click.Path(exists=True, dir_okay=False),
                      help="readid_taxpath TSV of viral reads excl. phiX and Retroviridae.")(fn)
    fn = click.option("--total", type=int, default=None,
                      help="Total simmet read count. Required unless --truth-all is given.")(fn)
    fn = click.option("--truth-all", default=None, type=click.Path(exists=True, dir_okay=False),
                      help="Master readid_taxpath TSV; data-line count used as --total.")(fn)
    return fn


def _resolve_truth(total, truth_all, truth_nohuman, truth_vi, truth_vinpnr):
    if total is None:
        if truth_all is None:
            raise click.ClickException("provide --total or --truth-all to set the read total")
        total = count_data_lines(truth_all)
    return load_truth(
        total=total,
        nohuman=truth_nohuman,
        vi=truth_vi,
        vi_nophix_noretro=truth_vinpnr,
    )


@cli.command()
@click.option("--removed-fastq", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Reads removed by the host filter (emtrace hflt.<host>.reads.fq.gz).")
@click.option("--label", required=True, help="Filter label for the table row (e.g. AS100).")
@click.option("--outdir", required=True, type=click.Path(file_okay=False))
@click.option("--strip-mate/--no-strip-mate", default=False, show_default=True,
              help="Drop a trailing /1 or /2 from read ids (paired reads vs per-fragment truth).")
@_truth_options
def bgfilter(removed_fastq, label, outdir, strip_mate, truth_all, total,
             truth_nohuman, truth_vi, truth_vinpnr):
    """Score one host-filter setting: fraction of reads removed per category (Module D)."""
    truth = _resolve_truth(total, truth_all, truth_nohuman, truth_vi, truth_vinpnr)
    row = bgfilter_run(BgfilterConfig(removed_fastq=removed_fastq, label=label), truth, outdir,
                       strip_mate=strip_mate)
    cells = " ".join(
        f"{k}={'NA' if row[k] is None else f'{row[k]:.3f}'}%"
        for k in ("homo_sapiens", "bacteria", "viruses", "viruses_excl_retro")
    )
    click.echo(f"Wrote table.tsv to {outdir}/. {cells}")


@cli.command(name="bgfilter-batch")
@click.option("--manifest", required=True, type=click.Path(exists=True, dir_okay=False),
              help="TSV manifest: cols filter, removed_fastq.")
@click.option("--outdir", required=True, type=click.Path(file_okay=False))
@click.option("--strip-mate/--no-strip-mate", default=False, show_default=True,
              help="Drop a trailing /1 or /2 from read ids (paired reads vs per-fragment truth).")
@_truth_options
def bgfilter_batch(manifest, outdir, strip_mate, truth_all, total,
                   truth_nohuman, truth_vi, truth_vinpnr):
    """Score many host-filter settings from a manifest (one combined table.tsv)."""
    man = pd.read_csv(manifest, sep="\t")
    required = {"filter", "removed_fastq"}
    missing = required - set(man.columns)
    if missing:
        raise click.ClickException(
            f"manifest is missing required column(s): {', '.join(sorted(missing))}"
        )
    truth = _resolve_truth(total, truth_all, truth_nohuman, truth_vi, truth_vinpnr)
    df = bgfilter_run_batch(man, truth, outdir=outdir, strip_mate=strip_mate)
    click.echo(f"Scored {len(df)} filter settings → {outdir}/table.tsv")


@cli.command(name="bgfilter-divbatch")
@click.option("--manifest", required=True, type=click.Path(exists=True, dir_okay=False),
              help="TSV manifest: cols filter, removed_fastq.")
@click.option("--truth-table", required=True, type=click.Path(exists=True, dir_okay=False),
              help="Per-read division truth TSV (cols readid, division, …).")
@click.option("--outdir", required=True, type=click.Path(file_okay=False))
@click.option("--strip-mate/--no-strip-mate", default=False, show_default=True,
              help="Drop a trailing /1 or /2 from read ids before matching truth.")
def bgfilter_divbatch(manifest, truth_table, outdir, strip_mate):
    """Score filter settings against a per-read division truth table (Module D).

    For the simulated-Nanopore datasets: one row per setting, columns = fraction
    of reads removed per division (Viruses / Bacteria / Eukaryota). Use to pick
    the -s threshold that removes the most host with the least viral loss.
    """
    man = pd.read_csv(manifest, sep="\t")
    missing = {"filter", "removed_fastq"} - set(man.columns)
    if missing:
        raise click.ClickException(
            f"manifest is missing required column(s): {', '.join(sorted(missing))}"
        )
    truth = load_division_truth(truth_table)
    df = bgfilter_run_division_batch(man, truth, outdir=outdir, strip_mate=strip_mate)
    click.echo(f"Scored {len(df)} filter settings → {outdir}/table.tsv "
               f"(divisions: {', '.join(truth.divisions)})")


@cli.command()
@click.option("--fastq-dir", required=True, type=click.Path(exists=True, file_okay=False),
              help="MinKNOW fastq_pass directory; per-read start_time headers required "
                   "unless --reads-only.")
@click.option("--outdir", required=True, type=click.Path(file_okay=False))
@click.option("--library", required=True, help="Name used in every output table.")
@click.option("--pipeline", type=click.Choice(["emtrace"]), default="emtrace",
              show_default=True)
@click.option("--truth-taxon", default=None,
              help="Expected taxon — enables first_call.tsv and first_detection.tsv.")
@click.option("--truth", "truth_table", default=None,
              type=click.Path(exists=True, dir_okay=False),
              help="Multi-taxon truth table — enables the L1 column (profiling runs).")
@click.option("--panel", default=None, type=click.Path(exists=True, dir_okay=False),
              help="Field-library panel TSV (taxon, evidence, source[, rank]): "
                   "taxa scored detected yes/no, no Pr/Rc. Excludes --truth-taxon/--truth.")
@click.option("--rank", default="species", show_default=True)
@click.option("--score-col", default="readn", show_default=True)
@click.option("--prob-col", default="prob", show_default=True)
@click.option("--read-min", default=10.0, show_default=True,
              help="Taxa below this many reads are ignored when scoring TP/FP/FN, "
                   "and count as not detected in a panel.")
@click.option("--tail", type=float, default=0.0, show_default=True,
              help="Cumulative-share tail trim, as a FRACTION (0.01 = 1%), "
                   "applied as a row filter before per-rank scoring. "
                   "0 = off, which is what the published replays used.")
@click.option("--tail-col", type=click.Choice(["readn", "readn_em"]),
              default="readn", show_default=True,
              help="Count column the tail trim ranks on; falls back to "
                   "--score-col. Does not affect --read-min.")
@click.option("--stabilise-on", default="L1", show_default=True,
              help="Curve defining stabilisation: L1 for profiling, F1 otherwise.")
@click.option("--stabilise-tol", default=0.02, show_default=True)
@click.option("--max-gap-h", default=24.0, show_default=True,
              help="Refuse a library whose sequencing time jumps more than this.")
@click.option("--hardlink", is_flag=True, default=False,
              help="Hard-link batches into the watch dir instead of copying.")
@click.option("--reads-only", is_flag=True, default=False,
              help="Library has no start_time headers (SRA strips them): replay "
                   "on read count only; seq_time_s and the ratio stay empty.")
@click.option("--sample", default="replay", show_default=True,
              help="Sample name passed to the pipeline.")
@click.option("--batch-timeout", default=3600.0, show_default=True,
              help="Seconds to wait for one batch before giving up.")
@click.option("--config", "pipeline_config", default=None,
              type=click.Path(exists=True, dir_okay=False),
              help="Pipeline config passed to `emtrace watch --config` "
                   "(e.g. the 16S profile).")
@click.option("--anns", default=None,
              help="Annotation strategy passed to `emtrace watch --anns`.")
@click.option("--threads", default=None, type=int,
              help="Threads passed to `emtrace watch -t`.")
def realtime(fastq_dir, outdir, library, pipeline, truth_taxon, truth_table, panel, rank,
             score_col, prob_col, read_min, tail, tail_col, stabilise_on, stabilise_tol,
             max_gap_h, hardlink, reads_only, sample, batch_timeout, pipeline_config, anns,
             threads) -> None:
    """Replay a MinKNOW run through a real-time pipeline and score it (Module C).

    Releases the library's own batch files one at a time — accelerated, so the
    next goes out as soon as the previous is reported — and writes:

    \b
      timeseries.tsv  per batch: both clocks, reads, top taxon, L1, Pr/Rc/F1
      first_call.tsv  first correct call that is never displaced (--truth-taxon)
      first_detection.tsv  first batch from which the truth taxon stays at or
                      above --read-min, whatever its rank (--truth-taxon)
      latency.tsv     median/p95 latency and the throughput ratio
      accuracy.tsv    TP/FP/FN/Pr/Rc/F1 at first call, first detection,
                      stabilisation and the end
      panel.tsv       per panel taxon: detected yes/no, first detection (--panel)
      panel_summary.tsv  PCR targets detected k/n, PCR negatives clear k/n,
                      study-reported taxa found k/n (--panel)
      abund/          every batch's cumulative abundance table (for rescoring)

    A ratio above 1 means the pipeline is faster than the sequencer, so the
    sequencing-time results would also hold in a live run. Accuracy answers the
    other half: how good the answer was at the moment it arrived.
    """
    cfg = RealtimeConfig(
        fastq_dir=fastq_dir, library=library, truth_taxon=truth_taxon,
        truth_table=truth_table, panel=panel, rank=rank, score_col=score_col, prob_col=prob_col,
        read_min=read_min, tail=tail, tail_col=tail_col,
        stabilise_on=stabilise_on, stabilise_tol=stabilise_tol,
        max_gap_h=max_gap_h, hardlink=hardlink, reads_only=reads_only,
    )
    driver = EMtraceDriver(
        sample=sample, timeout_s=batch_timeout,
        config=pipeline_config, anns=anns, threads=threads,
    )
    try:
        tables = realtime_run(cfg, driver, outdir)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    lat = tables["latency"].iloc[0]
    click.echo(f"Replayed {library}: {lat['batches']} batches, "
               f"median latency {lat['median_s']}s, ratio {lat['ratio']}")
    _echo_scores(tables, truth_taxon, read_min)
    click.echo(f"Tables → {outdir}/")


@cli.command(name="realtime-rescore")
@click.option("--outdir", required=True, type=click.Path(exists=True, file_okay=False),
              help="A finished `realtime` output directory (timeseries.tsv + abund/).")
@click.option("--library", required=True, help="Name used in every output table.")
@click.option("--truth-taxon", default=None)
@click.option("--truth", "truth_table", default=None,
              type=click.Path(exists=True, dir_okay=False))
@click.option("--panel", default=None, type=click.Path(exists=True, dir_okay=False))
@click.option("--rank", default="species", show_default=True)
@click.option("--score-col", default="readn", show_default=True)
@click.option("--prob-col", default="prob", show_default=True)
@click.option("--read-min", default=10.0, show_default=True)
@click.option("--tail", type=float, default=0.0, show_default=True)
@click.option("--tail-col", type=click.Choice(["readn", "readn_em"]),
              default="readn", show_default=True)
@click.option("--stabilise-on", default="L1", show_default=True)
@click.option("--stabilise-tol", default=0.02, show_default=True)
def realtime_rescore_cmd(outdir, library, truth_taxon, truth_table, panel, rank, score_col,
                         prob_col, read_min, tail, tail_col, stabilise_on, stabilise_tol):
    """Score a finished replay again without replaying it (Module C).

    Reads the replay's clocks from timeseries.tsv and its per-batch abundance
    tables from abund/, and rewrites the scoring tables in OUTDIR. The options
    mean what they mean for `realtime`.
    """
    cfg = RealtimeConfig(
        fastq_dir="", library=library, truth_taxon=truth_taxon, truth_table=truth_table,
        panel=panel, rank=rank, score_col=score_col, prob_col=prob_col, read_min=read_min,
        tail=tail, tail_col=tail_col, stabilise_on=stabilise_on, stabilise_tol=stabilise_tol,
    )
    try:
        tables = realtime_rescore(cfg, outdir)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    click.echo(f"Rescored {library}: {len(tables['timeseries'])} batches")
    _echo_scores(tables, truth_taxon, read_min)
    click.echo(f"Tables → {outdir}/")


def _echo_scores(tables: dict, truth_taxon: str | None, read_min: float) -> None:
    """Print a replay's headline scores: call, detection, accuracy or panel."""
    def _when(r) -> str:
        # --reads-only libraries have no sequencing clock; report compute time instead.
        if pd.isna(r["seq_time_s"]):
            return f"{r['wall_s']}s of compute"
        return f"{r['seq_time_s']}s of sequencing"

    fc = tables.get("first_call")
    if fc is not None and not fc.empty:
        r = fc.iloc[0]
        click.echo(f"First correct call: batch {r['batch']}, {r['reads']} reads, "
                   f"{_when(r)}")
    elif truth_taxon:
        click.echo(f"No persistent correct call for {truth_taxon!r}")
    fd = tables.get("first_detection")
    if fd is not None and not fd.empty:
        r = fd.iloc[0]
        click.echo(f"First detection (>= {read_min:g} reads): batch {r['batch']}, "
                   f"{r['reads']} reads, {_when(r)}")
    elif truth_taxon:
        click.echo(f"No persistent detection of {truth_taxon!r} at >= {read_min:g} reads")

    pt = tables.get("panel")
    if pt is not None:
        sm = tables["panel_summary"].iloc[0]
        click.echo(f"PCR targets detected {sm['pcr_pos_detected']}/{sm['pcr_pos_n']}; "
                   f"PCR negatives clear {sm['pcr_neg_clear']}/{sm['pcr_neg_n']}; "
                   f"study-reported taxa found {sm['reported_detected']}/{sm['reported_n']}")
        for _, r in pt.iterrows():
            first = ("" if r["detected"] == "no" else
                     f", from batch {r['first_detection_batch']} "
                     f"({r['first_detection_reads_cum']} reads)")
            click.echo(f"  {r['evidence']:<8} {r['taxon']}: {r['outcome']}, "
                       f"{r['reads_final']:g} reads{first}  [{r['source']}]")
        return

    acc = tables.get("accuracy")
    if acc is not None and not acc.empty:
        for _, r in acc.iterrows():
            click.echo(f"  {r['event']:<15} batch {r['batch']}: "
                       f"TP={r['TP']} FP={r['FP']} FN={r['FN']} F1={r['F1']}")


if __name__ == "__main__":
    cli()
