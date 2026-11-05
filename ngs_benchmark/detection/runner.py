"""Compose §3.1–§3.6 into a single end-to-end run.

The runner owns the deterministic pipeline order:

    load → normalise → [Buddle filter] → set_metrics + per_taxon
                     → PR curve / AUPR
                     → dilution metrics (when --dilution-series)
                     → write outputs

§3.4 ROC is deferred from the beta slice (raises if --tn-universe is given).
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from pathlib import Path
import pandas as pd

from ngs_benchmark.normalise import NormaliseOpts, normalise
from ngs_benchmark.io.adapters import get_adapter
from ngs_benchmark.io.adapters import truth as truth_io
from ngs_benchmark.io.writers import write_tsv, write_json, write_xlsx_summary

from .buddle import BuddleOpts, annotate as buddle_annotate, filter_by_rule
from .curves import pr_curve
from .dilution import SamplePred, dilution_metrics
from .metrics import per_taxon_table, set_metrics
from ngs_benchmark.schema import UNKNOWN


@dataclass
class DetectionConfig:
    pred: str
    truth: str
    outdir: str
    pred_fmt: str = "emtrace"
    ranks: tuple[str, ...] = ("species", "genus", "family")
    taxgroups: tuple[str, ...] = ("Viruses",)
    score_col: str = "readn"
    read_min: float = 0.0          # applied to `filter_col`, never `tail_col`
    min_rel: float = 0.0           # relative-abundance floor, FRACTION (0.0001 = 0.01%)
    min_rel_exclude: tuple[str, ...] = ()  # divisions left out of the floor denominator
    filter_col: str = "readn"      # count column read_min and min_rel act on
    tail_vi: float = 0.0           # FRACTION of the library to trim from viral taxa
    tail_ba: float = 0.0           # FRACTION to trim from everything else
    tail_col: str = "readn"        # count column the tail trim ranks on
    exclude_retro: bool = False
    keep_phages: bool = False
    # "unknown" at a rank (reads placed above it, e.g. MEGAN-LR's LCA) is no
    # taxon call at that rank: the row is left out of that rank's sets
    drop_unknown: bool = False
    # Buddle mode
    buddle: bool = False
    nc_pred: str | None = None
    total_reads_sample: float | None = None
    total_reads_nc: float | None = None
    # Dilution series
    dilution_series: bool = False
    # ROC — deferred
    tn_universe: str | None = None
    extra: dict = field(default_factory=dict)


def _make_norm_opts(cfg: DetectionConfig) -> NormaliseOpts:
    """Pred-side normalisation: applies all filters."""
    return NormaliseOpts(
        keep_phages=cfg.keep_phages,
        exclude_retro=cfg.exclude_retro,
        read_min=cfg.read_min,
        min_rel=cfg.min_rel,
        min_rel_exclude=cfg.min_rel_exclude,
        filter_col=cfg.filter_col,
        tail_vi=cfg.tail_vi,
        tail_ba=cfg.tail_ba,
        tail_col=cfg.tail_col,
    )


def _make_truth_norm_opts(cfg: DetectionConfig) -> NormaliseOpts:
    """Truth-side normalisation: only structural filters (phage/retro) apply.

    `read_min` and tail-trim are read-count thresholds appropriate for noisy
    pipeline output, not for ground-truth presence/absence assertions. Membership
    truth files usually carry no readn at all (the loader fills it with 1.0).
    """
    return NormaliseOpts(
        keep_phages=cfg.keep_phages,
        exclude_retro=cfg.exclude_retro,
        read_min=0.0,
        tail_vi=0.0,
        tail_ba=0.0,
        tail_col=cfg.tail_col,
    )


def drop_unknown_labels(pred: pd.DataFrame, ranks: tuple[str, ...]) -> pd.DataFrame:
    """Set "unknown" labels at the scored ranks to NA, so the set metrics, the
    per-taxon table and the PR curve skip them. Read counts are unchanged, so
    filters already applied (read_min, min_rel) are unaffected."""
    out = pred.copy()
    for rank in ranks:
        if rank in out.columns:
            out[rank] = out[rank].mask(out[rank] == UNKNOWN)
    return out


def run(cfg: DetectionConfig) -> dict:
    """Run a detection benchmark end-to-end. Returns the summary dict."""
    if cfg.tn_universe is not None:
        raise NotImplementedError(
            "--tn-universe (§3.4 ROC mode) is deferred from the beta slice. "
            "Run without it to get PR curves and AUPR."
        )

    outdir = Path(cfg.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load.
    adapter = get_adapter(cfg.pred_fmt)
    pred_raw = adapter(cfg.pred)
    truth_raw = truth_io.load(cfg.truth)

    # Normalise. Pred uses full filter set; truth uses only the structural filters
    # (read_min / tail-trim do not apply to ground-truth presence assertions).
    pred = normalise(pred_raw, _make_norm_opts(cfg))
    truth = normalise(truth_raw, _make_truth_norm_opts(cfg))
    norm_opts = _make_norm_opts(cfg)  # passed through for the optional NC pred below

    summary: dict = {"config": {k: v for k, v in asdict(cfg).items() if k != "extra"}}

    # §3.5 Buddle filter (optional).
    if cfg.buddle:
        if cfg.nc_pred is None:
            raise ValueError("--buddle requires --nc-pred (negative-control prediction file)")
        if cfg.total_reads_sample is None or cfg.total_reads_nc is None:
            raise ValueError("--buddle requires --total-reads-sample and --total-reads-nc")
        nc_raw = adapter(cfg.nc_pred)
        nc = normalise(nc_raw, norm_opts)
        # Annotate + filter at species rank (RPMR join key).
        pred_pre_buddle = pred.copy()
        annotated = buddle_annotate(
            pred,
            neg_ctrl=nc,
            total_reads_sample=cfg.total_reads_sample,
            total_reads_nc=cfg.total_reads_nc,
            opts=BuddleOpts(),
            rank="species",
        )
        write_tsv(annotated, outdir / "buddle_annotated.tsv")
        pred = filter_by_rule(annotated)
        # Also emit pre-Buddle metrics for the reviewer-facing comparison.
        pre_metrics = set_metrics(pred_pre_buddle, truth, ranks=cfg.ranks, taxgroups=cfg.taxgroups)
        write_tsv(pre_metrics, outdir / "metrics_pre_buddle.tsv")

    pred_labelled = pred  # the Excel summary lists the prediction as loaded
    if cfg.drop_unknown:
        pred = drop_unknown_labels(pred, cfg.ranks)

    # §3.2 set metrics + per-taxon table.
    metrics_df = set_metrics(pred, truth, ranks=cfg.ranks, taxgroups=cfg.taxgroups)
    write_tsv(metrics_df, outdir / "metrics.tsv")

    per_taxon_df = per_taxon_table(pred, truth, ranks=cfg.ranks, taxgroups=cfg.taxgroups)
    write_tsv(per_taxon_df, outdir / "per_taxon.tsv")

    # Legacy-R Excel summary (§2.3 `summary.xlsx`).
    write_xlsx_summary(
        metrics_df=metrics_df,
        pred=pred_labelled,
        truth=truth,
        path=outdir / "summary.xlsx",
        taxgroups=cfg.taxgroups,
    )

    # §3.3 PR curve at the primary (taxgroup, rank): first taxgroup × "species".
    primary_rank = cfg.ranks[0] if cfg.ranks else "species"
    primary_group = cfg.taxgroups[0] if cfg.taxgroups else "Viruses"
    curve = pr_curve(
        pred,
        truth,
        rank=primary_rank,
        taxgroup=primary_group,
        score_col=cfg.score_col,
    )
    write_tsv(curve.prc, outdir / "prc.tsv")
    summary.update(
        {
            "primary_taxgroup": primary_group,
            "primary_rank": primary_rank,
            "AUPR": curve.aupr,
            "F_max": curve.f_max,
            "threshold_at_Fmax": curve.threshold_at_fmax,
            "Pr_at_Fmax": curve.pr_fmax,
            "Rc_at_Fmax": curve.rc_fmax,
            "metrics": metrics_df.to_dict(orient="records"),
        }
    )

    # §3.6 dilution series — handled by the batch runner; the single-sample CLI cannot infer
    # dilution structure from one sample.
    if cfg.dilution_series:
        raise ValueError(
            "--dilution-series requires the batch runner. Use ngs_benchmark.detection.runner.run_batch()."
        )

    write_json(summary, outdir / "summary.json")
    return summary


def run_batch(
    cfg_template: DetectionConfig,
    manifest: pd.DataFrame,
    truth_path: str,
    outdir: str,
) -> dict:
    """Run many predictions × one truth (§5 detection-batch).

    `manifest` columns: pipeline, sample_id, path, fmt, [gc_per_ml], [matrix_id],
    [nc_sample_id], [total_reads], [total_reads_nc].

    Per-sample outputs go to <outdir>/<pipeline>/<sample_id>/; a top-level
    `combined_metrics.tsv` and `sensitivity_vs_gcml.tsv` (when dilution columns
    are present) are written at <outdir>/.
    """
    outroot = Path(outdir)
    outroot.mkdir(parents=True, exist_ok=True)

    truth_raw = truth_io.load(truth_path)
    truth = normalise(truth_raw, _make_truth_norm_opts(cfg_template))
    truth_has_sample_id = "sample_id" in truth_raw.columns

    combined_rows: list[dict] = []
    samples_for_dilution: list[SamplePred] = []
    for _, row in manifest.iterrows():
        sub_outdir = outroot / row["pipeline"] / row["sample_id"]
        # Per-sample truth. When the truth table is keyed by `sample_id`, each
        # sample must be scored against ONLY its own rows. Passing the whole
        # truth file scores every sample against the union of all samples' taxa,
        # which understates recall/F1 for datasets whose samples carry different
        # viruses (e.g. Petersen: 4 single-virus libraries → n_true=4 not 1).
        sample_truth_path = truth_path
        if truth_has_sample_id:
            sub_outdir.mkdir(parents=True, exist_ok=True)
            sample_truth = truth_raw[
                truth_raw["sample_id"].astype("string") == str(row["sample_id"])
            ]
            sample_truth_path = str(sub_outdir / "truth.sample.tsv")
            sample_truth.to_csv(sample_truth_path, sep="\t", index=False)
        sub_cfg = DetectionConfig(
            pred=row["path"],
            truth=sample_truth_path,
            outdir=str(sub_outdir),
            pred_fmt=row.get("fmt", cfg_template.pred_fmt),
            ranks=cfg_template.ranks,
            taxgroups=cfg_template.taxgroups,
            score_col=cfg_template.score_col,
            read_min=cfg_template.read_min,
            min_rel=cfg_template.min_rel,
            min_rel_exclude=cfg_template.min_rel_exclude,
            filter_col=cfg_template.filter_col,
            tail_vi=cfg_template.tail_vi,
            tail_ba=cfg_template.tail_ba,
            tail_col=cfg_template.tail_col,
            exclude_retro=cfg_template.exclude_retro,
            keep_phages=cfg_template.keep_phages,
            drop_unknown=cfg_template.drop_unknown,
            buddle=False,           # Buddle-per-sample TODO; safe default for batch
            dilution_series=False,
        )
        s_summary = run(sub_cfg)
        for m in s_summary["metrics"]:
            row_out = {"pipeline": row["pipeline"], "sample_id": row["sample_id"], **m}
            combined_rows.append(row_out)

        if cfg_template.dilution_series and "gc_per_ml" in row.index:
            adapter = get_adapter(row.get("fmt", cfg_template.pred_fmt))
            pred = normalise(adapter(row["path"]), _make_norm_opts(cfg_template))
            if cfg_template.drop_unknown:
                pred = drop_unknown_labels(pred, cfg_template.ranks)
            samples_for_dilution.append(
                SamplePred(
                    sample_id=row["sample_id"],
                    gc_per_ml=float(row["gc_per_ml"]),
                    matrix_id=str(row["matrix_id"]) if "matrix_id" in row.index and pd.notna(row.get("matrix_id")) else None,
                    pred=pred,
                )
            )

    combined_df = pd.DataFrame(combined_rows)
    write_tsv(combined_df, outroot / "combined_metrics.tsv")

    summary: dict = {"n_samples": len(manifest)}
    if cfg_template.dilution_series and samples_for_dilution:
        prim_rank = cfg_template.ranks[0] if cfg_template.ranks else "species"
        prim_group = cfg_template.taxgroups[0] if cfg_template.taxgroups else "Viruses"
        dil_df, dil_summary = dilution_metrics(
            samples_for_dilution,
            truth,
            rank=prim_rank,
            taxgroup=prim_group,
        )
        write_tsv(dil_df, outroot / "sensitivity_vs_gcml.tsv")
        summary.update(dil_summary)

    write_json(summary, outroot / "summary.json")
    return summary
