"""Module C — replay ordering, clock handling and the three metrics.

The metrics are pure functions over snapshots, so everything here runs without a
pipeline. Two behaviours are pinned because getting them wrong would produce a
plausible, wrong number rather than an error:

* batch files must replay in numeric order (``_10`` after ``_9``);
* the first correct call must be the first one that is *never displaced*.
"""

from __future__ import annotations

import gzip
from datetime import datetime, timedelta

import pandas as pd
import pytest

from ngs_benchmark.realtime.metrics import (
    Snapshot,
    first_call,
    first_detection,
    latency,
    stabilised_at,
    timeseries,
)
from ngs_benchmark.realtime.replay import (
    batch_files,
    contiguity_gaps,
    replay,
    scan_library,
)
from ngs_benchmark.realtime.runner import RealtimeConfig, EMtraceDriver, score_snapshot

T0 = datetime.fromisoformat("2025-04-16T07:42:09+03:00")


def _write_batch(d, n, reads, t_offset_s, gz=False):
    """Write a MinKNOW-shaped batch file with *reads* reads at a given offset."""
    lines = []
    for i in range(reads):
        ts = (T0 + timedelta(seconds=t_offset_s + i)).isoformat()
        lines += [f"@r{n}_{i} runid=x start_time={ts} barcode=barcode06", "ACGT", "+", "IIII"]
    body = "\n".join(lines) + "\n"
    path = d / (f"FAX_pass_barcode06_{n}.fastq" + (".gz" if gz else ""))
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body)
    return path


# --- replay ordering and clocks ---------------------------------------------

def test_batches_replay_in_numeric_not_lexicographic_order(tmp_path):
    """_10 must come after _9; sorted() on names would put it after _1."""
    for n in (0, 1, 2, 9, 10, 11):
        _write_batch(tmp_path, n, reads=2, t_offset_s=n * 60)
    order = [int(p.name.rsplit("_", 1)[1].split(".")[0]) for p in batch_files(tmp_path)]
    assert order == [0, 1, 2, 9, 10, 11]


def test_scan_library_reports_reads_and_relative_seq_time(tmp_path):
    _write_batch(tmp_path, 0, reads=3, t_offset_s=0)
    _write_batch(tmp_path, 1, reads=5, t_offset_s=600, gz=True)
    batches = scan_library(tmp_path)

    assert [b.reads for b in batches] == [3, 5]
    assert batches[0].index == 1                      # 1-based replay position
    assert batches[0].seq_time_s == pytest.approx(2)  # last read of a 3-read batch
    assert batches[1].seq_time_s == pytest.approx(604)


def test_library_without_timestamps_is_refused(tmp_path):
    (tmp_path / "FAX_pass_barcode06_0.fastq").write_text("@r1\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="start_time"):
        scan_library(tmp_path)


def test_library_without_timestamps_replays_on_read_count(tmp_path):
    """SRA deposits: reads_only keeps read counts and leaves sequencing time empty."""
    (tmp_path / "SRR1_pass_0.fastq").write_text("@SRR1.1 1/1\nACGT\n+\nIIII\n" * 3)
    (tmp_path / "SRR1_pass_1.fastq").write_text("@SRR1.4 1/1\nACGT\n+\nIIII\n" * 2)
    batches = scan_library(tmp_path, reads_only=True)
    assert [b.reads for b in batches] == [3, 2]
    assert all(b.seq_time_s is None for b in batches)


def test_reads_only_refuses_a_timed_library(tmp_path):
    _write_batch(tmp_path, 0, reads=2, t_offset_s=0)
    with pytest.raises(ValueError, match="sequencing time"):
        scan_library(tmp_path, reads_only=True)


def test_combined_runs_are_detected_as_a_gap(tmp_path):
    """MSD79-2's failure mode: two flow-cell runs merged into one directory."""
    _write_batch(tmp_path, 0, reads=2, t_offset_s=0)
    _write_batch(tmp_path, 1, reads=2, t_offset_s=48 * 24 * 3600)   # 48 days later
    gaps = contiguity_gaps(scan_library(tmp_path))
    assert gaps and gaps[0][0] == 2
    assert gaps[0][1] == pytest.approx(48 * 24, rel=0.01)


def test_replay_releases_one_batch_at_a_time(tmp_path):
    src, watch = tmp_path / "src", tmp_path / "watch"
    src.mkdir()
    for n in (0, 1, 2):
        _write_batch(src, n, reads=4, t_offset_s=n * 10)

    seen: list[int] = []

    def await_batch(index: int) -> None:
        # the file for this batch must be there, and no later one
        seen.append(len(list(watch.glob("*.fastq"))))

    releases = list(replay(scan_library(src), watch, await_batch))
    assert seen == [1, 2, 3]                       # released strictly one by one
    assert [r.reads_cum for r in releases] == [4, 8, 12]
    assert all(r.latency_s >= 0 for r in releases)


# --- metrics -----------------------------------------------------------------

def _snap(batch, top, reads, seq_t=None, prob=None, L1=None, lat=1.0):
    return Snapshot(
        batch=batch, seq_time_s=seq_t, wall_s=batch * 1.0, reads_cum=reads,
        latency_s=lat, top_taxon=top, em_prob=prob, L1=L1,
    )


def test_first_call_ignores_a_transient_correct_hit():
    """A lucky batch-1 hit that is later displaced is not the call."""
    snaps = [
        _snap(1, "Target virus", 100),     # correct, but displaced below
        _snap(2, "Homo sapiens", 200),
        _snap(3, "Target virus", 300),     # from here on it sticks
        _snap(4, "Target virus", 400),
    ]
    row = first_call(snaps, "Target virus", "L1").iloc[0]
    assert row["batch"] == 3
    assert row["reads"] == 300


def test_first_call_empty_when_never_settles():
    snaps = [_snap(1, "Target virus", 100), _snap(2, "Other", 200)]
    assert first_call(snaps, "Target virus", "L1").empty


def test_first_call_counts_unmapped_prefix_as_wrong():
    snaps = [_snap(1, None, 50), _snap(2, "Target virus", 100)]
    assert first_call(snaps, "Target virus", "L1").iloc[0]["batch"] == 2


def _det(batch, tp, top="Psychrobacter", reads=None):
    s = _snap(batch, top, reads or batch * 100)
    s.TP = tp
    return s


def test_first_detection_ignores_rank():
    """Field library: the virus passes read_min but a bacterium stays on top."""
    snaps = [_det(1, 0), _det(2, 1), _det(3, 1)]
    assert first_call(snaps, "Orthonairovirus", "ticks").empty
    row = first_detection(snaps, "ticks").iloc[0]
    assert row["batch"] == 2
    assert row["reads"] == 200


def test_first_detection_requires_persistence():
    snaps = [_det(1, 1), _det(2, 0), _det(3, 1), _det(4, 1)]
    assert first_detection(snaps, "ticks").iloc[0]["batch"] == 3


def test_first_detection_empty_when_lost_or_never_scored():
    assert first_detection([_det(1, 1), _det(2, 0)], "ticks").empty
    assert first_detection([_snap(1, "X", 100)], "ticks").empty   # TP None


def test_latency_ratio_says_whether_the_pipeline_keeps_up():
    # 2 batches, 10 s of compute, covering 100 s of sequencing -> ratio 10
    snaps = [_snap(1, "x", 1000, seq_t=50, lat=5.0), _snap(2, "x", 2000, seq_t=100, lat=5.0)]
    row = latency(snaps, "L1").iloc[0]
    assert row["compute_time_s"] == pytest.approx(10.0)
    assert row["ratio"] == pytest.approx(10.0)
    assert row["s_per_1k_reads"] == pytest.approx(5.0)


def test_latency_ratio_is_empty_without_sequencing_times():
    row = latency([_snap(1, "x", 1000, lat=2.0)], "L").iloc[0]
    assert row["ratio"] is None


def test_timeseries_has_one_row_per_batch_in_order():
    df = timeseries([_snap(1, "a", 10), _snap(2, "b", 20)], "L")
    assert list(df["batch"]) == [1, 2]
    assert list(df["library"]) == ["L", "L"]


def test_stabilised_at_finds_the_settling_point():
    snaps = [
        _snap(1, "a", 100, L1=0.9), _snap(2, "a", 200, L1=0.4),
        _snap(3, "a", 300, L1=0.26), _snap(4, "a", 400, L1=0.25),
    ]
    assert stabilised_at(snaps, tol=0.02).batch == 3
    assert stabilised_at(snaps, tol=0.001).batch == 4


# --- snapshot scoring --------------------------------------------------------

def test_score_snapshot_picks_the_dominant_taxon():
    abund = pd.DataFrame({
        "taxid": ["1", "2"], "species": ["Target virus", "Homo sapiens"],
        "readn": [900, 100], "prob": [0.99, 0.01],
    })
    cfg = RealtimeConfig(fastq_dir=".", library="L")
    got = score_snapshot(abund, cfg, truth=None)
    assert got["top_taxon"] == "Target virus"
    assert got["em_prob"] == pytest.approx(0.99)


def test_score_snapshot_tolerates_an_empty_table():
    cfg = RealtimeConfig(fastq_dir=".", library="L")
    assert score_snapshot(None, cfg, truth=None)["top_taxon"] is None
    assert score_snapshot(pd.DataFrame(), cfg, truth=None)["top_taxon"] is None


# --- accuracy at the milestones ----------------------------------------------

def test_set_scores_matches_framework_a_set_metrics():
    """Module C must not invent its own definition of precision and recall."""
    from ngs_benchmark.detection.metrics import set_metrics
    from ngs_benchmark.realtime.metrics import set_scores

    pred = pd.DataFrame({
        "species": ["a", "b", "x"], "division": ["Bacteria"] * 3,
    })
    truth = pd.DataFrame({
        "species": ["a", "b", "c"], "division": ["Bacteria"] * 3,
    })
    ref = set_metrics(pred, truth, ranks=("species",), taxgroups=("Bacteria",)).iloc[0]
    got = set_scores({"a", "b", "x"}, {"a", "b", "c"})

    assert (got["TP"], got["FP"], got["FN"]) == (ref["TP"], ref["FP"], ref["FN"])
    assert got["Pr"] == pytest.approx(ref["Pr"], abs=1e-4)
    assert got["Rc"] == pytest.approx(ref["Rc"], abs=1e-4)
    assert got["F1"] == pytest.approx(ref["F1"], abs=1e-4)


def test_accuracy_reports_the_three_milestones():
    snaps = [
        _snap(1, "X", 100, L1=0.9), _snap(2, "T", 200, L1=0.5),
        _snap(3, "T", 300, L1=0.26), _snap(4, "T", 400, L1=0.25),
    ]
    from ngs_benchmark.realtime.metrics import accuracy
    df = accuracy(snaps, "L1", truth_taxon="T")
    assert list(df["event"]) == ["first_call", "stabilised", "final"]
    assert list(df["batch"]) == [2, 3, 4]


def test_accuracy_omits_events_that_did_not_happen():
    from ngs_benchmark.realtime.metrics import accuracy
    snaps = [_snap(1, "X", 100), _snap(2, "X", 200)]     # never correct, no L1
    df = accuracy(snaps, "L", truth_taxon="T")
    assert list(df["event"]) == ["final"]


def test_accuracy_reports_first_detection_for_a_single_taxon():
    from ngs_benchmark.realtime.metrics import accuracy
    snaps = [_det(1, 0), _det(2, 1), _det(3, 1)]
    df = accuracy(snaps, "ticks", truth_taxon="Orthonairovirus", stabilise_on="F1")
    assert "first_call" not in set(df["event"])
    assert df.loc[df["event"] == "first_detection", "batch"].tolist() == [2]


def test_read_min_keeps_noise_taxa_out_of_precision():
    """Without a threshold every one-read taxon is a false positive."""
    abund = pd.DataFrame({
        "taxid": ["1", "2", "3"],
        "species": ["Target virus", "Noise one", "Noise two"],
        "readn": [500, 3, 1],
    })
    cfg = RealtimeConfig(fastq_dir=".", library="L", truth_taxon="Target virus", read_min=10)
    got = score_snapshot(abund, cfg, truth=None)
    assert (got["TP"], got["FP"], got["FN"]) == (1, 0, 0)
    assert got["Pr"] == 1.0

    loose = RealtimeConfig(fastq_dir=".", library="L", truth_taxon="Target virus", read_min=0)
    assert score_snapshot(abund, loose, truth=None)["FP"] == 2


def test_single_taxon_truth_gives_coarse_but_correct_scores():
    """The single-virus regime: one true taxon, so recall is 0 or 1."""
    abund = pd.DataFrame({
        "taxid": ["1", "2"], "species": ["Other virus", "Target virus"],
        "readn": [800, 400],
    })
    cfg = RealtimeConfig(fastq_dir=".", library="L", truth_taxon="Target virus")
    got = score_snapshot(abund, cfg, truth=None)
    assert got["n_true"] == 1
    assert got["Rc"] == 1.0          # the truth taxon is reported ...
    assert got["Pr"] == 0.5          # ... beside one false positive
    assert got["top_taxon"] == "Other virus"   # but it is not yet the top call


# --- EMtraceDriver: the contract with `emtrace watch` --------------------

def test_driver_command_passes_benchmark_settings(tmp_path):
    drv = EMtraceDriver(sample="lib", config="cfg.yaml", anns="16S.rrnDB.dedup", threads=8)
    cmd = drv.command(tmp_path / "watch", tmp_path / "pipeline")
    assert cmd[:2] == ["emtrace", "watch"]
    assert "--out" not in cmd
    # start-up handshake: log the "Watching ..." line, and catch early files
    assert "-v" in cmd and "--scan-existing" in cmd
    for flag, value in (("--dir", str(tmp_path / "watch")), ("-r", str(tmp_path / "pipeline")),
                        ("-s", "lib"), ("--config", "cfg.yaml"),
                        ("--anns", "16S.rrnDB.dedup"), ("-t", "8")):
        assert cmd[cmd.index(flag) + 1] == value


def test_driver_command_uses_only_real_watch_options(tmp_path):
    main = pytest.importorskip("emtrace.__main__")
    watch_opts = {o for p in main.cli.commands["watch"].params for o in p.opts}
    cmd = EMtraceDriver(config="c.yaml", anns="a", threads=2).command(tmp_path, tmp_path)
    flags = [tok for tok in cmd[2:] if tok.startswith("-")]
    assert flags and set(flags) <= watch_opts, set(flags) - watch_opts


def test_driver_reads_outputs_under_the_sample_dir(tmp_path):
    drv = EMtraceDriver(sample="lib")
    drv._outdir = tmp_path
    res = tmp_path / "lib" / "results"
    res.mkdir(parents=True)
    pd.DataFrame({"taxid": ["1"], "readn": [5]}).to_csv(res / "abund_table.tsv", sep="\t", index=False)
    assert drv.abundance()["readn"].tolist() == [5]
    assert drv._batch_dir(3) == tmp_path / "lib" / "batches" / "batch_0003"


def test_driver_waits_for_the_batch_callback_line(tmp_path):
    drv = EMtraceDriver(poll_s=0.01, timeout_s=0.2)
    drv._outdir = tmp_path
    drv.log_path.write_text("\n--- Batch 1 ---\nFiles: 1\n--- Batch 10 ---\n")
    drv.await_batch(1)
    drv.await_batch(10)
    with pytest.raises(TimeoutError):
        drv.await_batch(2)  # "--- Batch 2" must not match a prefix of another batch


def test_driver_raises_at_once_when_a_batch_fails(tmp_path):
    drv = EMtraceDriver(poll_s=0.01, timeout_s=5)
    drv._outdir = tmp_path
    drv.log_path.write_text("ERROR | Batch 4 pipeline failed: boom\n")
    with pytest.raises(RuntimeError, match="batch 4 pipeline failed"):
        drv.await_batch(4)


def test_zero_abundance_truth_rows_are_not_expected_taxa():
    truth = pd.DataFrame({
        "genus": ["Bacillus", "Listeria", "Saccharomyces"],
        "division": "Bacteria",
        "pct_16S_only": [60.0, 40.0, 0.0],
    })
    abund = pd.DataFrame({"genus": ["Bacillus", "Listeria"], "readn": [60, 40], "prob": [0.6, 0.4]})
    cfg = RealtimeConfig(fastq_dir=".", library="x", rank="genus")
    out = score_snapshot(abund, cfg, truth)
    assert (out["TP"], out["FP"], out["FN"]) == (2, 0, 0)


# --- field-library panel -----------------------------------------------------

from ngs_benchmark.realtime.metrics import panel, panel_summary  # noqa: E402
from ngs_benchmark.realtime.runner import load_panel, rescore  # noqa: E402

PANEL_TSV = (
    "# 58-37, Ergunay 2023\n"
    "taxon\tevidence\tsource\n"
    "Orthonairovirus\tpcr_pos\tTable S3 (PCR)\n"
    "Jingmenvirus\tpcr_pos\tTable S2 (PCR)\n"
    "Mivirus\treported\tTable S4 (BLAST)\n"
    "Phlebovirus\tpcr_neg\tnot listed\n"
)


def _panel(tmp_path, text=PANEL_TSV):
    path = tmp_path / "panel.tsv"
    path.write_text(text)
    return path


def test_load_panel_defaults_rank_and_skips_comments(tmp_path):
    df = load_panel(_panel(tmp_path), rank="genus")
    assert list(df["taxon"]) == ["Orthonairovirus", "Jingmenvirus", "Mivirus", "Phlebovirus"]
    assert set(df["rank"]) == {"genus"}
    assert df.loc[2, "source"] == "Table S4 (BLAST)"


def test_load_panel_refuses_unknown_evidence_and_duplicates(tmp_path):
    with pytest.raises(ValueError, match="unknown evidence"):
        load_panel(_panel(tmp_path, "taxon\tevidence\nX\tblast\n"), rank="genus")
    with pytest.raises(ValueError, match="twice"):
        load_panel(_panel(tmp_path, "taxon\tevidence\nX\tpcr_pos\nX\treported\n"), rank="genus")


def test_score_snapshot_counts_panel_taxa_without_set_metrics(tmp_path):
    abund = pd.DataFrame({
        "taxid": ["1", "2", "3", "4"],
        "species": ["Aigai virus", "CCHFV", "Jingmen tick virus", "Psychrobacter sp."],
        "genus": ["Orthonairovirus", "Orthonairovirus", "Jingmenvirus", "Psychrobacter"],
        "readn": [30, 5, 3, 500], "prob": [0.05, 0.01, 0.01, 0.9],
    })
    cfg = RealtimeConfig(fastq_dir=".", library="L", rank="genus")
    got = score_snapshot(abund, cfg, truth=None, targets=load_panel(_panel(tmp_path), "genus"))
    assert got["panel_reads"] == {"Orthonairovirus": 35.0, "Jingmenvirus": 3.0,
                                  "Mivirus": 0.0, "Phlebovirus": 0.0}
    assert got["top_taxon"] == "Psychrobacter"
    assert got["Pr"] is None and got["Rc"] is None and got["FP"] is None


def _psnap(batch, **reads):
    s = _snap(batch, "Psychrobacter", batch * 4000)
    s.panel_reads = reads
    return s


def test_panel_scores_each_evidence_class(tmp_path):
    targets = load_panel(_panel(tmp_path), rank="genus")
    snaps = [
        _psnap(1, Orthonairovirus=15, Jingmenvirus=4, Mivirus=12, Phlebovirus=0),
        _psnap(2, Orthonairovirus=40, Jingmenvirus=11, Mivirus=9, Phlebovirus=2),
        _psnap(3, Orthonairovirus=60, Jingmenvirus=20, Mivirus=14, Phlebovirus=3),
    ]
    df = panel(snaps, targets, "L", read_min=10).set_index("taxon")
    assert df.loc["Orthonairovirus", ["detected", "outcome", "first_detection_batch"]].tolist() \
        == ["yes", "detected", 1]
    assert df.loc["Jingmenvirus", "first_detection_batch"] == 2
    # Mivirus dips below read_min at batch 2: detection counts from batch 3
    assert df.loc["Mivirus", ["outcome", "first_detection_batch"]].tolist() == ["found", 3]
    assert df.loc["Phlebovirus", ["detected", "outcome"]].tolist() == ["no", "clear"]
    assert pd.isna(df.loc["Phlebovirus", "first_detection_batch"])
    assert df.loc["Phlebovirus", "reads_final"] == 3

    sm = panel_summary(df.reset_index(), "L").iloc[0]
    assert (sm["pcr_pos_detected"], sm["pcr_pos_n"]) == (2, 2)
    assert (sm["pcr_neg_clear"], sm["pcr_neg_n"]) == (1, 1)
    assert (sm["reported_detected"], sm["reported_n"]) == (1, 1)


def test_panel_flags_a_pcr_negative_that_is_detected(tmp_path):
    targets = load_panel(_panel(tmp_path, "taxon\tevidence\nOrthonairovirus\tpcr_neg\n"), "genus")
    row = panel([_psnap(1, Orthonairovirus=12)], targets, "L", read_min=10).iloc[0]
    assert (row["detected"], row["outcome"]) == ("yes", "false_detection")
    assert panel_summary(pd.DataFrame([row]), "L").iloc[0]["pcr_neg_clear"] == 0


def test_panel_excludes_a_truth(tmp_path):
    cfg = RealtimeConfig(fastq_dir=".", library="L", truth_taxon="X", panel=str(_panel(tmp_path)))
    with pytest.raises(ValueError, match="--panel"):
        rescore(cfg, tmp_path)


def test_rescore_rebuilds_tables_from_abund_snapshots(tmp_path):
    out = tmp_path / "run"
    (out / "abund").mkdir(parents=True)
    pd.DataFrame({
        "library": "L", "batch": [1, 2], "seq_time_s": [None, None], "wall_s": [60.0, 120.0],
        "reads_cum": [4000, 8000], "latency_s": [60.0, 60.0],
    }).to_csv(out / "timeseries.tsv", sep="\t", index=False)
    for b, n in ((1, 5), (2, 25)):
        pd.DataFrame({
            "taxid": ["1", "2"], "genus": ["Orthonairovirus", "Psychrobacter"],
            "readn": [n, 300], "prob": [0.1, 0.9],
        }).to_csv(out / "abund" / f"batch_{b:04d}.tsv", sep="\t", index=False)
    (out / "first_call.tsv").write_text("stale\n")     # from an earlier truth-taxon scoring

    cfg = RealtimeConfig(fastq_dir="", library="L", rank="genus",
                         panel=str(_panel(tmp_path, "taxon\tevidence\nOrthonairovirus\tpcr_pos\n")))
    tables = rescore(cfg, out)
    row = tables["panel"].iloc[0]
    assert (row["detected"], row["first_detection_batch"], row["reads_final"]) == ("yes", 2, 25)
    assert tables["timeseries"]["Pr"].isna().all()
    assert (out / "panel.tsv").is_file() and not (out / "first_call.tsv").exists()
