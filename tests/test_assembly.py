"""Module B (assembly quality) — unit + integration + CLI tests.

Pure-Python units (FASTA parsing, identity maths, summary, viral_reads) run
everywhere. Tests that actually shell out to MAFFT are skipped when the binary
is absent.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from ngs_benchmark.assembly import metrics
from ngs_benchmark.assembly import mafft
from ngs_benchmark.assembly.runner import AssemblyConfig, count_viral_reads, run
from ngs_benchmark.cli import cli

HAVE_MAFFT = shutil.which("mafft") is not None
needs_mafft = pytest.mark.skipif(not HAVE_MAFFT, reason="mafft not on PATH")


# --- pure-Python units -------------------------------------------------------

def test_read_fasta_lengths_multisegment(tmp_path: Path) -> None:
    fa = tmp_path / "ref.fa"
    fa.write_text(">segL desc here\nACGTACGT\nACGT\n>segS\nGGGG\n")
    assert metrics.read_fasta_lengths(str(fa)) == {"segL": 12, "segS": 4}


@needs_mafft
def test_coverage_table_seqid_mismatch_single_segment(tmp_path: Path) -> None:
    # Reference seqid (synthetic mutant id) differs from consensus seqid
    # (original accession) — single-segment positional fallback must kick in.
    seq = "ACGTACGTACGTTGCATGCA" * 30
    ref = tmp_path / "ref.fa"
    con = tmp_path / "con.fa"
    ref.write_text(f">MUTANT_00 synthetic\n{seq}\n")
    con.write_text(f">OX637002.1 real accession\n{seq}\n")
    df = metrics.coverage_table(str(con), str(ref))
    assert df.iloc[0]["segment"] == "MUTANT_00"
    assert df.iloc[0]["coverage_pct"] == 100.0
    assert df.iloc[0]["positions_called_pct"] == 100.0


def test_read_fasta_seqs_multisegment(tmp_path: Path) -> None:
    fa = tmp_path / "ref.fa"
    fa.write_text(">segL desc\nacgt\nACGT\n>segS\ngggg\n")
    assert metrics.read_fasta_seqs(str(fa)) == {"segL": "ACGTACGT", "segS": "GGGG"}


def test_identity_counts_excludes_gaps_and_n() -> None:
    # 8 columns, one per case:
    #   col:  1    2    3    4    5     6     7     8
    #   ref:  A    C    G    T    A     N     -     A
    #   cons: A    C    G    A    N     C     C     -
    aln_ref = "ACGTAN-A"
    aln_cons = "ACGANCC-"
    covered, aligned, matching = mafft.identity_counts(aln_ref, aln_cons)
    # covered (both non-gap, cons non-N): cols 1-4 + col 6 (ref=N,cons=C) = 5
    #   col 5 (ref=A, cons=N) not covered — N treated as gap for coverage
    # aligned (both non-gap, non-N):  cols 1-4 = 4
    # matching:                       cols 1-3 = 3
    assert (covered, aligned, matching) == (5, 4, 3)


def test_identity_counts_perfect() -> None:
    covered, aligned, matching = mafft.identity_counts("ACGTACGT", "ACGTACGT")
    assert covered == 8 and aligned == 8 and matching == 8


def test_summarise_single_segment() -> None:
    cov = pd.DataFrame([{
        "segment": "s", "ref_len": 1000, "covered_bp": 700,
        "aligned_bp": 700, "matching_bp": 693,
        "coverage_pct": 70.0, "identity_pct": 99.0,
    }])
    s = metrics.summarise(cov)
    assert s["n_segments"] == 1
    assert s["weighted_coverage_pct"] == 70.0
    assert s["weighted_identity_pct"] == 99.0
    assert s["mean_coverage_pct"] == 70.0


def test_summarise_weighted_vs_mean_multisegment() -> None:
    # Long segment well covered, short segment poorly covered: weighted != mean.
    cov = pd.DataFrame([
        {"segment": "L", "ref_len": 900, "covered_bp": 900, "aligned_bp": 900,
         "matching_bp": 900, "coverage_pct": 100.0, "identity_pct": 100.0},
        {"segment": "S", "ref_len": 100, "covered_bp": 0, "aligned_bp": 0,
         "matching_bp": 0, "coverage_pct": 0.0, "identity_pct": 0.0},
    ])
    s = metrics.summarise(cov)
    assert s["mean_coverage_pct"] == 50.0            # (100 + 0) / 2
    assert s["weighted_coverage_pct"] == 90.0        # 900 / 1000


def test_summarise_na_when_no_metrics() -> None:
    df = pd.DataFrame([{
        "segment": "s", "ref_len": 1000, "covered_bp": 0, "aligned_bp": 0,
        "matching_bp": 0, "coverage_pct": None, "identity_pct": None,
    }])
    s = metrics.summarise(df)
    assert s["mean_coverage_pct"] is None
    assert s["weighted_coverage_pct"] is None
    assert s["total_ref_len"] == 1000


# --- viral_reads -------------------------------------------------------------

def test_count_viral_reads_states(tmp_path: Path) -> None:
    assert count_viral_reads(None) is None                      # column omitted
    assert count_viral_reads(str(tmp_path / "absent.csv")) == "no hits"

    f = tmp_path / "stats_div.csv"
    f.write_text("division,read_count\nViruses,1007\nBacteria,42\n")
    assert count_viral_reads(str(f)) == 1007

    no_vi = tmp_path / "stats_div_novi.csv"
    no_vi.write_text("division,read_count\nBacteria,42\n")
    assert count_viral_reads(str(no_vi)) == 0                   # Viruses row absent


# --- integration (mafft) --------------------------------------------------

def _write_seq(path: Path, name: str, seq: str) -> None:
    path.write_text(f">{name}\n{seq}\n")


@needs_mafft
def test_align_pair_one_substitution() -> None:
    # 60 bp; consensus has a single substitution → aligned 60, matching 59.
    seq = "ACGTACGTACGTTGCATGCA" * 3
    cons = seq[:30] + ("A" if seq[30] != "A" else "C") + seq[31:]
    aln_ref, aln_cons = mafft.align_pair(seq, cons)
    assert len(aln_ref) == len(aln_cons)
    covered, aligned, matching = mafft.identity_counts(aln_ref, aln_cons)
    assert covered == 60 and aligned == 60 and matching == 59


@needs_mafft
def test_align_pair_reverse_complemented_consensus() -> None:
    # A consensus on the opposite strand is flipped by MAFFT and scores 100%.
    seq = "ATGGCGTACGTTAGCCTAGGATCCGATCGTTAACGGCATGCAAGTCGACTTAAGCTTG" * 3
    rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    aln_ref, aln_cons = mafft.align_pair(seq, rc)
    covered, aligned, matching = mafft.identity_counts(aln_ref, aln_cons)
    assert covered == len(seq) and aligned == matching == len(seq)


@needs_mafft
def test_run_perfect_consensus(tmp_path: Path) -> None:
    # A ~600 bp sequence; consensus identical to reference → ~100% cov & identity.
    seq = ("ACGTACGTACGTTGCATGCA" * 30)
    ref = tmp_path / "ref.fa"
    con = tmp_path / "con.fa"
    _write_seq(ref, "refX", seq)
    _write_seq(con, "refX", seq)

    cfg = AssemblyConfig(consensus=str(con), reference=str(ref),
                         outdir=str(tmp_path / "out"), sample="perfect")
    summary = run(cfg)
    assert summary["sample"] == "perfect"
    # Global MAFFT alignment of identical sequences → exactly 100% both.
    assert summary["weighted_coverage_pct"] == 100.0
    assert summary["weighted_identity_pct"] == 100.0
    assert (tmp_path / "out" / "coverage.tsv").exists()
    assert (tmp_path / "out" / "summary.tsv").exists()


@needs_mafft
def test_run_missing_consensus_is_na(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    _write_seq(ref, "refX", "ACGTACGTACGT" * 50)
    cfg = AssemblyConfig(consensus=str(tmp_path / "nope.fa"), reference=str(ref),
                         outdir=str(tmp_path / "out"), sample="na")
    summary = run(cfg)
    assert summary["weighted_coverage_pct"] is None
    assert summary["weighted_identity_pct"] is None
    cov = pd.read_csv(tmp_path / "out" / "coverage.tsv", sep="\t")
    assert cov.iloc[0]["ref_len"] == 600


def test_run_missing_reference_raises(tmp_path: Path) -> None:
    cfg = AssemblyConfig(consensus=str(tmp_path / "c.fa"),
                         reference=str(tmp_path / "absent_ref.fa"),
                         outdir=str(tmp_path / "out"))
    with pytest.raises(FileNotFoundError):
        run(cfg)


@needs_mafft
def test_cli_assembly_e2e(tmp_path: Path) -> None:
    seq = "ACGTACGTACGTTGCATGCA" * 30
    ref = tmp_path / "ref.fa"
    con = tmp_path / "con.fa"
    stats_div = tmp_path / "stats_div.csv"
    _write_seq(ref, "refX", seq)
    _write_seq(con, "refX", seq)
    stats_div.write_text("division,read_count\nViruses,42\n")

    outdir = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "assembly",
        "--consensus", str(con),
        "--reference", str(ref),
        "--stats-div", str(stats_div),
        "--sample", "cliX",
        "--outdir", str(outdir),
    ])
    assert result.exit_code == 0, result.output
    summary = pd.read_csv(outdir / "summary.tsv", sep="\t")
    assert summary.iloc[0]["viral_reads"] == 42
    assert summary.iloc[0]["weighted_coverage_pct"] >= 98.0
    cov = pd.read_csv(outdir / "coverage.tsv", sep="\t")
    assert cov.iloc[0]["viral_reads"] == 42


@needs_mafft
def test_cli_assembly_batch_e2e(tmp_path: Path) -> None:
    seq = "ACGTACGTACGTTGCATGCA" * 30
    ref = tmp_path / "ref.fa"
    con = tmp_path / "con.fa"
    _write_seq(ref, "refX", seq)
    _write_seq(con, "refX", seq)

    manifest = tmp_path / "manifest.tsv"
    pd.DataFrame([
        {"sample_id": "s1", "consensus": str(con), "reference": str(ref)},
        {"sample_id": "s2", "consensus": str(tmp_path / "absent.fa"), "reference": str(ref)},
    ]).to_csv(manifest, sep="\t", index=False)

    outdir = tmp_path / "batch"
    result = CliRunner().invoke(cli, [
        "assembly-batch", "--manifest", str(manifest), "--outdir", str(outdir),
    ])
    assert result.exit_code == 0, result.output

    # Exactly two combined TSVs, no per-sample subdirectories.
    assert (outdir / "summary.tsv").exists()
    assert (outdir / "coverage.tsv").exists()
    assert not (outdir / "s1").exists()
    assert not (outdir / "s2").exists()

    combined = pd.read_csv(outdir / "summary.tsv", sep="\t")
    assert list(combined["sample"]) == ["s1", "s2"]
    assert combined.set_index("sample").loc["s1", "weighted_coverage_pct"] >= 98.0
    assert pd.isna(combined.set_index("sample").loc["s2", "weighted_coverage_pct"])

    # coverage.tsv carries a leading sample column, one row per segment per sample.
    cov = pd.read_csv(outdir / "coverage.tsv", sep="\t")
    assert cov.columns[0] == "sample"
    assert set(cov["sample"]) == {"s1", "s2"}


# --- positions_called (consensus N-counting) ---------------------------------

def test_read_fasta_ncounts_multirecord(tmp_path: Path) -> None:
    fa = tmp_path / "cons.fa"
    # segL: 12 bp, 3 N (incl lowercase); segS: 4 bp, 0 N
    fa.write_text(">segL desc\nACGTNN\nACGTnA\n>segS\nGGGG\n")
    assert metrics.read_fasta_ncounts(str(fa)) == {"segL": (12, 3), "segS": (4, 0)}


def test_positions_called_pct_helper() -> None:
    assert metrics._positions_called(1000, 800) == 80.0
    assert metrics._positions_called(0, 0) == 0.0


def test_summarise_positions_called_multisegment() -> None:
    import pandas as pd
    # 800/1000 called + 400/500 called → weighted = 1200/1500 = 80%; mean = 80%
    df = pd.DataFrame([
        {"segment": "L", "ref_len": 1000, "covered_bp": 900, "aligned_bp": 900,
         "matching_bp": 900, "coverage_pct": 90.0, "identity_pct": 100.0,
         "called_bp": 800, "positions_called_pct": 80.0},
        {"segment": "S", "ref_len": 500, "covered_bp": 400, "aligned_bp": 400,
         "matching_bp": 400, "coverage_pct": 80.0, "identity_pct": 100.0,
         "called_bp": 400, "positions_called_pct": 80.0},
    ], columns=metrics.COVERAGE_COLUMNS)
    s = metrics.summarise(df)
    assert s["mean_positions_called_pct"] == 80.0
    assert s["weighted_positions_called_pct"] == 80.0


def test_summarise_positions_called_na_without_column() -> None:
    import pandas as pd
    # called_bp all None (no consensus read) → positions_called summaries are None
    df = pd.DataFrame([{
        "segment": "s", "ref_len": 1000, "covered_bp": 0, "aligned_bp": 0,
        "matching_bp": 0, "coverage_pct": None, "identity_pct": None,
        "called_bp": None, "positions_called_pct": None,
    }], columns=metrics.COVERAGE_COLUMNS)
    s = metrics.summarise(df)
    assert s["mean_positions_called_pct"] is None
    assert s["weighted_positions_called_pct"] is None


@needs_mafft
def test_coverage_table_positions_called_from_consensus(tmp_path: Path) -> None:
    # 600 bp genome, consensus has a 60 bp N-run → 540/600 = 90% called.
    seq = "ACGTACGTACGTTGCATGCA" * 30          # 600 bp
    cons_seq = "N" * 60 + seq[60:]              # mask first 60 bp
    ref = tmp_path / "ref.fa"
    con = tmp_path / "con.fa"
    _write_seq(ref, "refX", seq)
    _write_seq(con, "refX", cons_seq)

    df = metrics.coverage_table(str(con), str(ref))
    assert df.iloc[0]["called_bp"] == 540
    assert df.iloc[0]["positions_called_pct"] == 90.0


@needs_mafft
def test_run_positions_called_end_to_end(tmp_path: Path) -> None:
    seq = "ACGTACGTACGTTGCATGCA" * 30
    cons_seq = "N" * 60 + seq[60:]
    ref = tmp_path / "ref.fa"
    con = tmp_path / "con.fa"
    _write_seq(ref, "refX", seq)
    _write_seq(con, "refX", cons_seq)

    cfg = AssemblyConfig(consensus=str(con), reference=str(ref),
                         outdir=str(tmp_path / "out"), sample="s")
    summary = run(cfg)
    assert summary["weighted_positions_called_pct"] == 90.0
    assert summary["mean_positions_called_pct"] == 90.0


@needs_mafft
def test_run_missing_consensus_positions_called_na(tmp_path: Path) -> None:
    ref = tmp_path / "ref.fa"
    _write_seq(ref, "refX", "ACGTACGTACGT" * 50)
    cfg = AssemblyConfig(consensus=str(tmp_path / "nope.fa"), reference=str(ref),
                         outdir=str(tmp_path / "out"), sample="na")
    summary = run(cfg)
    assert summary["weighted_positions_called_pct"] is None
    cov = pd.read_csv(tmp_path / "out" / "coverage.tsv", sep="\t")
    assert pd.isna(cov.iloc[0]["called_bp"])


@needs_mafft
def test_cli_assembly_positions_called_echo(tmp_path: Path) -> None:
    seq = "ACGTACGTACGTTGCATGCA" * 30
    cons_seq = "N" * 60 + seq[60:]
    ref = tmp_path / "ref.fa"
    con = tmp_path / "con.fa"
    _write_seq(ref, "refX", seq)
    _write_seq(con, "refX", cons_seq)

    outdir = tmp_path / "out"
    result = CliRunner().invoke(cli, [
        "assembly", "--consensus", str(con), "--reference", str(ref),
        "--sample", "x", "--outdir", str(outdir),
    ])
    assert result.exit_code == 0, result.output
    assert "positions_called=90.00%" in result.output
    summary = pd.read_csv(outdir / "summary.tsv", sep="\t")
    assert summary.iloc[0]["weighted_positions_called_pct"] == 90.0


def test_cli_assembly_batch_missing_column(tmp_path: Path) -> None:
    manifest = tmp_path / "bad.tsv"
    pd.DataFrame([{"sample_id": "s1", "consensus": "c.fa"}]).to_csv(
        manifest, sep="\t", index=False)
    result = CliRunner().invoke(cli, [
        "assembly-batch", "--manifest", str(manifest), "--outdir", str(tmp_path / "o"),
    ])
    assert result.exit_code != 0
    assert "reference" in result.output


# ---------------------------------------------------------------------------
# indel_counts — gap runs of the consensus against the reference
# ---------------------------------------------------------------------------

def test_indel_counts_runs_counted_once() -> None:
    # 2-bp deletion (one run) and 1-bp insertion (one run)
    ref  = "ACGTACGT-ACGT"
    cons = "AC--ACGTTACGT"
    assert mafft.indel_counts(ref, cons) == (1, 1)


def test_indel_counts_ignores_unaligned_ends_and_n() -> None:
    # leading/trailing consensus gaps (partial consensus) are outside the called
    # span; an N-masked stretch is skipped, never a deletion
    ref  = "ACGTACGTACGTACGT"
    cons = "--GTANNNACG-AC--"
    ins, dele = mafft.indel_counts(ref, cons)
    assert (ins, dele) == (0, 1)


def test_indel_counts_n_inside_insertion_does_not_split() -> None:
    # one 6-bp insertion whose middle is masked counts once; an all-N insertion
    # is not counted
    ref  = "ACGT------ACGT---ACGT"
    cons = "ACGTAANNAAACGTNNNACGT"
    assert mafft.indel_counts(ref, cons) == (1, 0)


def test_indel_counts_no_called_bases() -> None:
    assert mafft.indel_counts("ACGT", "NNNN") == (0, 0)


def test_summarise_indel_totals() -> None:
    import pandas as pd
    df = pd.DataFrame([
        {"segment": "L", "ref_len": 1000, "covered_bp": 900, "aligned_bp": 900,
         "matching_bp": 897, "coverage_pct": 90.0, "identity_pct": 99.67,
         "called_bp": 900, "positions_called_pct": 90.0, "insertions": 1, "deletions": 2},
        {"segment": "S", "ref_len": 500, "covered_bp": 400, "aligned_bp": 400,
         "matching_bp": 399, "coverage_pct": 80.0, "identity_pct": 99.75,
         "called_bp": 400, "positions_called_pct": 80.0, "insertions": 0, "deletions": 1},
    ], columns=metrics.COVERAGE_COLUMNS)
    s = metrics.summarise(df)
    assert (s["substitutions"], s["insertions"], s["deletions"]) == (4, 1, 3)
