"""Module D (background-filtering efficiency) — unit + CLI tests.

A tiny synthetic simmet (10 reads) exercises the set-arithmetic categories and
the streaming fraction-removed maths without any external data:

    human:    h1 h2 h3              (3)   = Total − nohuman
    bacteria: b1 b2 b3              (3)   = nohuman − vi
    viruses:  v1 v2 (euk) r1 (retro) p1 (phiX)   (4) = vi
      vi_nophix_noretro: v1 v2      (2)
      phiX:              p1         (1)   species_id 10847
      HERV/retro:        r1         (1)   = vi − vinpnr − phix
    Viruses excl. Retroviridae: v1 v2 p1 (3)  = vinpnr ∪ phix  (phiX retained)
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from ngs_benchmark.bgfilter import metrics, truth
from ngs_benchmark.bgfilter.runner import BgfilterConfig, run, run_batch, run_division_batch
from ngs_benchmark.cli import cli

TOTAL = 10


def _write_truth(tmp_path: Path) -> dict[str, Path]:
    nohuman = tmp_path / "nohuman.tsv"
    nohuman.write_text(
        "readid\tsuperkingdom\n"
        "b1\tBacteria\nb2\tBacteria\nb3\tBacteria\n"
        "v1\tViruses\nv2\tViruses\nr1\tViruses\np1\tViruses\n"
    )
    vi = tmp_path / "vi.tsv"
    vi.write_text(
        "readid\tspecies_id\n"
        "v1\t111\nv2\t111\nr1\t222\np1\t10847\n"  # p1 is phiX
    )
    vinpnr = tmp_path / "vi_nophix_noretro.tsv"
    vinpnr.write_text("readid\nv1\nv2\n")
    master = tmp_path / "all.tsv"
    all_ids = ["h1", "h2", "h3", "b1", "b2", "b3", "v1", "v2", "r1", "p1"]
    master.write_text("readid\n" + "".join(f"{r}\n" for r in all_ids))
    return {"nohuman": nohuman, "vi": vi, "vinpnr": vinpnr, "master": master}


def _load(tmp_path: Path) -> truth.BgfilterTruth:
    f = _write_truth(tmp_path)
    return truth.load_truth(
        total=TOTAL, nohuman=f["nohuman"], vi=f["vi"], vi_nophix_noretro=f["vinpnr"],
    )


def _write_fastq(path: Path, ids: list[str], gz: bool = False) -> None:
    body = "".join(f"@{rid} some desc\nACGT\n+\nIIII\n" for rid in ids)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(body)
    else:
        path.write_text(body)


# --- truth loading + denominators -------------------------------------------

def test_truth_denominators(tmp_path: Path) -> None:
    t = _load(tmp_path)
    d = t.denominators
    assert d == {
        "total": 10, "homo_sapiens": 3, "bacteria": 3,
        "viruses": 4, "viruses_excl_retro": 3, "herv": 1,
    }


def test_phix_detected_by_species_id(tmp_path: Path) -> None:
    t = _load(tmp_path)
    assert t.phix == {"p1"}
    assert t.vi == {"v1", "v2", "r1", "p1"}


def test_count_data_lines(tmp_path: Path) -> None:
    f = _write_truth(tmp_path)
    assert truth.count_data_lines(f["master"]) == 10


# --- fastq id streaming ------------------------------------------------------

def test_iter_fastq_ids_first_token_plain_and_gz(tmp_path: Path) -> None:
    plain = tmp_path / "r.fastq"
    gzf = tmp_path / "r.fastq.gz"
    _write_fastq(plain, ["x1", "x2"], gz=False)
    _write_fastq(gzf, ["x1", "x2"], gz=True)
    assert list(metrics.iter_fastq_ids(plain)) == ["x1", "x2"]
    assert list(metrics.iter_fastq_ids(gzf)) == ["x1", "x2"]


def test_iter_fastq_ids_strip_mate(tmp_path: Path) -> None:
    fq = tmp_path / "r1.fastq"
    # paired-style ids with /1 suffix; base id (no suffix) should be yielded
    _write_fastq(fq, ["gi|1|x-35000/1", "gi|1|x-34998/1", "plain_id"])
    assert list(metrics.iter_fastq_ids(fq, strip_mate=False)) == [
        "gi|1|x-35000/1", "gi|1|x-34998/1", "plain_id"]
    assert list(metrics.iter_fastq_ids(fq, strip_mate=True)) == [
        "gi|1|x-35000", "gi|1|x-34998", "plain_id"]


def test_strip_mate_matches_per_fragment_truth(tmp_path: Path) -> None:
    # per-fragment truth (no /1); R1-only removed reads carry /1 → need strip_mate
    t = _load(tmp_path)
    fq = tmp_path / "removed_r1.fastq"
    _write_fastq(fq, ["h1/1", "h2/1", "b1/1"])     # R1 mate suffixes
    no_strip = metrics.count_removed(fq, t, strip_mate=False)
    assert no_strip["total"] == 3 and no_strip["homo_sapiens"] == 3  # all "non-nohuman"
    assert no_strip["bacteria"] == 0                                 # b1/1 ∉ truth ids
    strip = metrics.count_removed(fq, t, strip_mate=True)
    assert strip["homo_sapiens"] == 2 and strip["bacteria"] == 1     # h1,h2 human; b1 bacteria


# --- fraction-removed maths --------------------------------------------------

def test_count_removed_and_fractions(tmp_path: Path) -> None:
    t = _load(tmp_path)
    removed_fq = tmp_path / "removed.fastq"
    # filter removed all human + all bacteria + retro + phiX, kept euk viruses
    _write_fastq(removed_fq, ["h1", "h2", "h3", "b1", "b2", "b3", "r1", "p1"])

    removed = metrics.count_removed(removed_fq, t)
    assert removed == {
        "total": 8, "homo_sapiens": 3, "bacteria": 3,
        "viruses": 2, "viruses_excl_retro": 1, "herv": 1,
    }
    fr = metrics.fractions_removed(removed, t.denominators)
    assert fr["homo_sapiens"] == 100.0
    assert fr["bacteria"] == 100.0
    assert fr["viruses"] == 50.0
    assert fr["viruses_excl_retro"] == 33.333
    assert fr["herv"] == 100.0
    assert fr["total"] == 80.0


def test_classify_removed_row_shape(tmp_path: Path) -> None:
    t = _load(tmp_path)
    removed_fq = tmp_path / "removed.fastq"
    _write_fastq(removed_fq, ["h1", "b1"])
    row = metrics.classify_removed(removed_fq, t, "AS100")
    assert row["filter"] == "AS100"
    assert set(row) == {"filter", *metrics.CATEGORY_LABELS.keys()}


def test_fractions_none_when_denominator_zero(tmp_path: Path) -> None:
    # herv denom 0 → None rather than ZeroDivisionError
    t = truth.BgfilterTruth(total=4, nohuman={"v1"}, vi={"v1"},
                            vi_nophix_noretro={"v1"}, phix=set())
    fr = metrics.fractions_removed(
        {"total": 0, "homo_sapiens": 0, "bacteria": 0, "viruses": 0,
         "viruses_excl_retro": 0, "herv": 0},
        t.denominators,
    )
    assert fr["herv"] is None          # vi − vinpnr − phix = 0
    assert fr["bacteria"] is None      # nohuman − vi = 0


# --- runner ------------------------------------------------------------------

def test_run_batch_table(tmp_path: Path) -> None:
    t = _load(tmp_path)
    fq0 = tmp_path / "AS0.fastq"
    fq100 = tmp_path / "AS100.fastq"
    _write_fastq(fq0, ["h1", "h2", "h3", "b1", "b2", "b3", "v1", "v2", "r1", "p1"])  # all
    _write_fastq(fq100, ["h1", "h2", "h3"])                                          # human only
    man = pd.DataFrame([
        {"filter": "AS0", "removed_fastq": str(fq0)},
        {"filter": "AS100", "removed_fastq": str(fq100)},
    ])
    out = tmp_path / "out"
    df = run_batch(man, t, outdir=str(out))
    assert list(df["filter"]) == ["AS0", "AS100"]
    assert list(df.columns) == ["filter", *metrics.CATEGORY_LABELS.keys()]
    # AS0 removes everything → 100% in every category
    as0 = df.set_index("filter").loc["AS0"]
    assert as0["homo_sapiens"] == 100.0 and as0["viruses"] == 100.0
    # AS100 removes only human
    as100 = df.set_index("filter").loc["AS100"]
    assert as100["homo_sapiens"] == 100.0 and as100["viruses"] == 0.0
    assert (out / "table.tsv").exists()


def test_run_single(tmp_path: Path) -> None:
    t = _load(tmp_path)
    fq = tmp_path / "AS100.fastq"
    _write_fastq(fq, ["h1", "h2", "h3"])
    out = tmp_path / "single"
    row = run(BgfilterConfig(removed_fastq=str(fq), label="AS100"), t, str(out))
    assert row["homo_sapiens"] == 100.0
    assert (out / "table.tsv").exists()


# --- CLI ---------------------------------------------------------------------

def test_cli_bgfilter_batch_with_total(tmp_path: Path) -> None:
    f = _write_truth(tmp_path)
    fq0 = tmp_path / "AS0.fastq"
    _write_fastq(fq0, ["h1", "b1", "v1", "r1", "p1"])
    man = tmp_path / "manifest.tsv"
    pd.DataFrame([{"filter": "AS0", "removed_fastq": str(fq0)}]).to_csv(
        man, sep="\t", index=False)
    out = tmp_path / "cli_out"
    res = CliRunner().invoke(cli, [
        "bgfilter-batch", "--manifest", str(man), "--outdir", str(out),
        "--total", str(TOTAL),
        "--truth-nohuman", str(f["nohuman"]),
        "--truth-vi", str(f["vi"]),
        "--truth-vi-nophix-noretro", str(f["vinpnr"]),
    ])
    assert res.exit_code == 0, res.output
    df = pd.read_csv(out / "table.tsv", sep="\t")
    assert df.iloc[0]["filter"] == "AS0"
    assert df.iloc[0]["homo_sapiens"] == 33.333   # 1 of 3 humans removed (h1)
    assert df.iloc[0]["viruses"] == 75.0          # v1, r1, p1 of 4 viral removed


def test_cli_bgfilter_total_from_truth_all(tmp_path: Path) -> None:
    f = _write_truth(tmp_path)
    fq = tmp_path / "AS100.fastq"
    _write_fastq(fq, ["h1", "h2", "h3"])
    out = tmp_path / "cli_single"
    res = CliRunner().invoke(cli, [
        "bgfilter", "--removed-fastq", str(fq), "--label", "AS100", "--outdir", str(out),
        "--truth-all", str(f["master"]),
        "--truth-nohuman", str(f["nohuman"]),
        "--truth-vi", str(f["vi"]),
        "--truth-vi-nophix-noretro", str(f["vinpnr"]),
    ])
    assert res.exit_code == 0, res.output
    df = pd.read_csv(out / "table.tsv", sep="\t")
    assert df.iloc[0]["homo_sapiens"] == 100.0


def test_cli_bgfilter_batch_missing_column(tmp_path: Path) -> None:
    f = _write_truth(tmp_path)
    man = tmp_path / "bad.tsv"
    pd.DataFrame([{"filter": "AS0"}]).to_csv(man, sep="\t", index=False)
    res = CliRunner().invoke(cli, [
        "bgfilter-batch", "--manifest", str(man), "--outdir", str(tmp_path / "o"),
        "--total", str(TOTAL),
        "--truth-nohuman", str(f["nohuman"]),
        "--truth-vi", str(f["vi"]),
        "--truth-vi-nophix-noretro", str(f["vinpnr"]),
    ])
    assert res.exit_code != 0
    assert "removed_fastq" in res.output


# --- division-labeled truth (simulated-Nanopore datasets) --------------------

def _write_divtruth(tmp_path: Path) -> Path:
    p = tmp_path / "combined.truth.tsv"
    p.write_text(
        "readid\tdivision\tspecies\ttaxid\tsource\tmutation_rate\n"
        "v1\tViruses\tx\t1\tCOVID\t20\n"
        "v2\tViruses\tx\t1\tCOVID\t20\n"
        "b1\tBacteria\tx\t2\tSaureus\t20\n"
        "b2\tBacteria\tx\t2\tSaureus\t20\n"
        "b3\tBacteria\tx\t2\tSaureus\t20\n"
        "e1\tEukaryota\tx\t9606\thuman\t20\n"
        "e2\tEukaryota\tx\t9606\thuman\t20\n"
        "e3\tEukaryota\tx\t9606\thuman\t20\n"
        "e4\tEukaryota\tx\t9606\thuman\t20\n"
    )
    return p


def test_division_truth_denominators_and_order(tmp_path: Path) -> None:
    t = truth.load_division_truth(_write_divtruth(tmp_path))
    assert t.total == 9
    assert t.denominators == {"total": 9, "Viruses": 2, "Bacteria": 3, "Eukaryota": 4}
    # preferred biological order: Viruses, Bacteria, Eukaryota
    assert t.divisions == ["Viruses", "Bacteria", "Eukaryota"]


def test_count_removed_division(tmp_path: Path) -> None:
    t = truth.load_division_truth(_write_divtruth(tmp_path))
    removed_fq = tmp_path / "removed.fastq"
    # filter removed all eukaryota + all bacteria + 1 of 2 viruses
    _write_fastq(removed_fq, ["e1", "e2", "e3", "e4", "b1", "b2", "b3", "v1"])
    removed = metrics.count_removed_division(removed_fq, t)
    assert removed == {"total": 8, "Eukaryota": 4, "Bacteria": 3, "Viruses": 1}


def test_run_division_batch_table(tmp_path: Path) -> None:
    t = truth.load_division_truth(_write_divtruth(tmp_path))
    fq_lo = tmp_path / "AS0.fastq"      # removes everything
    fq_hi = tmp_path / "AS100.fastq"    # removes host only (no viruses)
    _write_fastq(fq_lo, ["v1", "v2", "b1", "b2", "b3", "e1", "e2", "e3", "e4"])
    _write_fastq(fq_hi, ["b1", "b2", "b3", "e1", "e2", "e3", "e4"])
    man = pd.DataFrame([
        {"filter": "AS0", "removed_fastq": str(fq_lo)},
        {"filter": "AS100", "removed_fastq": str(fq_hi)},
    ])
    out = tmp_path / "out"
    df = run_division_batch(man, t, outdir=str(out))
    assert list(df.columns) == ["filter", "total", "Viruses", "Bacteria", "Eukaryota"]
    as0 = df.set_index("filter").loc["AS0"]
    assert as0["Viruses"] == 100.0 and as0["Bacteria"] == 100.0 and as0["Eukaryota"] == 100.0
    as100 = df.set_index("filter").loc["AS100"]
    assert as100["Viruses"] == 0.0 and as100["Bacteria"] == 100.0 and as100["Eukaryota"] == 100.0
    assert (out / "table.tsv").exists()


def test_cli_bgfilter_divbatch(tmp_path: Path) -> None:
    tt = _write_divtruth(tmp_path)
    fq = tmp_path / "AS100.fastq"
    _write_fastq(fq, ["b1", "b2", "b3", "e1", "e2", "e3", "e4"])   # host only
    man = tmp_path / "manifest.tsv"
    pd.DataFrame([{"filter": "AS100", "removed_fastq": str(fq)}]).to_csv(
        man, sep="\t", index=False)
    out = tmp_path / "cli_out"
    res = CliRunner().invoke(cli, [
        "bgfilter-divbatch", "--manifest", str(man),
        "--truth-table", str(tt), "--outdir", str(out),
    ])
    assert res.exit_code == 0, res.output
    df = pd.read_csv(out / "table.tsv", sep="\t")
    assert df.iloc[0]["Viruses"] == 0.0
    assert df.iloc[0]["Eukaryota"] == 100.0


def test_cli_bgfilter_requires_total_or_truth_all(tmp_path: Path) -> None:
    f = _write_truth(tmp_path)
    fq = tmp_path / "AS100.fastq"
    _write_fastq(fq, ["h1"])
    res = CliRunner().invoke(cli, [
        "bgfilter", "--removed-fastq", str(fq), "--label", "AS100",
        "--outdir", str(tmp_path / "o"),
        "--truth-nohuman", str(f["nohuman"]),
        "--truth-vi", str(f["vi"]),
        "--truth-vi-nophix-noretro", str(f["vinpnr"]),
    ])
    assert res.exit_code != 0
    assert "total" in res.output.lower()
