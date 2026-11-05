"""run_batch must score each sample against its own sample_id truth rows.

Regression test for the Petersen scoring bug: a truth table keyed by
`sample_id` was applied whole to every sample, so single-virus libraries were
scored against the union of all samples' viruses (n_true inflated, recall/F1
understated).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ngs_benchmark.detection.runner import DetectionConfig, run_batch


def _write_pred(path: Path, species: str, taxid: int) -> None:
    pd.DataFrame(
        {
            "readn": [100],
            "taxid": [taxid],
            "species": [species],
            "genus": [f"G_{species}"],
            "family": [f"F_{species}"],
            "order": ["OV"],
            "class": ["CV"],
            "division": ["Viruses"],
            "prob": [0.99],
        }
    ).to_csv(path, sep="\t", index=False)


def test_run_batch_filters_truth_by_sample_id(tmp_path: Path) -> None:
    # Two single-virus samples: s1 -> vA, s2 -> vB.
    p1 = tmp_path / "s1.tsv"
    p2 = tmp_path / "s2.tsv"
    _write_pred(p1, "vA", 1)
    _write_pred(p2, "vB", 2)

    # Truth keyed by sample_id; union across the file = {vA, vB}.
    truth_path = tmp_path / "truth.tsv"
    pd.DataFrame(
        {
            "sample_id": ["s1", "s2"],
            "taxid": [1, 2],
            "species": ["vA", "vB"],
            "genus": ["G_vA", "G_vB"],
            "family": ["F_vA", "F_vB"],
            "order": ["OV", "OV"],
            "class": ["CV", "CV"],
            "division": ["Viruses", "Viruses"],
        }
    ).to_csv(truth_path, sep="\t", index=False)

    manifest = pd.DataFrame(
        {
            "pipeline": ["emtrace", "emtrace"],
            "sample_id": ["s1", "s2"],
            "path": [str(p1), str(p2)],
            "fmt": ["emtrace", "emtrace"],
        }
    )

    out = tmp_path / "out"
    cfg = DetectionConfig(
        pred="",
        truth=str(truth_path),
        outdir=str(out),
        ranks=("species",),
        taxgroups=("Viruses",),
    )
    run_batch(cfg, manifest, truth_path=str(truth_path), outdir=str(out))

    combined = pd.read_csv(out / "combined_metrics.tsv", sep="\t")
    sp = combined[combined["rank"] == "species"].set_index("sample_id")

    # Each sample is scored against ITS OWN single virus -> n_true == 1, perfect.
    for sid in ("s1", "s2"):
        assert sp.loc[sid, "n_true"] == 1, combined.to_string()
        assert sp.loc[sid, "TP"] == 1
        assert sp.loc[sid, "FN"] == 0
        assert sp.loc[sid, "Pr"] == 1.0
        assert sp.loc[sid, "Rc"] == 1.0
        assert sp.loc[sid, "F1"] == 1.0

    # The per-sample truth subset is written for provenance.
    assert (out / "emtrace" / "s1" / "truth.sample.tsv").exists()
