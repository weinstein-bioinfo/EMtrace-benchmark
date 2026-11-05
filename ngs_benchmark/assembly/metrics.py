"""Genome coverage and identity from a consensus aligned to a reference.

Per reference record ("segment") the consensus is globally aligned to the
reference with MAFFT (ginsi), and we report:

    coverage_pct = covered_bp / ref_len      × 100   (breadth: ref positions spanned)
    identity_pct = matching_bp / aligned_bp  × 100   (accuracy of called bases)

`covered_bp` is the number of reference positions aligned to a non-gap, non-N
consensus base. `insertions` / `deletions` count indel events (gap runs) of the
consensus against the reference in the same alignment (`mafft.indel_counts`);
substitutions are `aligned_bp - matching_bp`. `aligned_bp` is the number of columns where both sequences are
non-gap and non-N — the identity denominator — and `matching_bp` is how many of
those agree. Consensus N positions are excluded from identity (completeness is
reported separately as positions_called); the alignment + identity logic is
ported from virefset (see `mafft.py`).

MAFFT replaced minimap2 here: minimap2 PAF column 10 without `-c` is a seed-based
estimate, not true identity, and badly under-reported it for divergent consensus.
"""

from __future__ import annotations

import pandas as pd

from .mafft import align_pair, identity_counts, indel_counts

COVERAGE_COLUMNS = [
    "segment", "ref_len", "covered_bp", "aligned_bp", "matching_bp",
    "coverage_pct", "identity_pct", "called_bp", "positions_called_pct",
    "insertions", "deletions",
]


def read_fasta_lengths(path: str) -> dict[str, int]:
    """Return ``{seqid: length}`` for each FASTA record.

    ``seqid`` is the first whitespace-delimited token of the header (so reference
    and consensus records are matched by accession).
    """
    lengths: dict[str, int] = {}
    name: str | None = None
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name = line[1:].split()[0] if len(line) > 1 else ""
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def read_fasta_ncounts(path: str) -> dict[str, tuple[int, int]]:
    """Return ``{seqid: (length, n_count)}`` for each FASTA record.

    ``seqid`` is the first whitespace-delimited token of the header (matching
    ``read_fasta_lengths``). ``n_count`` counts ``N``/``n`` (not-called bases);
    every other character counts toward ``length``. Used to derive
    positions_called directly from a consensus FASTA — the authoritative final
    output — rather than an intermediate coverage mask.
    """
    counts: dict[str, tuple[int, int]] = {}
    name: str | None = None
    length = 0
    ncount = 0
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                counts[name] = (length, ncount)
            name = line[1:].split()[0] if len(line) > 1 else ""
            length = 0
            ncount = 0
        else:
            seq = line.strip()
            length += len(seq)
            ncount += seq.count("N") + seq.count("n")
    if name is not None:
        counts[name] = (length, ncount)
    return counts


def read_fasta_seqs(path: str) -> dict[str, str]:
    """Return ``{seqid: sequence}`` (uppercased) for each FASTA record."""
    seqs: dict[str, list[str]] = {}
    name: str | None = None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0] if len(line) > 1 else ""
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def _positions_called(rlen: int, called_bp: int) -> float:
    """positions_called_pct for one segment: called bases / reference length."""
    return round(called_bp / rlen * 100, 2) if rlen else 0.0


def coverage_table(consensus: str, reference: str,
                   num_threads: int = 1) -> pd.DataFrame:
    """Per-segment coverage/identity/positions-called table.

    Every reference record gets a row; a segment absent from the consensus reports
    zeros so the table always reflects the full reference (multi-segment genomes
    included). Each present segment is globally aligned to its reference with
    MAFFT and scored by `identity_counts` (consensus N excluded from identity).

    ``called_bp`` and ``positions_called_pct`` come straight from the consensus
    FASTA: ``called_bp`` is the count of non-``N`` bases in the consensus record
    (ref-based assembly produces a reference-length consensus, so
    ``positions_called_pct = called_bp / ref_len × 100``).
    """
    ref_lengths = read_fasta_lengths(reference)
    ref_seqs = read_fasta_seqs(reference)
    cons_seqs = read_fasta_seqs(consensus)
    cons_ncounts = read_fasta_ncounts(consensus)

    # Positional fallback: when reference and consensus use different seqids
    # (e.g. mutant FASTA uses a synthetic id while the assembler writes the
    # original accession), and the reference has exactly one segment, map the
    # sole consensus sequence onto it regardless of name.
    _cons_ids = list(cons_seqs)
    _single_fallback = cons_seqs[_cons_ids[0]] if len(_cons_ids) == 1 else None
    _ncounts_fallback = cons_ncounts[_cons_ids[0]] if len(_cons_ids) == 1 else None

    rows: list[dict] = []
    for seg, rlen in ref_lengths.items():
        # Prefer exact seqid match; fall back to sole-sequence positional match.
        cons_seq = cons_seqs.get(seg) or (
            _single_fallback if len(ref_lengths) == 1 else None
        )
        nc_entry = cons_ncounts.get(seg) or (
            _ncounts_fallback if len(ref_lengths) == 1 else None
        )
        clen, ncount = nc_entry if nc_entry else (0, 0)
        called = clen - ncount
        pos_pct = _positions_called(rlen, called)

        if not cons_seq:
            rows.append({
                "segment": seg, "ref_len": rlen, "covered_bp": 0,
                "aligned_bp": 0, "matching_bp": 0,
                "coverage_pct": 0.0, "identity_pct": 0.0,
                "called_bp": called, "positions_called_pct": pos_pct,
                "insertions": 0, "deletions": 0,
            })
            continue
        aln_ref, aln_cons = align_pair(ref_seqs[seg], cons_seq, num_threads)
        covered, aligned, matching = identity_counts(aln_ref, aln_cons)
        insertions, deletions = indel_counts(aln_ref, aln_cons)
        rows.append({
            "segment": seg, "ref_len": rlen, "covered_bp": covered,
            "aligned_bp": aligned, "matching_bp": matching,
            "coverage_pct": round(covered / rlen * 100, 2) if rlen else 0.0,
            "identity_pct": round(matching / aligned * 100, 2) if aligned else 0.0,
            "called_bp": called, "positions_called_pct": pos_pct,
            "insertions": insertions, "deletions": deletions,
        })
    return pd.DataFrame(rows, columns=COVERAGE_COLUMNS)


def empty_coverage_table(reference: str) -> pd.DataFrame:
    """Coverage table with NA metrics — used when the consensus is absent/empty.

    With no consensus to read, ``called_bp``/``positions_called_pct`` are ``None``
    alongside the alignment metrics.
    """
    ref_lengths = read_fasta_lengths(reference)
    rows = [{
        "segment": seg, "ref_len": rlen, "covered_bp": 0,
        "aligned_bp": 0, "matching_bp": 0,
        "coverage_pct": None, "identity_pct": None,
        "called_bp": None, "positions_called_pct": None,
        "insertions": None, "deletions": None,
    } for seg, rlen in ref_lengths.items()]
    return pd.DataFrame(rows, columns=COVERAGE_COLUMNS)


def summarise(cov_df: pd.DataFrame) -> dict:
    """Collapse a per-segment table to one summary row.

    `mean_*` is the unweighted across-segment mean (each segment counts equally);
    `weighted_*` is length-weighted (the genome-wide figure). They coincide for a
    single-segment genome. NA per-segment metrics (no consensus) propagate to NA.

    `positions_called_pct` columns are populated when ``called_bp`` is present
    (i.e. a consensus was read); otherwise both summary columns are ``None``.
    """
    n = len(cov_df)
    total_ref = int(cov_df["ref_len"].sum()) if n else 0
    total_cov = int(cov_df["covered_bp"].sum()) if n else 0
    total_aln = int(cov_df["aligned_bp"].sum()) if n else 0
    total_match = int(cov_df["matching_bp"].sum()) if n else 0

    has_metrics = n > 0 and cov_df["coverage_pct"].notna().any()
    has_positions = (
        n > 0 and
        "called_bp" in cov_df.columns and
        cov_df["called_bp"].notna().any()
    )
    total_called = int(cov_df["called_bp"].sum()) if has_positions else None
    has_indels = (
        has_metrics and
        {"insertions", "deletions"} <= set(cov_df.columns) and
        cov_df["insertions"].notna().any()
    )

    return {
        "n_segments": n,
        "total_ref_len": total_ref,
        "mean_coverage_pct": round(cov_df["coverage_pct"].mean(), 2) if has_metrics else None,
        "mean_identity_pct": round(cov_df["identity_pct"].mean(), 2) if has_metrics else None,
        "weighted_coverage_pct": round(total_cov / total_ref * 100, 2)
        if (has_metrics and total_ref) else None,
        "weighted_identity_pct": round(total_match / total_aln * 100, 2)
        if (has_metrics and total_aln) else None,
        "mean_positions_called_pct": round(cov_df["positions_called_pct"].mean(), 2)
        if has_positions else None,
        "weighted_positions_called_pct": round(total_called / total_ref * 100, 2)
        if (has_positions and total_ref) else None,
        # error counts over all segments (events; substitutions per aligned base)
        "substitutions": (total_aln - total_match) if has_metrics else None,
        "insertions": int(cov_df["insertions"].sum()) if has_indels else None,
        "deletions": int(cov_df["deletions"].sum()) if has_indels else None,
    }
