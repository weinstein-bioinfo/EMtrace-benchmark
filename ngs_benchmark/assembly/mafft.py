"""Pairwise nucleotide identity via MAFFT — ported/adapted from virefset.

Standalone: no virefset or Biopython dependency. A consensus is compared to its
reference by a global pairwise MAFFT alignment (ginsi), and identity is computed
over confidently-called columns (gaps and consensus N positions excluded). This
replaces minimap2 for consensus↔reference scoring, which under-reported identity
because PAF column 10 without ``-c`` is a seed-based estimate, not true identity.

Both coverage and identity exclude consensus N positions (N treated as a gap):
``covered_bp`` counts only non-gap, non-N consensus positions; ``aligned_bp``
likewise excludes N from the identity denominator, so identity reflects accuracy
of confidently-called bases only (completeness reported separately via
``positions_called_pct``).  The identity logic follows virefset
``steps/distmat._calc_identity``; the wrapper
follows ``virefset/wrappers/mafft.py`` (ginsi = ``--globalpair --maxiterate 1000``),
trimmed to the single-pair case needed here.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

#: ginsi — global pairwise with iterative refinement (matches virefset genome mode).
GINSI_ARGS = ("--globalpair", "--maxiterate", "1000")

#: Reverse-complement the consensus when that aligns better to the reference
#: (written first, so never flipped). A consensus can sit on the other strand,
#: e.g. an ambisense Lassa S segment assembled on a RefSeq record stored in the
#: opposite orientation; aligned as given it would score ~40% identity. MAFFT
#: renames a flipped record with an ``_R_`` prefix (see `_parse_two`).
#: The fast 6-mer check (not --adjustdirectionaccurately): on 15 Framework B
#: pairs, reverse-strand LASV_S and 25-30%-divergent forward genomes alike, both
#: gave identical results; the accurate mode about doubled the time per pair.
DIRECTION_ARGS = ("--adjustdirection",)


def have_mafft() -> bool:
    """True if the ``mafft`` binary is on PATH."""
    return shutil.which("mafft") is not None


def _mafft_exe() -> str:
    exe = shutil.which("mafft")
    if exe is None:
        raise RuntimeError(
            "mafft not found on PATH; it is required for assembly identity "
            "(consensus↔reference alignment). "
            "Install: https://mafft.cbrc.jp/alignment/software/"
        )
    return exe


def _parse_two(fasta_text: str) -> tuple[str, str]:
    """Parse MAFFT FASTA stdout (records ``__ref__`` and ``__cons__``).

    A record reverse-complemented by ``--adjustdirection*`` comes back as
    ``_R___cons__``; it is read as the consensus.
    """
    seqs: dict[str, list[str]] = {}
    name: str | None = None
    for line in fasta_text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0] if len(line) > 1 else ""
            if name.startswith("_R_"):
                name = name[3:]
            seqs[name] = []
        elif name is not None:
            seqs[name].append(line.strip())
    ref = "".join(seqs.get("__ref__", [])).upper()
    cons = "".join(seqs.get("__cons__", [])).upper()
    return ref, cons


def align_pair(ref_seq: str, cons_seq: str, num_threads: int = 1) -> tuple[str, str]:
    """Globally align two nucleotide sequences with MAFFT (ginsi).

    Returns ``(aligned_ref, aligned_cons)`` — gapped, uppercased, equal length.
    The records get distinct ids (``__ref__`` / ``__cons__``) so they can be told
    apart even when the consensus inherits the reference accession as its name.
    """
    exe = _mafft_exe()
    with tempfile.TemporaryDirectory() as td:
        in_fa = Path(td) / "pair.fa"
        in_fa.write_text(f">__ref__\n{ref_seq}\n>__cons__\n{cons_seq}\n")
        cmd = [exe, "--quiet", *GINSI_ARGS, *DIRECTION_ARGS,
               "--thread", str(num_threads), str(in_fa)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"MAFFT failed (exit {proc.returncode}):\n{proc.stderr.strip()}"
        )
    return _parse_two(proc.stdout)


def identity_counts(aligned_ref: str, aligned_cons: str) -> tuple[int, int, int]:
    """Coverage + identity counts over a pairwise alignment.

    Adapted from virefset ``distmat._calc_identity``, with consensus N positions
    treated like gaps (excluded from numerator and denominator) so identity
    reflects accuracy of confidently-called bases only — completeness is reported
    separately as positions_called.

    Returns ``(covered_bp, aligned_bp, matching_bp)``:
        covered_bp  — reference positions aligned to a non-gap, non-N consensus
                      base: reference positions confidently called by the assembly.
        aligned_bp  — columns where both sequences are non-gap and non-N: the
                      identity denominator.
        matching_bp — aligned_bp columns where the two bases are equal.
    """
    covered = aligned = matching = 0
    for cr, cc in zip(aligned_ref, aligned_cons):
        ref_gap = cr == "-"
        cons_gap = cc == "-"
        if not ref_gap and not cons_gap and cc not in "Nn":
            covered += 1
        if ref_gap or cons_gap or cr in "Nn" or cc in "Nn":
            continue
        aligned += 1
        if cr == cc:
            matching += 1
    return covered, aligned, matching


def indel_counts(aligned_ref: str, aligned_cons: str) -> tuple[int, int]:
    """Insertions and deletions of the consensus relative to the reference.

    Read from the same pairwise alignment as `identity_counts`:

        insertion  — a run of gap columns in the reference (opposite consensus
                     bases); N bases inside it do not split it, and a run made
                     only of N is not counted
        deletion   — a run of gap columns in the consensus (opposite reference
                     bases); a consensus N column ends it

    Each run counts once, whatever its length. Only columns between the first
    and the last called (non-gap, non-N) consensus base are examined, and
    consensus N columns are skipped, so masked (uncalled) regions and the
    unaligned ends of a partial consensus never count as deletions.

    Returns ``(insertions, deletions)``.
    """
    ref, cons = aligned_ref.upper(), aligned_cons.upper()
    called = [i for i, cc in enumerate(cons) if cc not in "-N"]
    if not called:
        return 0, 0
    lo, hi = called[0], called[-1]
    insertions = deletions = 0
    i = lo
    while i <= hi:
        if cons[i] == "N":
            i += 1
        elif ref[i] == "-":                 # insertion: the whole reference-gap run
            insertions += 1
            while i <= hi and ref[i] == "-":
                i += 1
        elif cons[i] == "-":                # deletion: the consensus-gap run
            deletions += 1
            while i <= hi and cons[i] == "-":
                i += 1
        else:
            i += 1
    return insertions, deletions
