#!/usr/bin/env python3
"""Score the coverage-over-time assemblies of framework_C.coverage_over_time.bash.

For every assembly point <lib>/coverage_over_time/b<NNNN>[.medaka]/consensus.fa,
each reference segment is scored against the consensus record of the same name
by ngs-benchmark (`ngs_benchmark.assembly.metrics.coverage_table`: MAFFT G-INS-i,
--adjustdirection, one alignment at a time, indels as gap runs), and joined to
the batch's sequencing time and cumulative reads from <lib>/timeseries.tsv.

One row per library x batch x segment; whole-genome coverage is
sum(covered_bp) / sum(ref_len) over the 10 segments (done in the R figure).
A segment with no consensus record scores covered_bp 0.

    python3 jobs/framework_C/framework_C.coverage_over_time.score.py \\
        REF.fasta LIB [LIB ...] > metrics/framework_C/coverage_over_time.tsv
"""
import csv
import glob
import re
import sys

from ngs_benchmark.assembly.metrics import coverage_table

FC = "results/framework_C"


def main():
    ref = sys.argv[1]
    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(["library", "batch", "polished", "seq_time_s", "reads_cum", "segment",
                  "ref_len", "covered_bp", "subs", "ins", "del", "identity"])
    for lib in sys.argv[2:]:
        timing = {}
        with open(f"{FC}/{lib}/timeseries.tsv") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                timing[int(row["batch"])] = (row["seq_time_s"], row["reads_cum"])
        for cons in sorted(glob.glob(f"{FC}/{lib}/coverage_over_time/b*/consensus.fa")):
            m = re.search(r"/b(\d+)(\.medaka)?/consensus\.fa$", cons)
            batch, polished = int(m[1]), bool(m[2])
            seq_time, reads_cum = timing[batch]
            print(f"[{lib} b{batch}{' polished' if polished else ''}]", file=sys.stderr, flush=True)
            for r in coverage_table(cons, ref).itertuples():
                ident = 100 * r.matching_bp / r.aligned_bp if r.aligned_bp else float("nan")
                out.writerow([lib, batch, int(polished), seq_time, reads_cum, r.segment,
                              r.ref_len, r.covered_bp, r.aligned_bp - r.matching_bp,
                              r.insertions, r.deletions, f"{ident:.3f}"])


if __name__ == "__main__":
    main()
