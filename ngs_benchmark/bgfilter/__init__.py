"""Module D — background-filtering efficiency (LazypipeX Table 1 analogue).

Scores how a host filter partitions reads by true origin: the fraction of reads
removed per category (Homo sapiens, Bacteria, Viruses, HERV, Viruses excl.
Retroviridae) at each filter setting. Public entry points:

    from ngs_benchmark.bgfilter.truth import load_truth, count_data_lines
    from ngs_benchmark.bgfilter.runner import BgfilterConfig, run, run_batch
"""
