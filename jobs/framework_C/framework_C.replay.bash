#!/usr/bin/env bash
# ============================================================================
# Framework C — real-time replay of MinKNOW runs through `emtrace watch`
# ============================================================================
#
# Design: the project notes. Each library's own MinKNOW batch files
# are released one at a time into a watch directory (accelerated: the next file
# goes out as soon as the previous batch is reported), and Module C
# (`ngs-benchmark realtime`) scores every batch:
#   timeseries.tsv  first_call.tsv  latency.tsv  accuracy.tsv
#
# Runs LOCALLY on the laptop, by design: real-time profiling happens next to the
# sequencer (the project notes). Run from the repo root:
#   caffeinate -ims bash jobs/framework_C/framework_C.replay.bash
#
# Hit mode: ALL HITS (EMtrace default). Both configs set
# annotation.filter_tophits: false explicitly, and the script refuses a config
# that turns it on.
#
# Run conditions recorded per library in <outdir>/provenance.txt:
#   host, CPU, RAM, threads, config, anns, hit mode, git commits, tool versions.
#
# Orbivirus (Orbivirus_kevo; config emtrace.config.frameC.orbivirus.yaml):
#   the published viral library. Kevo orbivirus (KEVOOV, isolate fi663), FLAC
#   library from a C6/36 isolate, R9.4.1 / LSK110, 464,000 reads over 19.0 h
#   (VZRU; doi 10.1099/jgv.0.002101). Host filter kitome + Aedes albopictus + Homo sapiens.
#   NOVEL-VIRUS CASE: KEVOOV (PQ165789-98) is not in minimap.gbvi_abv1g, whose
#   closest relatives sit at ~71 % nt, so the call is scored at genus
#   (Orbivirus). Three cases, each in its own output dir, all kept, all on a
#   prebuilt index (the project notes):
#     RUN_ORBI_RefSeq    Orbivirus_kevo.refseq_abv1g_k11w5  minimap.refseq.abv.1g.k11w5
#                        (prebuilt k11 w5 index, 29 GB)
#     RUN_ORBI_RefSeqDmnd Orbivirus_kevo.refseq_abv1g_k11w5_dmnd
#                        minimap.refseq.abv.1g.k11w5,diamond.uniref.vi
#                        (the same, chained: reads minimap2 leaves unassigned go
#                        to DIAMOND blastp on UniRef100 viral, bitscore >= 80,
#                        ORFs >= 90 nt: the manuscript's chain2 settings.
#                        2026/10/22 test on batches 1-5 at bitscore 100 / ORF
#                        100: +87 % Orbivirus reads, 3,177 vs 1,697)
#     RUN_ORBI_GBRefSeq  Orbivirus_kevo.gbvi_abv1g          minimap.gbvi_abv1g
#                        (GenBank vi + RefSeq abv 1g, prebuilt index at minimap2's
#                        default k15 w10, 37.6 GB: the larger database, no
#                        re-indexing per batch; asks whether the CHeRI relatives
#                        it adds outweigh the larger seeds in the first batches)
#   2026/10/22: the GenBank case previously ran on the FASTA at k11 w5
#   (Orbivirus_kevo.gbvi_abv1g_k11w5, minimap.gbvi_abv1g.fa); minimap2 re-indexed
#   it on every batch (~18 min/batch), so it was stopped after 36 of 116 batches.
#   Per case the script derives a config whose anns strategy is
#   `--ann1 <db>`, because Module C only passes --anns through. All cases
#   share one set of synthetic batch files.
#   SYNTHETIC BATCHES: the library is archived as one concatenated FASTQ, so
#   stage_synthetic re-cuts it into SYNTH_BATCH_READS-read files sorted by
#   start_time, keeping the MinKNOW headers (as for the gut mock, below).
#
# TBEV (TBEV_FBC92293; config emtrace.config.frameC.tbev.yaml): published
#   in-house viral library (VZRU). One run, flow cell
#   FBC92293, library tbev_teemu_231026 (23 Oct 2026), R10.4.1, Dorado
#   hac@v5.2.0, pass reads; 08:21-13:20 +03:00, 4.98 h, largest gap 168 s.
#   Barcodes 10, 12, 13, 14 (IL209-11-15, IL232-5, IL239-1, IL256-1; 41-56 %
#   TBEV in their first 5,000 reads) are concatenated into one library of
#   313,867 reads; barcode 11 (IL232-4, 1 TBEV read in 5,000) is dropped.
#   Preprocessing for these short reads (N50 177-299 bp):
#   porechop_abi no_split true, chopper Q8 / 100 bp.
#   Dorado left the native-barcode adapters on, so porechop_abi trims ~118 bp
#   from nearly every read start and drops ~35 % of reads (81-270 bp raw,
#   essentially adapter only). Host filter kitome + Ixodes ricinus + Ixodes
#   persulcatus + Homo sapiens, annotation minimap.refseq.abv.1g (prebuilt k15
#   index; TBEV NC_001672 is in it), scored at species.
#   SYNTHETIC BATCHES: Dorado st:Z: tags, so stage_synthetic sorts and re-cuts
#   it like the gut mock (~79 batches of 4,000).
#   TBEV_FIRST_BATCH=true replays batch 0 alone (TBEV_FBC92293.refseq_abv1g.first_batch)
#   to time the four-host bgfilter before committing to the full replay.
#
# Ticks (ticks_<pool>; config emtrace.config.frameC.ticks.yaml): the PUBLIC
#   FIELD library. PRJNA878809, Ergunay et al. 2023 Front Microbiol 14:1177651
#   (doi 10.3389/fmicb.2023.1177651): 40 pools of 301 field-collected ticks from
#   Anatolia, Turkiye, 2020-2021; untargeted random-hexamer cDNA on a GridION
#   Mk1C. Unlike every other viral library here the virus is NOT amplified, so the
#   load is 0.1-0.5 % of reads and the call has to accrue over batches.
#   Pools (CCHFV status from the paper's supplementary Table S3,
#   LIT/NanoporeMetagenomics/Ergunay-2023-supp/Table 1.XLS):
#     80-43  SRR21504608  303,354 reads   PCR-NEGATIVE, 4 CCHFV reads in the paper.
#                                         Replayed 28 Oct 2026: 3 Orthonairovirus reads,
#                                         no call. A negative case.
#     35-33  SRR21504592  216,828 reads   PCR-positive, 783 CCHFV reads; L and S covered
#     58-37  SRR21504617   44,503 reads   PCR-positive, 1,752 CCHFV reads; L, M and S
#                                         covered (the only pool with all three segments)
#     32-31  SRR21504594  820,288 reads   PCR-positive, 914 CCHFV reads
#   (The six PCR-positive pools: 58-37, 32-31, 35-33, 62-38, 66-41, 46-35.)
#   TRUTH TAXON IS THE GENUS Orthonairovirus, not CCHFV at species. The paper's
#   sequences cluster with Aigai virus / genogroup VI (AP92-like), and NCBI files
#   Aigai as its own species (Orthonairovirus parahaemorrhagiae) beside CCHFV
#   (O. haemorrhagiae), so the reads split between two species exactly as the Kevo
#   orbivirus reads split between relatives. Scored at genus, as Orbivirus_kevo is.
#   Two annotation variants, as for Orbivirus_kevo:
#     RUN_TICKS_RefSeq     ticks_<pool>.refseq_abv1g_k11w5       abv.1g.k11w5
#     RUN_TICKS_RefSeqDmnd ticks_<pool>.refseq_abv1g_k11w5_dmnd  abv.1g.k11w5.dmnd
#                          (minimap2 k11 w5, then DIAMOND blastp on UniRef100 viral
#                          for the unassigned reads; bitscore >= 80, ORFs >= 90 nt)
#   SMALL SEEDS ARE REQUIRED: on pool 35-33, k11 w5 recovers 416 CCHFV reads
#   against 55 at the k15 default (7.6x), the Kevo effect again.
#   NO TIMESTAMPS: SRA strips start_time AND the read IDs (deflines are
#   `@SRR21504608.1 1/1`), so stage_ordered cuts the library in file order and the
#   replay has NO sequencing-time axis - read count is the only x-axis. Module C
#   runs with --reads-only: seq_time_s, seq_time_covered_s and the throughput
#   ratio are empty; latency per batch and per 1k reads are still measured.
#   SCORED WITH A PANEL, NOT A TRUTH (--panel panels/ticks_<pool>.panel.tsv): a
#   field library has no complete truth, so Module C computes no Pr/Rc/FP. Each
#   taxon the paper reports is scored detected yes/no (>= 10 reads, persistent
#   first detection, whatever its rank) by its evidence: pcr_pos / pcr_neg (the
#   in-house generic-nairovirus and JMTV PCRs, Tables S2-S3) and reported (the
#   paper's BLAST, Tables S4-S5; agreement only). Outputs: panel.tsv and
#   panel_summary.tsv ("PCR targets detected k/n", "PCR negatives clear k/n").
#   The ENA FASTQ is gzip-tested after download (cf. the corrupt MPXV archive).
#   READS ARRIVE PREPROCESSED: the deposit is post-Porechop and post-NanoFilt
#   (q <= 9, length <= 100 bp removed), host reads still in. Our preprocessing is
#   therefore near a no-op and the read funnel is NOT comparable with the in-house
#   libraries.
#
# 16S mock (barcode79): the configuration frozen for 16S benchmark Part 2
#   (the project notes): --config emtrace.config.16S.yaml,
#   --anns 16S.rrnDB.dedup, minimap2 -N 25, em.temp 10, all hits.
#   `emtrace watch` has no -N option, so the script derives a config from the
#   16S one with annotation.minimap.max_target_seqs 25, and with
#   preprocess.bgfilter.bgfilter '0': watch falls back to the config's host list
#   when --hostgen is absent, and Part 1 ran without host filtering on the mock.
#   Scored at genus against the 16S-only truth, with the 10-read cutoff and the
#   profile weighted on readn_em (EM-expected reads), as in the 16S benchmark.
#
# 16S gut mock (gutmock_sup_r1..r4): ONT open data zymo_16s_2025.09, R10.4.1,
#   ZymoBIOMICS Gut Microbiome Standard, four replicates (two barcodes on each of
#   two flow cells, ~2 h of sequencing each), replayed to measure how reliable
#   time-to-stabilisation is across replicates. Same derived config and anns as
#   barcode79 (the frozen 16S settings; the R10.4.1 cells also ran without host
#   filtering), scored at genus against gutmock_16S_truth.tsv.
#   SYNTHETIC BATCHES: ONT publishes one re-basecalled BAM per replicate, not
#   the MinKNOW batch files, so stage_synthetic re-cuts each replicate's FASTQ
#   (data/sra/16S/<library>.fastq.gz, from polish.16S.comparators.bash
#   --part fetch) into MinKNOW-style files: reads sorted by Dorado's st (start
#   time) tag, SYNTH_BATCH_READS per file (MinKNOW's default 4,000), header
#   rewritten to `start_time=`. Unlike the other libraries, batch boundaries are
#   therefore the benchmark's, not MinKNOW's — state this beside the results.
#
# Run:
#   cd EMtrace-benchmark
#   caffeinate -ims bash jobs/framework_C/framework_C.replay.bash
# ============================================================================

set -uo pipefail

# =============================================================================
# Control switches
# =============================================================================
RUN_ORBI_RefSeq=false   # Orbivirus_kevo vs RefSeq abv 1g k11w5 index       116 batches, 464,000 reads
RUN_ORBI_RefSeqDmnd=false # Orbivirus_kevo vs RefSeq abv 1g k11w5 index -> DIAMOND uniref.vi (bits 80, ORF 90; chain2)
RUN_ORBI_GBRefSeq=false # Orbivirus_kevo vs GenBank vi + RefSeq abv 1g prebuilt k15 index
RUN_TBEV=false           # TBEV_FBC92293 vs RefSeq abv 1g k15 index         ~79 batches, 313,867 reads
TBEV_FIRST_BATCH=false   # true: replay batch 0 only, to time the four-host bgfilter
RUN_TICKS_RefSeq=true     # ticks vs RefSeq abv 1g k11w5 index              PUBLIC FIELD library
RUN_TICKS_RefSeqDmnd=false  # ticks vs RefSeq abv 1g k11w5 -> DIAMOND uniref.vi (chain2)
TICK_POOLS="35-33 58-37"         # space separated: 80-43 (PCR-negative); 35-33, 58-37, 32-31 (PCR-positive)
RUN_BARCODE79=false   # 16S mock barcode79      26 batches, 100,128 reads (frozen Part 2 settings)
RUN_GUTMOCK=false     # 16S gut mock R10.4.1    ~18-22 batches, 72,830-84,233 reads per replicate
GUTMOCK_REPS="1 2 3 4"

STAGE_HOSTS=true      # download + index missing host genomes before replaying
FORCE=false           # true: re-run libraries that already have timeseries.tsv
# =============================================================================

THREADS=8
DATA="data"
OUT_ROOT="results/framework_C"
LOG_DIR="logs/framework_C"
ORBI_LIB="Orbivirus_kevo"
ORBI_FQ="$DATA/fastq/Orbivirus_kevo.fastq.gz"   # reads available from the authors on request
ORBI_BATCH_DIR="$DATA/orbivirus.kevo"           # synthetic batch dir <lib>/ is staged here
ORBI_BASE_CFG="jobs/framework_C/emtrace.config.frameC.orbivirus.yaml"
ORBI_ANNS="frameC.orbi"                         # strategy written into each derived config
# Cases: run name (output dir + sample), --ann1 db key under databases: in ORBI_BASE_CFG
ORBI_REFSEQ_RUN="$ORBI_LIB.refseq_abv1g_k11w5";  ORBI_REFSEQ_ANN1="minimap.refseq.abv.1g.k11w5"
ORBI_DMND_RUN="$ORBI_LIB.refseq_abv1g_k11w5_dmnd"; ORBI_DMND_ANN1="minimap.refseq.abv.1g.k11w5,diamond.uniref.vi"
ORBI_DMND_SET=(annotation.diamond_blastp.bitscore=80 annotation.orfipy.min_len=90)   # manuscript chain2
ORBI_GB_RUN="$ORBI_LIB.gbvi_abv1g";              ORBI_GB_ANN1="minimap.gbvi_abv1g"
ORBI_ANY=$([[ "$RUN_ORBI_RefSeq" == true || "$RUN_ORBI_RefSeqDmnd" == true || "$RUN_ORBI_GBRefSeq" == true ]] && echo true || echo false)
ORBI_TAXON="Orbivirus"                          # genus: KEVOOV itself is not in the db
TBEV_LIB="TBEV_FBC92293"
TBEV_BARCODES="IL209-11-15 IL232-5 IL239-1 IL256-1"   # IL232-4 dropped (no TBEV)
TBEV_FQ_DIR="$DATA/fastq/TBEV"                  # reads available from the authors on request
TBEV_FQ="$DATA/fastq/$TBEV_LIB.fastq.gz"        # the barcodes concatenated (built below)
TBEV_BATCH_DIR="$DATA/tbev.inhouse"             # synthetic batch dir <lib>/ is staged here
TBEV_CFG="jobs/framework_C/emtrace.config.frameC.tbev.yaml"
TBEV_ANNS="abv.1g"                              # minimap.refseq.abv.1g: RefSeq abv 1-genome (23 GB, k15)
TBEV_RUN="$TBEV_LIB.refseq_abv1g"
TBEV_TAXON="Orthoflavivirus encephalitidis"     # TBEV, taxid 11084
TICK_CFG="jobs/framework_C/emtrace.config.frameC.ticks.yaml"
TICK_FQ_DIR="$DATA/fastq/ticks.ergunay"     # <SRR>.fastq.gz, fetched from ENA
TICK_BATCH_DIR="$DATA/ticks.ergunay"        # synthetic batch dir <lib>/ is staged here
TICK_REFSEQ_ANNS="abv.1g.k11w5"             # minimap.refseq.abv.1g.k11w5 (prebuilt, 29 GB)
TICK_DMND_ANNS="abv.1g.k11w5.dmnd"          # the same, chained to DIAMOND uniref.vi
TICK_PANELS="jobs/framework_C/panels" # ticks_<pool>.panel.tsv, one per pool; GENUS
                                            # rank: the CCHFV is Aigai-like, so its reads
                                            # split between O. haemorrhagiae and
                                            # O. parahaemorrhagiae (see header)
TICK_ANY=$([[ "$RUN_TICKS_RefSeq" == true || "$RUN_TICKS_RefSeqDmnd" == true ]] && echo true || echo false)
# tick_srr <pool> -> SRA run accession
tick_srr() {
  case "$1" in
    80-43) echo SRR21504608 ;;   # 303,354 reads, CCHFV PCR-NEGATIVE (4 reads, Table S3)
    35-33) echo SRR21504592 ;;   # 216,828 reads, CCHFV PCR-positive (783 reads; L, S)
    58-37) echo SRR21504617 ;;   #  44,503 reads, CCHFV PCR-positive (1,752 reads; L, M, S)
    32-31) echo SRR21504594 ;;   # 820,288 reads, CCHFV PCR-positive (914 reads)
    *) echo "" ;;
  esac
}
S16_BASE_CFG="jobs/framework_16S/emtrace.config.16S.yaml"
S16_CFG="$LOG_DIR/emtrace.config.16S.frameC.N25.yaml"   # derived below
S16_N=25                    # minimap2 -N (annotation.minimap.max_target_seqs), frozen in 16S Part 2
S16_ANNS="16S.rrnDB.dedup"
S16_TRUTH="data/truth/zymobiomics_16S_truth.tsv"
GUT_TRUTH="data/truth/gutmock_16S_truth.tsv"
GUT_MODEL="sup"
GUT_FQ_DIR="data/sra/16S"            # <library>.fastq.gz, st:Z: tags kept
GUT_BATCH_DIR="$DATA/gutmock.16S"          # synthetic batch dirs, one per replicate
SYNTH_BATCH_READS=4000                     # MinKNOW's default reads per file (gut mock, Orbivirus)
MIN_FREE_GB=25        # refuse to start below this much free disk

# Per-batch timeout: generous, because each batch reloads ~40 GB of indexes.
BATCH_TIMEOUT=7200

mkdir -p "$OUT_ROOT" "$LOG_DIR"

# --- Preconditions ------------------------------------------------------------
[[ -f jobs/framework_C/framework_C.replay.bash ]] || { echo "Run from the EMtrace-benchmark repo root." >&2; exit 1; }
: "${EMTRACE_DB_DIR:?EMTRACE_DB_DIR must point at the EMtrace databases}"

for tool in emtrace ngs-benchmark minimap2 chopper; do
  command -v "$tool" > /dev/null || { echo "Missing tool: $tool" >&2; exit 1; }
done
# Python for the inline helpers: the interpreter of the emtrace launcher, so
# `import emtrace` (and yaml) work whatever python is first on PATH.
RT_PY=$(head -1 "$(command -v emtrace)" | sed -E 's/^#! *//; s/^\/usr\/bin\/env +//')
if ! "$RT_PY" -c 'import emtrace, yaml' 2> /dev/null; then
  echo "Cannot import emtrace + yaml with $RT_PY (from the emtrace shebang)." >&2
  exit 1
fi
if ! ngs-benchmark realtime --help | grep -q -- "--config"; then
  echo "ngs-benchmark is too old: its realtime command has no --config (driver fix missing)." >&2
  exit 1
fi
if [[ "$TICK_ANY" == true ]] && ! ngs-benchmark realtime --help | grep -q -- "--panel"; then
  echo "ngs-benchmark is too old: its realtime command has no --panel (tick libraries)." >&2
  exit 1
fi
if ! emtrace run --help | grep -q -- "--tophits"; then
  echo "EMtrace is too old: all-hits annotation (--tophits) is missing." >&2
  exit 1
fi

free_gb=$(df -g "$EMTRACE_DB_DIR" | awk 'NR==2 {print $4}')
if (( free_gb < MIN_FREE_GB )); then
  echo "Only ${free_gb} GB free on the data volume (need ${MIN_FREE_GB})." >&2
  exit 1
fi

# Derive the 16S replay config: the base 16S config with -N 25 and no host filter.
if [[ "$RUN_BARCODE79" == true || "$RUN_GUTMOCK" == true ]]; then
  "$RT_PY" - "$S16_BASE_CFG" "$S16_CFG" "$S16_N" <<'CFGPY' || { echo "[FATAL] cannot derive $S16_CFG" >&2; exit 1; }
import sys, yaml
src, dst, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
cfg = yaml.safe_load(open(src))
cfg["annotation"]["minimap"]["max_target_seqs"] = n
cfg["preprocess"]["bgfilter"]["bgfilter"] = "0"
with open(dst, "w") as fh:
    fh.write(f"# Derived by jobs/framework_C/framework_C.replay.bash from {src}:\n")
    fh.write(f"#   annotation.minimap.max_target_seqs {n}; preprocess.bgfilter.bgfilter '0'\n")
    yaml.safe_dump(cfg, fh, sort_keys=False, default_flow_style=False)
CFGPY
fi

# Derive an Orbivirus replay config: the base config with one anns strategy,
# `--ann1 <db>`, so the database is set here rather than in the YAML. Seeding
# (annotation.minimap k 11 / w 5) is in the base YAML. Optional <path>=<value>
# arguments override existing keys (dotted path from the YAML root).
# derive_orbi_cfg <ann1> <dst> [path=value ...]
derive_orbi_cfg() {
  "$RT_PY" - "$ORBI_BASE_CFG" "$2" "$ORBI_ANNS" "$1" "${@:3}" <<'CFGPY' || { echo "[FATAL] cannot derive $2" >&2; exit 1; }
import sys, yaml
src, dst, key, ann1 = sys.argv[1:5]
cfg = yaml.safe_load(open(src))
missing = [db for db in ann1.split(",") if db not in cfg["databases"]]
if missing:
    sys.exit(f"--ann1 {ann1}: not under databases: in {src}: {', '.join(missing)}")
cfg["annotation"]["anns.options"][key] = f"--ann1 {ann1}"
cfg["annotation"]["anns"] = key
for kv in sys.argv[5:]:
    path, val = kv.split("=", 1)
    *parents, leaf = path.split(".")
    node = cfg
    for p in parents:
        node = node[p]
    if leaf not in node:
        sys.exit(f"{path}: no such key in {src}")
    node[leaf] = yaml.safe_load(val)
with open(dst, "w") as fh:
    fh.write(f"# Derived by jobs/framework_C/framework_C.replay.bash from {src}:\n")
    fh.write(f"#   annotation.anns {key} = '--ann1 {ann1}'\n")
    for kv in sys.argv[5:]:
        fh.write(f"#   {kv}\n")
    yaml.safe_dump(cfg, fh, sort_keys=False, default_flow_style=False)
CFGPY
}
[[ "$RUN_ORBI_RefSeq" == true ]] && derive_orbi_cfg "$ORBI_REFSEQ_ANN1" "$LOG_DIR/emtrace.config.$ORBI_REFSEQ_RUN.yaml"
[[ "$RUN_ORBI_RefSeqDmnd" == true ]] && derive_orbi_cfg "$ORBI_DMND_ANN1" "$LOG_DIR/emtrace.config.$ORBI_DMND_RUN.yaml" "${ORBI_DMND_SET[@]}"
[[ "$RUN_ORBI_GBRefSeq" == true ]] && derive_orbi_cfg "$ORBI_GB_ANN1" "$LOG_DIR/emtrace.config.$ORBI_GB_RUN.yaml"

# The GenBank case runs on the .fa.gz; EMtrace looks its taxids up as
# <db>.acc2taxid, so the symlink to the .fa.acc2taxid must exist (the project notes).
if [[ "$RUN_ORBI_GBRefSeq" == true ]]; then
  gb_fa="$EMTRACE_DB_DIR/minimap.genbank_vi+refseq_abv1g.fa.gz"
  [[ -s "$gb_fa" ]] || { echo "[FATAL] $gb_fa missing" >&2; exit 1; }
  [[ -s "$gb_fa.acc2taxid" ]] || { echo "[FATAL] $gb_fa.acc2taxid missing — ln -s minimap.genbank_vi+refseq_abv1g.fa.acc2taxid $gb_fa.acc2taxid" >&2; exit 1; }
fi
if [[ "$TICK_ANY" == true ]]; then
  k11_mmi="$EMTRACE_DB_DIR/minimap.refseq.abv.1g.k11w5.fa.mmi"
  [[ -s "$k11_mmi" ]] || { echo "[FATAL] $k11_mmi missing" >&2; exit 1; }
fi
if [[ "$RUN_ORBI_RefSeqDmnd" == true || "$RUN_TICKS_RefSeqDmnd" == true ]]; then
  dmnd_db="$EMTRACE_DB_DIR/diamond.uniref100.vi.dmnd"
  [[ -s "$dmnd_db" ]] || { echo "[FATAL] $dmnd_db missing" >&2; exit 1; }
fi

# Refuse a config that switches top-hit filtering on: this replay is all-hits.
for cfg in "$S16_BASE_CFG" "$ORBI_BASE_CFG" "$TBEV_CFG" "$TICK_CFG"; do
  if grep -Eq '^\s*filter_tophits:\s*true' "$cfg"; then
    echo "$cfg sets annotation.filter_tophits: true — this replay is all-hits." >&2
    exit 1
  fi
done

# Refuse a config that trims the abundance tail in the pipeline. Scoring
# applies its own tail (framework_C.barcode79.ranks.py, framework_16S.score.py);
# trimming twice would silently change what the published numbers mean.
for cfg in "$S16_BASE_CFG" "$ORBI_BASE_CFG" "$TBEV_CFG" "$TICK_CFG"; do
  if grep -Eq '^\s*tail:\s*0*\.?0*[1-9]' "$cfg"; then
    echo "$cfg sets report.tail > 0 — the benchmark trims at scoring time, not here." >&2
    exit 1
  fi
done

PREPROC_TOOL=$(command -v porechop_abi > /dev/null && echo porechop_abi || echo "porechop (porechop_abi not installed)")

# --- Stage host genomes (viral libraries) -----------------------------------------
# stage_hosts <config> <host...>: download + minimap2-index any host not yet staged.
stage_hosts() {
  local cfg="$1" host fasta; shift
  for host in "$@"; do
    fasta=$("$RT_PY" - "$cfg" "$host" <<'PY'
import sys
from emtrace.core.config import LazypipeConfig
cfg = LazypipeConfig.from_yaml(sys.argv[1])
print(cfg.resolve_db_path(sys.argv[2]) or "")
PY
)
    if [[ -z "$fasta" ]]; then
      echo "[FATAL] host $host is not defined in $cfg" >&2; exit 1
    fi
    if [[ ! -f "$fasta" ]]; then
      echo "[stage] downloading $host"
      emtrace databases --config "$cfg" --install "$host" -v \
        2>&1 | tee -a "$LOG_DIR/stage_hosts.log" || { echo "[FATAL] install $host failed" >&2; exit 1; }
    fi
    if [[ ! -f "${fasta}.mmi" ]]; then
      echo "[stage] building minimap2 index for $host"
      emtrace databases --config "$cfg" --build-index "$host" -v \
        2>&1 | tee -a "$LOG_DIR/stage_hosts.log" || { echo "[FATAL] index $host failed" >&2; exit 1; }
    fi
    echo "[stage] $host ready: ${fasta}.mmi"
  done
}

# stage_ordered <lib> <fastq.gz> <dir>
# Cut a library into SYNTH_BATCH_READS-read files in FILE ORDER, for deposits that
# carry no start_time (SRA strips it). Headers are written out exactly as deposited:
# no timestamp is invented, so there is no sequencing-time axis and Module C has to
# run in its read-count-only mode. Batch boundaries AND batch order are the
# benchmark's here, which must be stated beside the results.
stage_ordered() {
  local lib="$1" fq="$2" dir="$3"
  [[ -f "$dir/batches.tsv" ]] && { echo "[stage] $lib: $dir ready"; return 0; }
  [[ -s "$fq" ]] || { echo "[FATAL] $fq missing" >&2; return 1; }
  rm -rf "$dir"; mkdir -p "$dir"
  "$RT_PY" - "$fq" "$dir" "$lib" "$SYNTH_BATCH_READS" <<'ORDPY'
import gzip, sys
from pathlib import Path
fq, out, lib, per = sys.argv[1], Path(sys.argv[2]), sys.argv[3], int(sys.argv[4])
n = batch = 0
chunk = []
man = open(out / "batches.tsv.part", "w")
man.write("batch\treads\n")
def flush():
    global chunk, batch
    if not chunk:
        return
    with gzip.open(out / f"{lib}_pass_{batch}.fastq.gz", "wt") as oh:
        oh.writelines(chunk)
    man.write(f"{batch}\t{len(chunk) // 4}\n")
    batch += 1
    chunk = []
with gzip.open(fq, "rt") as fh:
    while True:
        rec = [fh.readline() for _ in range(4)]
        if not rec[0]:
            break
        chunk.extend(rec)
        n += 1
        if n % per == 0:
            flush()
flush()
man.close()
(out / "batches.tsv.part").rename(out / "batches.tsv")
print(f"[stage] {lib}: {n} reads -> {batch} batches of <= {per} (file order, no start_time)")
ORDPY
}

# fetch_tick <pool>: ensure the pool's ENA fastq is on disk; prints its path.
fetch_tick() {
  local pool="$1" srr fq url
  srr=$(tick_srr "$pool")
  [[ -n "$srr" ]] || { echo "[FATAL] unknown tick pool: $pool" >&2; return 1; }
  fq="$TICK_FQ_DIR/$srr.fastq.gz"
  if [[ ! -s "$fq" ]]; then
    mkdir -p "$TICK_FQ_DIR"
    url=$(curl -s --max-time 120 -G "https://www.ebi.ac.uk/ena/portal/api/filereport" \
            --data-urlencode "accession=$srr" --data-urlencode "result=read_run" \
            --data-urlencode "fields=fastq_ftp" --data-urlencode "format=tsv" \
          | awk -F'\t' 'NR == 2 {print $2}' | cut -d';' -f1)
    [[ -n "$url" ]] || { echo "[FATAL] no ENA fastq listed for $srr" >&2; return 1; }
    echo "[stage] tick pool $pool ($srr): downloading" >&2
    curl -sS --max-time 3600 --retry 3 --retry-delay 5 -o "$fq.part" "https://$url" \
      && gzip -t "$fq.part" \
      && mv "$fq.part" "$fq" || { echo "[FATAL] download of $srr failed or is corrupt" >&2; return 1; }
  fi
  echo "$fq"
}

if [[ "$STAGE_HOSTS" == true ]]; then
  [[ "$ORBI_ANY" == true ]] && stage_hosts "$ORBI_BASE_CFG" kitome Aedes_albopictus Homo_sapiens
  [[ "$RUN_TBEV" == true ]] && stage_hosts "$TBEV_CFG" kitome Ixodes_ricinus Ixodes_persulcatus Homo_sapiens
  [[ "$TICK_ANY" == true ]] && stage_hosts "$TICK_CFG" kitome Rhipicephalus_sanguineus Homo_sapiens
fi

# --- Replay one library --------------------------------------------------------
# replay <library> <fastq_dir> <config> <anns> [scoring options...]
replay() {
  local lib="$1" fastq_dir="$2" cfg="$3" anns="$4"; shift 4
  local outdir="$OUT_ROOT/$lib"

  if [[ -f "$outdir/timeseries.tsv" && "$FORCE" != true ]]; then
    echo "[skip] $lib (timeseries.tsv exists; FORCE=true to re-run)"
    return 0
  fi
  # A half-finished replay cannot be resumed: its watch directory already holds
  # released batches. Keep it for inspection and start clean.
  if [[ -d "$outdir" ]]; then
    local aside="${outdir}.incomplete.$(date +%Y%m%d-%H%M%S)"
    echo "[move] previous $lib output -> $aside"
    mv "$outdir" "$aside"
  fi
  mkdir -p "$outdir"

  {
    echo "library:      $lib"
    echo "date:         $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "host:         $(hostname)"
    echo "cpu:          $(sysctl -n machdep.cpu.brand_string 2>/dev/null)"
    echo "ram_gb:       $(( $(sysctl -n hw.memsize) / 1024 / 1024 / 1024 ))"
    echo "threads:      $THREADS"
    echo "fastq_dir:    $fastq_dir"
    echo "config:       $cfg"
    echo "anns:         $anns"
    echo "hit_mode:     all hits (annotation.filter_tophits false)"
    echo "preprocess:   $PREPROC_TOOL"
    echo "emtrace:   $(git -C "$("$RT_PY" -c 'import emtrace,pathlib;print(pathlib.Path(emtrace.__file__).resolve().parents[2])')" log -1 --format='%h %cd' 2>/dev/null)"
    echo "benchmark:    $(git log -1 --format='%h %cd')"
    echo "minimap2:     $(minimap2 --version)"
    echo "scoring:      $*"
  } > "$outdir/provenance.txt"

  echo "[run ] $lib  config=$(basename "$cfg")  anns=$anns"
  local t0=$SECONDS
  # caffeinate keeps the Mac awake for the whole replay (-i idle sleep, -m disk
  # sleep, -s system sleep on AC power); it exits with the wrapped command's status.
  caffeinate -ims ngs-benchmark realtime \
    --fastq-dir "$fastq_dir" \
    --outdir "$outdir" \
    --library "$lib" \
    --sample "$lib" \
    --config "$cfg" \
    --anns "$anns" \
    --threads "$THREADS" \
    --hardlink \
    --batch-timeout "$BATCH_TIMEOUT" \
    "$@" \
    2>&1 | tee "$LOG_DIR/${lib}.log"
  local status=${PIPESTATUS[0]}
  echo "wall_s:       $(( SECONDS - t0 ))" >> "$outdir/provenance.txt"
  echo "exit_status:  $status" >> "$outdir/provenance.txt"

  if [[ $status -ne 0 || ! -f "$outdir/timeseries.tsv" ]]; then
    echo "[WARN] $lib failed (exit $status) — see $LOG_DIR/${lib}.log and $outdir/pipeline/watch.log" >&2
  fi
}

# --- Synthetic MinKNOW batches ----------------------------------------------------
# stage_synthetic <library> <fastq.gz> <dir>: re-cut a concatenated FASTQ into
# MinKNOW-style batch files of SYNTH_BATCH_READS reads, sorted by read start
# time. Start time comes from a MinKNOW `start_time=` header field (header kept
# as is) or a Dorado `st:Z:` tag (header rewritten to `@<id> start_time=<st>`).
# Idempotent: a complete batch dir (manifest written last) is reused; anything
# else is rebuilt.
stage_synthetic() {
  local lib="$1" fq="$2" dir="$3"
  [[ -f "$dir/batches.tsv" ]] && { echo "[stage] $lib: $dir ready"; return 0; }
  [[ -s "$fq" ]] || { echo "[FATAL] $fq missing" >&2; return 1; }
  rm -rf "$dir"; mkdir -p "$dir"
  "$RT_PY" - "$fq" "$dir" "$lib" "$SYNTH_BATCH_READS" <<'STAGEPY'
import gzip, re, sys
from datetime import datetime
from pathlib import Path
fq, out, lib, per = sys.argv[1], Path(sys.argv[2]), sys.argv[3], int(sys.argv[4])
minknow_re = re.compile(r"\bstart_time=(\S+)")
dorado_re = re.compile(r"\bst:Z:(\S+)")
reads, missing = [], 0
with gzip.open(fq, "rt") as fh:
    while True:
        head = fh.readline()
        if not head:
            break
        seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
        if m := minknow_re.search(head):
            st, header = m.group(1), head
        elif m := dorado_re.search(head):
            st = m.group(1)
            header = f"@{head[1:].split()[0]} start_time={st}\n"
        else:
            missing += 1
            continue
        t = datetime.fromisoformat(st.replace("Z", "+00:00"))
        reads.append((t, st, header, seq, qual))
if missing:
    sys.exit(f"{missing} reads in {fq} have no start_time= field or st:Z: tag")
reads.sort(key=lambda r: r[0])   # stable: ties keep file order
with open(out / "batches.tsv.part", "w") as man:
    man.write("batch\treads\tfirst_start_time\tlast_start_time\n")
    for n, i in enumerate(range(0, len(reads), per)):
        chunk = reads[i:i + per]
        with gzip.open(out / f"{lib}_pass_{n}.fastq.gz", "wt") as oh:
            for _, st, header, seq, qual in chunk:
                oh.write(f"{header}{seq}+\n{qual}")
        man.write(f"{n}\t{len(chunk)}\t{chunk[0][1]}\t{chunk[-1][1]}\n")
(out / "batches.tsv.part").rename(out / "batches.tsv")
print(f"[stage] {lib}: {len(reads)} reads -> {n + 1} batches of <= {per}")
STAGEPY
}

# --- Orbivirus (published viral library) --------------------------------------------
# One true taxon at genus (KEVOOV is not in the db); stabilisation on the EM
# probability (±0.01). All cases replay the same batch files.
ORBI_SCORING=(--truth-taxon "$ORBI_TAXON" --rank genus --stabilise-on em_prob --stabilise-tol 0.01)
if [[ "$ORBI_ANY" == true ]]; then
  stage_synthetic "$ORBI_LIB" "$ORBI_FQ" "$ORBI_BATCH_DIR/$ORBI_LIB" || exit 1
fi
[[ "$RUN_ORBI_RefSeq" == true ]] && replay "$ORBI_REFSEQ_RUN" "$ORBI_BATCH_DIR/$ORBI_LIB" \
  "$LOG_DIR/emtrace.config.$ORBI_REFSEQ_RUN.yaml" "$ORBI_ANNS" "${ORBI_SCORING[@]}"
[[ "$RUN_ORBI_RefSeqDmnd" == true ]] && replay "$ORBI_DMND_RUN" "$ORBI_BATCH_DIR/$ORBI_LIB" \
  "$LOG_DIR/emtrace.config.$ORBI_DMND_RUN.yaml" "$ORBI_ANNS" "${ORBI_SCORING[@]}"
[[ "$RUN_ORBI_GBRefSeq" == true ]] && replay "$ORBI_GB_RUN" "$ORBI_BATCH_DIR/$ORBI_LIB" \
  "$LOG_DIR/emtrace.config.$ORBI_GB_RUN.yaml" "$ORBI_ANNS" "${ORBI_SCORING[@]}"

# --- TBEV (published in-house viral library) ---------------------------------------
# One true taxon at species; stabilisation on the EM probability (±0.01).
# The barcodes are concatenated first (gzip members, so plain cat).
TBEV_SCORING=(--truth-taxon "$TBEV_TAXON" --rank species --stabilise-on em_prob --stabilise-tol 0.01)
if [[ "$RUN_TBEV" == true ]]; then
  if [[ ! -s "$TBEV_FQ" ]]; then
    for bc in $TBEV_BARCODES; do
      [[ -s "$TBEV_FQ_DIR/$bc.fastq.gz" ]] || { echo "[FATAL] $TBEV_FQ_DIR/$bc.fastq.gz missing" >&2; exit 1; }
    done
    echo "[stage] $TBEV_LIB: concatenating $TBEV_BARCODES"
    (cd "$TBEV_FQ_DIR" && cat $(printf '%s.fastq.gz ' $TBEV_BARCODES)) > "$TBEV_FQ.part" \
      && mv "$TBEV_FQ.part" "$TBEV_FQ" || { echo "[FATAL] cannot build $TBEV_FQ" >&2; exit 1; }
  fi
  stage_synthetic "$TBEV_LIB" "$TBEV_FQ" "$TBEV_BATCH_DIR/$TBEV_LIB" || exit 1
  if [[ "$TBEV_FIRST_BATCH" == true ]]; then
    first="$TBEV_BATCH_DIR/$TBEV_LIB.first_batch"
    mkdir -p "$first"
    ln -f "$TBEV_BATCH_DIR/$TBEV_LIB/${TBEV_LIB}_pass_0.fastq.gz" "$first/"
    replay "$TBEV_RUN.first_batch" "$first" "$TBEV_CFG" "$TBEV_ANNS" "${TBEV_SCORING[@]}"
  else
    replay "$TBEV_RUN" "$TBEV_BATCH_DIR/$TBEV_LIB" "$TBEV_CFG" "$TBEV_ANNS" "${TBEV_SCORING[@]}"
  fi
fi

# --- Ticks (public field library, PRJNA878809) ------------------------------------
# No truth: each pool is scored against its panel (see header), at GENUS rank: the
# pools carry Aigai-like (genogroup VI) CCHFV, which NCBI files as its own species,
# so the reads split across two species the way the Kevo orbivirus reads split
# across relatives. Stabilisation on the EM probability (+-0.01), as for the other
# viral libraries. Two annotation variants per pool, as for Orbivirus_kevo.
# --reads-only: SRA stripped start_time, so there is no sequencing-time clock.
if [[ "$TICK_ANY" == true ]]; then
  for pool in $TICK_POOLS; do
    [[ -f "$TICK_PANELS/ticks_$pool.panel.tsv" ]] \
      || { echo "[FATAL] no panel for pool $pool: $TICK_PANELS/ticks_$pool.panel.tsv" >&2; exit 1; }
    TICK_SCORING=(--panel "$TICK_PANELS/ticks_$pool.panel.tsv" --rank genus
                  --stabilise-on em_prob --stabilise-tol 0.01 --reads-only)
    tick_fq=$(fetch_tick "$pool") || exit 1
    tick_lib="ticks_$pool"
    stage_ordered "$tick_lib" "$tick_fq" "$TICK_BATCH_DIR/$tick_lib" || exit 1
    [[ "$RUN_TICKS_RefSeq" == true ]] && replay "$tick_lib.refseq_abv1g_k11w5" \
      "$TICK_BATCH_DIR/$tick_lib" "$TICK_CFG" "$TICK_REFSEQ_ANNS" "${TICK_SCORING[@]}"
    [[ "$RUN_TICKS_RefSeqDmnd" == true ]] && replay "$tick_lib.refseq_abv1g_k11w5_dmnd" \
      "$TICK_BATCH_DIR/$tick_lib" "$TICK_CFG" "$TICK_DMND_ANNS" "${TICK_SCORING[@]}"
  done
fi

# --- 16S mock community ----------------------------------------------------------
# Profiling: scored at genus against the 16S-only truth; stabilisation on L1 (±0.02).
# --score-col readn_em: the 10-read cutoff and the profile use EM-expected reads,
# the default of the 16S benchmark scorer (framework_16S.score.py --filter-on).
[[ "$RUN_BARCODE79" == true ]] && replay "barcode79" \
  "$DATA/mock_community.16S/barcode79_fastq_pass" "$S16_CFG" "$S16_ANNS" \
  --truth "$S16_TRUTH" --rank genus --score-col readn_em \
  --stabilise-on L1 --stabilise-tol 0.02

# --- 16S gut mock, R10.4.1 replicates ---------------------------------------------
# Synthetic batches from <library>.fastq.gz (see header and stage_synthetic).
if [[ "$RUN_GUTMOCK" == true ]]; then
  for r in $GUTMOCK_REPS; do
    lib="gutmock_${GUT_MODEL}_r${r}"
    [[ -s "$GUT_FQ_DIR/$lib.fastq.gz" ]] || { echo "[FATAL] $GUT_FQ_DIR/$lib.fastq.gz missing — run polish.16S.comparators.bash --part fetch" >&2; exit 1; }
    stage_synthetic "$lib" "$GUT_FQ_DIR/$lib.fastq.gz" "$GUT_BATCH_DIR/$lib" || exit 1
    replay "$lib" "$GUT_BATCH_DIR/$lib" "$S16_CFG" "$S16_ANNS" \
      --truth "$GUT_TRUTH" --rank genus --score-col readn_em \
      --stabilise-on L1 --stabilise-tol 0.02
  done
fi

echo
echo "Tables: $OUT_ROOT/<library>/{timeseries,first_call,latency,accuracy}.tsv"
