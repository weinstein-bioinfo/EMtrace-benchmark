#
# PROJECT:  Benchmarking EMtrace Beta on viral Nanopore data
# DESC:     Framework B — Consensus Assembly Quality on Petersen 2026 simulated reads
#           Run commands line-by-line from the repo root (EMtrace-benchmark/).
# DATE:     2026/05

LOG_DIR="logs"
mkdir -p "$LOG_DIR"


# ============================================================================
# Framework B — Genome Simulation  (Petersen 2026 simulated dataset)
# ============================================================================
#
# Produces:
#   data/petersen_2026/weights/   — positional variability weights (sim_02)
#   data/petersen_2026/mutants/   — 28 artificial genome FASTAs  (sim_03)
#   data/petersen_2026/fastq/     — 28 simulated FASTQ sets       (sim_04)
#   data/petersen_2026/mixed/     — background-embedded FASTQs    (sim_05)
#
# Targets: COVID (SARS-CoV-2), DENV, LASV_L, LASV_S
# Deviations: 0 / 5 / 10 / 15 / 20 / 25 / 30 % mutation from reference
# Dependencies: python3 + biopython + numpy + pyyaml, seqkit, pbsim3
# Status: DONE 2026 May 20
#
# Run the full pipeline in one command (skips sim_01 — using ViMOP alignments):
bash jobs/framework_B/sim_pipeline.sh --targets COVID DENV LASV_L LASV_S --source vimop --mode both \
    --seed 42 --threads 8 --depth 30 2>&1 | tee "$LOG_DIR/sim_pipeline.log"

# --- OR run individual steps ---

# B.1. Compute positional variability weights from ViMOP curated alignments.
#      Reads the pre-aligned multifastas shipped in the ViMOP benchmark;
#      skips the MAFFT alignment step (--source vimop).
#      Output: data/petersen_2026/weights/<VIRUS>.weights.{npy,tsv}
#      Status: DONE 2026 May 20
python3 jobs/framework_B/B_sim02_align_and_weights.py --source vimop --targets COVID DENV LASV_L LASV_S -v \
    2>&1 | tee "$LOG_DIR/B1_sim02_weights.log"

# B.2. Generate 28 artificial mutant genomes (7 deviations × 4 viruses).
#      Applies weighted random substitutions to each reference genome.
#      Output: data/petersen_2026/mutants/<VIRUS>/<VIRUS>_<DEV>.fasta
#      Status: DONE 2026 May 20
python3 jobs/framework_B/B_sim03_generate_mutants.py --targets COVID DENV LASV_L LASV_S --seed 42 \
    2>&1 | tee "$LOG_DIR/B2_sim03_mutants.log"

# B.3. Download 9 background reference genomes from NCBI RefSeq (run once).
#      Resolves GCF accessions via Entrez esearch → esummary → FtpPath_RefSeq,
#      downloads *_genomic.fna.gz, and writes local paths into background_refs.yaml.
#      Must run before B.4 (sim_04 background mode reads the fasta paths written here).
#      See jobs/assets/background_refs.yaml for the 9-organism composition (ViMOP Table S8).
#      Status: DONE 2026 May 20
python3 jobs/framework_B/B_sim03_bg_download.py -v \
    2>&1 | tee "$LOG_DIR/B3_sim_bg_download.log"

# B.4. Simulate nanopore reads with PBSIM3 (viral mutants + background genomes).
#      --mode both: viral reads (7 deviations × 4 viruses at depth 30) and background
#        reads (9 organisms; depth computed per genome so ≥ n_reads × 2 reads produced).
#      Writes fastq paths for background organisms back into background_refs.yaml.
#      Output: data/petersen_2026/fastq/<VIRUS>_<DEV>/<VIRUS>_<DEV>.fastq.gz  (28 files)
#              data/petersen_2026/fastq/<BG_NAME>/<BG_NAME>.fastq.gz            (9 files)
python3 jobs/framework_B/B_sim04_simulate_reads.py --targets COVID DENV LASV_L LASV_S --mode viral --seed 42 --depth 30 --threads 8 \
    2>&1 | tee "$LOG_DIR/B4_sim04_vi_reads.log"
python3 jobs/framework_B/B_sim04_simulate_reads.py --targets COVID DENV LASV_L LASV_S --mode background --seed 42 --depth 30 --threads 8 \
    2>&1 | tee "$LOG_DIR/B4_sim04_bghuman_reads.log"

# B.5. Build background-only pool (500 000 reads) and 28 virus+background mixed FASTQs.
#      Samples exactly n_reads from each background FASTQ (seqkit), then pools.
#      Each mixed FASTQ: 500 virus reads + 500 000 background reads.
#      Output: data/petersen_2026/mixed/background_only_500000.fastq.gz
#              data/petersen_2026/mixed/<VIRUS>_<DEV>_mixed.fastq.gz  (28 files)
python3 jobs/framework_B/B_sim05_mix_background.py --targets COVID DENV LASV_L LASV_S --mode mixed --seed 42 \
    2>&1 | tee "$LOG_DIR/B5_sim05_mix.log"


# clean script-level logs:
rm -fr jobs/logs/sim_*.log
