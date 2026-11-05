#!/usr/bin/env bash
# ============================================================================
# Framework 16S — Spaghetti comparator runner
# ============================================================================
#
# Derived from module1/spaghetti.sh of github.com/adlape95/Spaghetti
# (Latorre-Pérez et al. 2021, Front Microbiol 12:768240; upstream commit
# ac9d2ec). Upstream states no single licence: "Spaghetti is based on many
# open-source tools. Please check the specific licenses of each individual
# tool."
#
# Vendored here rather than patched in place because the clone lives under
# data/, which is gitignored — the exact script that produced a result
# has to be version-controlled.
#
# Changes from upstream, all recorded in the project notes:
#   1. shebang added (upstream opens with a comment block, so exec fails with
#      "Exec format error")
#   2. -t 12 -> "$SPAG_THREADS" (default 8: 12 would starve the OS on a
#      12-core laptop during the long minimap2 steps)
#   3. the authors' hardcoded ~/Descargas/... database paths -> "$SPAG_DB"
#   4. chimera removal made switchable — see SKIP_YACRD below
#   5. ${i%%.*} -> ${i%.fastq} in the chimera block: upstream strips from the
#      FIRST dot, so any dot in the run directory's path (e.g. the Part 2 root
#      results/framework_16S.part2/) sent .scrubb.fastq, .paf and .yacrd
#      outside the run directory. Identical output for dot-free paths, so the
#      Part 1 barcode79 run is unaffected.
#
# Usage (from the repo root):
#   export PATH="$HOME/miniconda3/envs/yacrd/bin:\
#     $PWD/data/comparators/spaghetti-venv/bin:\
#     $PWD/data/comparators/spaghetti-src/module1:$PATH"
#   bash jobs/framework_16S/framework_16S.spaghetti.sh results/framework_16S/spaghetti/
#
# The argument is a DIRECTORY and needs its trailing slash — the script
# concatenates it with a glob. It must contain uncompressed *.fastq (the glob
# is `$1*fastq`), and every intermediate is written beside the input, so give
# the run its own directory rather than pointing at data/.
#
# Requires on PATH: porechop NanoFilt NanoStat minimap2 yacrd
#                   filterPAF.py merfePAF.py taxonomyTable.py
# ============================================================================

# Every step globs the run directory and writes its output back into it, so a
# second run over a dirty directory feeds the pipeline its own intermediates:
# porechop re-chops -porechop.fastq, and NanoFilt/NanoStat get handed the
# -NanoStat.txt report (upstream's ${i::-6} then truncates the name to
# "-NanoSt"). Refuse that rather than produce a plausible mess.
_stale=$(ls "$1"*-porechop* "$1"*-nanofilt* "$1"*scrubb* "$1"*.paf 2>/dev/null | wc -l | tr -d " ")
if [[ "${_stale:-0}" -gt 0 && "${FORCE:-false}" != true ]]; then
  echo "ERROR: $1 already holds $_stale intermediate file(s) from a previous run." >&2
  echo "       Every step globs this directory, so re-running would re-process them." >&2
  echo "       Keep only the input FASTQ, or re-run with FORCE=true." >&2
  exit 1
fi

SPAG_DB="${SPAG_DB:-$PWD/data/comparators/spaghetti_db}"
SPAG_THREADS="${SPAG_THREADS:-8}"
[[ -f "$SPAG_DB/dna-sequences.mmi" && -f "$SPAG_DB/taxonomy.tsv" ]] || {
  echo "Spaghetti DB incomplete: $SPAG_DB (want dna-sequences.mmi + taxonomy.tsv)" >&2
  exit 1; }

#############
# ARGUMENTS #
#############

# 1: input directory (required)

########################
# Read pre-processment #
########################

echo "Starting analysis"
echo "Pre-processing the reads"
echo "------------------------"
echo ""

###################
# Adapter removal #
###################

echo ""
echo "1. Adapter removal"
echo "------------------"

# Only the raw inputs: anything carrying a step suffix is our own output.
time (for i in $(ls $1*fastq 2>/dev/null | grep -Ev -- '-(porechop|nanofilt)|scrubb')
do
	porechop -t "$SPAG_THREADS" -i $i -o ${i::-6}-porechop.fastq
done)

#################################
# Length filtering 1200-1800 pb #
#################################

echo ""
echo "2. Length filtering"
echo "-------------------"

time (for i in $1*-porechop.fastq
do
	cat $i | NanoFilt -l 1200 --maxlength 1800 > ${i::-6}-nanofilt.fastq
done)

######
# QC #
######

echo ""
echo "3. Quality check"
echo "----------------"

time (for i in $1*-porechop-nanofilt.fastq
do
	NanoStat -t "$SPAG_THREADS" --fastq $i > ${i::-6}-NanoStat.txt
done)

###################
# Chimera removal #
###################

echo ""
echo "4. Chimera removal"
echo "------------------"

# SKIP_YACRD=true  : chimera removal is skipped; the nanofilt reads are passed
#                    through unchanged (a symlink, so no copy). Default.
# SKIP_YACRD=false : run the published step, streaming the all-vs-all PAF through
#                    a FIFO so it is never written to disk (~490 GB if it were).
# Rationale and measurements: the project notes
SKIP_YACRD="${SKIP_YACRD:-true}"

time (for i in $1*nanofilt.fastq
do
	if [[ "$SKIP_YACRD" == true ]]; then
		echo "  chimera removal SKIPPED (SKIP_YACRD=true) for $i"
		ln -sf "$(basename $i)" ${i%.fastq}.scrubb.fastq
	else
		mkfifo ${i%.fastq}.paf
		minimap2 -x ava-ont -g 500 -t "$SPAG_THREADS" $i $i > ${i%.fastq}.paf &
		yacrd -i ${i%.fastq}.paf -o ${i%.fastq}.yacrd -c 4 -n 0.4 scrubb -i $i -o ${i%.fastq}.scrubb.fastq
		wait
		rm -f ${i%.fastq}.paf
	fi
done)

###########################
# Read mapping (minimap2) #
###########################

# The script is based on minimap2's page (https://github.com/lh3/minimap2)
# Minimap2 has been used by other authors for 16S rRNA sequencing with ONT
# See: https://f1000research.com/articles/7-1755
# Also: https://elifesciences.org/articles/61504
# Also: https://www.sciencedirect.com/science/article/pii/S2001037019303745

# Silva 138 has been already indexed
# Silva 138 was download in Qiime2 format directly from: https://docs.qiime2.org/2020.8/data-resources/#marker-gene-reference-databases

#################
# Map sequences #
#################

echo ""
echo "5. Mapping"
echo "----------"

# Silva 138
time (for i in $1*scrubb.fastq
do echo "Mapping" $i; minimap2 -x map-ont -t "$SPAG_THREADS" --secondary=no -K 10M "$SPAG_DB/dna-sequences.mmi" $i > $i.paf
done)

# For understanding PAF format: https://github.com/lh3/miniasm/blob/master/PAF.md

#############
# Filtering #
#############

echo ""
echo "6. Filtering the PAF files"
echo "--------------------------"

# Although secondary alignments are turned off, some query reads have been mapped to multiple database seqs.
# Let's filter the output file to just keep one alignment per read.
# We will keep the largest alignment.
# We will remove alignments <500 pb

time (for i in $1*.paf
do filterPAF.py -i $i > ${i::-4}-f.paf
done)

rm -r $1filteredPAFs
mkdir $1filteredPAFs
mv $1*-f.paf $1filteredPAFs

########################################
# Summarize & Merge filtered PAF files #
########################################

echo ""
echo "7. Summarizing the PAF files"
echo "----------------------------"

# It creates a table with the number of sequences assigned to each Database's ID for each sample
# It's something like a QIIME (v1) OTU table.

merfePAF.py -i $1filteredPAFs/ > $1otu_table.csv

####################################
# Create a phyloseq taxonomy table #
####################################

echo ""
echo "8. Creating the taxonomy table"
echo "------------------------------"

# This script creates the taxonomy table needed for loading the data into the phyloseq package.

taxonomyTable.py -i $1otu_table.csv -t "$SPAG_DB/taxonomy.tsv" > $1phyloseq_taxonomy.csv
