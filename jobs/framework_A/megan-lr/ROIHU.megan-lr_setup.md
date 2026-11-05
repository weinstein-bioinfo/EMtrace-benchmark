# Running MEGAN-LR on Roihu — Framework A

LAST+MEGAN-LR pipeline (Huson et al. 2018): align long reads/contigs against
a protein reference DB with LAST (MAF output) -> convert to DAA with MAF2DAA ->
meganize the DAA (taxonomic + functional binning) -> export/inspect.

## Environment (Roihu)

| what | path |
|---|---|
| conda env (bioconda: `last` + `megan`) | `/projappl/project_XXXXXXX/venv/megan-lr/bin` |
| megan-map mapping DB (accession -> taxonomy/function) | `/scratch/project_XXXXXXX/emtrace-benchmark/db/megan-map-Jul2020-2-X.db` |
| LAST protein index (built by step 0 below) | `$databases/megan_lr/refseq.abvh.1g.234.*` |
| protein reference FASTA (shared with metaMix/EMtrace) | `$databases/2026/NCBI/refseq.abvh.1g.234/diamond.refseq.abvh.1g.faa` (in the release subdir root, NOT its `diamond/` subdirectory — that holds only the `.dmnd` index) |

## Protein reference DB

Reuse the same cross-pipeline-consistent RefSeq protein FASTA that metaMix's
DIAMOND DB is now built from — `framework_A/db_refseq_abvh_protein.build.sbatch`
(viral + archaea/bacteria 1/genus + human, RefSeq release 234).

## Pipeline

One file per dataset for step 1 — no RUN_LABEL, no case statement, each file
already has its own `--array` count and accession list baked in (`buddle_dna`
| `buddle_rna` | `gauthier` | `peterson`):

| step | script | produces |
|---|---|---|
| 0 (once) | `framework_A/db_refseq_abvh_protein.{dload,build}.sbatch` then `framework_A/megan-lr/builddb.sbatch` | `$databases/megan_lr/refseq.abvh.1g.234.*` (lastdb protein index) |
| 1 (array, per dataset) | `framework_A/megan-lr/<dataset>.sbatch` | `.../megan_lr/<dataset>/<ACC>/<ACC>_megan_taxcount.tsv` |

Requires `framework_A/db_refseq_abvh_protein.build.sbatch` (for the protein
FASTA) and, per dataset, `framework_A/metamix/preprocess/<dataset>_taxprofiler.sbatch`
(for the host-depleted reads — nf-core/taxprofiler's host-removal
subworkflow, reused as-is so both comparators see identical input; same
rationale as `ROIHU.metamix_setup.md` — MEGAN-LR does NOT need
`framework_A/metamix/preprocess/<dataset>_diamond.sbatch` or
`framework_A/metamix/<dataset>.sbatch`, only the taxprofiler step) to
have completed first.


```bash
sbatch framework_A/db_refseq_abvh_protein.dload.sbatch    # if not already run for metaMix
sbatch framework_A/db_refseq_abvh_protein.build.sbatch
sbatch framework_A/megan-lr/builddb.sbatch                # once, builds LAST index
# wait for it, then, per dataset (needs metaMix's taxprofiler step done for each):
sbatch framework_A/metamix/preprocess/buddle_dna_taxprofiler.sbatch  # if not already run for metaMix
sbatch --array=1 framework_A/megan-lr/buddle_dna.sbatch    # smoke-test one task first
sbatch framework_A/megan-lr/buddle_dna.sbatch              # full array
```

## Per-sample steps (inside each `framework_A/megan-lr/<dataset>.sbatch`)

1. `seqkit fq2fa` — trimmed/host-filtered FASTQ -> FASTA (both `lastal` and
   `maf2daa` are given the same FASTA; `maf2daa` needs the
   original query file — FASTA or FASTQ, gzip ok — alongside the MAF, via
   `-r`/`--reads`, "to fix alignment order in the DAA").
2. `lastal -F15 -P <threads> <lastdb-index> reads.fa > reads.maf` — frameshift-
   aware alignment of DNA reads against the protein index (frameshift cost 15,
   the standard LAST+MEGAN-LR long-read setting per the paper's methods.
   Default output format is already MAF.
3. `maf2daa -i reads.maf -r reads.fa -o reads.daa -p <threads>` — MAF -> DAA.
4. `daa-meganizer -i reads.daa -lg -ram readCount -class -mdb $MEGAN_MAPDB -t <threads>` —
   meganize (taxonomic + functional binning), flags per the `-h` dump below.
   `-ram readCount` is required: `-lg` otherwise defaults to `alignedBases`, so
   exported counts would be aligned bases rather than reads.
5. `daa2info -i reads.daa -c2c Taxonomy` — export taxid-to-count. **A
   `taxid\tcount` header line is prepended in the script** before this reaches
   ngs-benchmark: `daa2info`'s own output has no header row by default (its
   `-h` confirms `-l`/`-m`, which add header lines, both default off), but the
   `megan_lr` adapter's loader does `pd.read_csv(path, sep="\t", header=0)`,
   which would otherwise silently treat the first taxon's row as the header.



 


## Reference: `lastdb -h` (captured on Roihu)

```
Prepare sequences for subsequent alignment with lastal.

Main Options:
 -h, --help  show all options and their default settings, and exit
 -p  interpret the sequences as proteins
 -c  soft-mask lowercase letters (in reference *and* query sequences)
 -u  seeding scheme (default: YASS if DNA, else PSEUDO if -q, else exact-match)
 -P  number of parallel threads (default: 1)

Advanced Options (default settings):
 -q  interpret the sequences as proteins and append */STOP
 -S  strand: 0=reverse, 1=forward, 2=both (default: 1)
 -R  lowercase & simple-sequence options (default: 03 for -q, else 01)
 -U  maximum tandem repeat unit length (protein: 50, DNA: 100 or 400)
 -w  use initial matches starting at every w-th position in each sequence (1)
 -W  use "minimum" positions in sliding windows of W consecutive positions (1)
 -Q  input format: fastx, keep, sanger, solexa, illumina (default: fasta)
 -s  volume size (default: unlimited)
 -m  seed patterns (1=match, 0=anything, @=transition)
 -d  DNA seed patterns (N=match, n=anything, R=purine match, etc.)
 -a  user-defined alphabet
 -i  minimum limit on initial matches per query position (default: 0)
 -b  maximum length for buckets
 -B  use max bucket length with memory <= (memory for stored positions) / B (4)
 -C  child table type: 0=none, 1=byte-size, 2=short-size, 3=full (default: 0)
 -x  just count sequences and letters
 -D  print all sequences in lastdb files
 --bits=N  use this many bits per base for DNA sequence (default: 8)
 --circular  these sequences are circular
 -v  be verbose: write messages about what lastdb is doing
 -V, --version  show version information, and exit
```

## Reference: `daa-meganizer -h` (captured on Roihu)

```
Meganizer [options]
DESCRIPTION
	Prepares ('meganizes') a DIAMOND .daa file for use with MEGAN
OPTIONS
 Files  
	-i, --in [string(s)]                 Input DAA file(s). Each is meganized separately. Mandatory option.
	-mdf, --metaDataFile [string(s)]     Files containing metadata to be included in files. 
 Mode  
	-lg, --longReads                     Parse and analyse as long reads. Default value: false.
 Parameters  
	-class, --classify                   Run classification algorithm. Default value: true.
	-ms, --minScore [number]             Min score. Default value: 50.0.
	-me, --maxExpected [number]          Max expected. Default value: 0.01.
	-mpi, --minPercentIdentity [number]   Min percent identity. Default value: 0.0.
	-top, --topPercent [number]          Top percent. Default value: 10.0.
	-supp, --minSupportPercent [number]   Min support as percent of assigned reads (0==off). Default value: 0.01.
	-sup, --minSupport [number]          Min support (0==off). Default value: 0.
	-mrc, --minPercentReadCover [number]   Min percent of read length to be covered by alignments. Default value: 0.0.
	-mrefc, --minPercentReferenceCover [number]   Min percent of reference length to be covered by alignments. Default value: 0.0.
	-mrl, --minReadLength [number]       Minimum read length. Default value: 0.
	-alg, --lcaAlgorithm [string]        Set the LCA algorithm to use for taxonomic assignment. Default value: naive. Legal values: naive, weighted, longReads
	-lcp, --lcaCoveragePercent [number]   Set the percent for the LCA to cover. Default value: 100.0.
	-ram, --readAssignmentMode [string]   Set the read assignment mode. Default value: alignedBases in long read mode, readCount else.
	-cf, --conFile [string]              File of contaminant taxa (one Id or name per line). 
 Classification support:
	-mdb, --mapDB [string]               MEGAN mapping db (file megan-map.db). 
	-on, --only [string(s)]              Use only named classifications (if not set: use all). 
 Deprecated classification support:
	-tn, --parseTaxonNames               Parse taxon names. Default value: true.
	-a2t, --acc2taxa [string]            Accession-to-Taxonomy mapping file. 
	-s2t, --syn2taxa [string]            Synonyms-to-Taxonomy mapping file. 
	-t4t, --tags4taxonomy [string]       Tags for taxonomy id parsing (must set to activate id parsing). 
	-a2ec, --acc2ec [string]             Accession-to-EC mapping file. 
	-s2ec, --syn2ec [string]             Synonyms-to-EC mapping file. 
	-t4ec, --tags4ec [string]            Tags for EC id parsing (must set to activate id parsing). 
	-a2eggnog, --acc2eggnog [string]     Accession-to-EGGNOG mapping file. 
	-s2eggnog, --syn2eggnog [string]     Synonyms-to-EGGNOG mapping file. 
	-t4eggnog, --tags4eggnog [string]    Tags for EGGNOG id parsing (must set to activate id parsing). 
	-a2gtdb, --acc2gtdb [string]         Accession-to-GTDB mapping file. 
	-s2gtdb, --syn2gtdb [string]         Synonyms-to-GTDB mapping file. 
	-t4gtdb, --tags4gtdb [string]        Tags for GTDB id parsing (must set to activate id parsing). 
	-a2interpro2go, --acc2interpro2go [string]   Accession-to-INTERPRO2GO mapping file. 
	-s2interpro2go, --syn2interpro2go [string]   Synonyms-to-INTERPRO2GO mapping file. 
	-t4interpro2go, --tags4interpro2go [string]   Tags for INTERPRO2GO id parsing (must set to activate id parsing). 
	-a2kegg, --acc2kegg [string]         Accession-to-KEGG mapping file. 
	-s2kegg, --syn2kegg [string]         Synonyms-to-KEGG mapping file. 
	-t4kegg, --tags4kegg [string]        Tags for KEGG id parsing (must set to activate id parsing). 
	-a2pgpt, --acc2pgpt [string]         Accession-to-PGPT mapping file. 
	-s2pgpt, --syn2pgpt [string]         Synonyms-to-PGPT mapping file. 
	-t4pgpt, --tags4pgpt [string]        Tags for PGPT id parsing (must set to activate id parsing). 
	-a2seed, --acc2seed [string]         Accession-to-SEED mapping file. 
	-s2seed, --syn2seed [string]         Synonyms-to-SEED mapping file. 
	-t4seed, --tags4seed [string]        Tags for SEED id parsing (must set to activate id parsing). 
	-fwa, --firstWordIsAccession         First word in reference header is accession number (set to 'true' for NCBI-nr downloaded Sep 2016 or later). Default value: true.
	-atags, --accessionTags [string(s)]   List of accession tags. Default value(s): 'gb|' 'ref|'.
 Other:
	-t, --threads [number]               Number of threads. Default value: 8.
	-cs, --cacheSize [number]            Cache size for SQLITE (use with care). Default value: -10000.
	-P, --propertiesFile [string]        Properties file. Default value: ~/.MEGAN.def.
	-v, --verbose                        Echo commandline options and be verbose. Default value: false.
	-h, --help                           Show program usage and quit.
AUTHOR(s)
	Daniel H. Huson.
```

## Reference: `maf2daa` (captured on Roihu)

```
MAF2DAA Converter by Benjamin Albrecht
Copyright (C) 2017 Benjamin Albrecht. This program comes with ABSOLUTELY NO WARRANTY.
Input
	-i, --in                 sets path to MAF-File (can also be piped in, no gzip allowed here)
	-r, -- reads             sets path to query-file in FASTA or FASTQ format (can also be gzipped)
Output
	-o, --out                sets path of the reported DAA-File
Parameter
	-top, --topPercent       sets top percent of reads kept during filtering (default: 10.0)
	-p, --procs              sets number of used processors (default: maximal number)
	-ps, --streamingProcs    sets number of used processors while input is piped-in (default: 1)
	-cs, --chunkSize         sets chunk-size of temporary MAF files (default: 500mb)
	-t, --tmp                sets folder for temporary files (default: parent folder of the resulting DAA-File)
Other
	-v, --verbose            sets verbose mode reporting numbers of reads/references/alignments being analyzed)
	-h, --help               shows program usage and quits
AUTHOR
	Benjamin Albrecht
VERSION
	v0.8.5
```

## Reference: `daa2info -h` (captured on Roihu)

```
SYNOPSIS
	DAA2Info [options]
DESCRIPTION
	Analyses a DIAMOND file
OPTIONS
 Input and Output  
	-i, --in [string]                    Input DAA file. Mandatory option.
	-o, --out [string]                   Output file (stdout or .gz ok). Default value: stdout.
 Commands  
	-l, --list                           List general info about file. Default value: false.
	-m, --listMore                       List more info about file (if meganized). Default value: false.
	-c2c, --class2count [string(s)]      List class to count for named classification(s) (Possible values: EC EGGNOG GTDB INTERPRO2GO KEGG PGPT SEED Taxonomy). 
	-r2c, --read2class [string(s)]       List read to class assignments for named classification(s) (Possible values: EC EGGNOG GTDB INTERPRO2GO KEGG PGPT SEED Taxonomy). 
	-n, --names                          Report class names rather than class Id numbers. Default value: false.
	-p, --paths                          Report class paths rather than class Id numbers. Default value: false.
	-r, --prefixRank                     When reporting class paths for taxonomy, prefix single letter to indicate taxonomic rank. Default value: false.
	-mro, --majorRanksOnly               Only use major taxonomic ranks. Default value: false.
	-bo, --bacteriaOnly                  Only report bacterial reads and counts in taxonomic report. Default value: false.
	-vo, --virusOnly                     Only report viral reads and counts in taxonomic report. Default value: false.
	-u, --ignoreUnassigned               Don't report on reads that are unassigned. Default value: true.
	-s, --sum                            Use summarized rather than assigned counts when listing class to count. Default value: false.
	-es, --extractSummaryFile [string]   Output a MEGAN summary file (contains all classifications, but no reads or alignments). 
AUTHOR(s)
	Daniel H. Huson.
VERSION
	MEGAN Community Edition (version 6.25.10, built 27 Jun 2024).
```

## Reference: `lastal -h` (captured on Roihu)

```
Find and align similar sequences.
 -f  output format: TAB, MAF, BlastTab, BlastTab+ (default: MAF)
 -P  number of parallel threads (default: 1)
 -F  frameshift cost(s) (off)
 -p  match/mismatch score matrix (DNA-DNA: HUMSUM, protein-protein: BLOSUM62, DNA-protein: BLOSUM80)
 -a  gap existence cost (24 if HUMSUM, else 11 if protein reference, else 7)
 -J  score type: 0=ordinary, 1=full (1 for new-style frameshifts, else 0)
 -Q  input format: fastx, keep, sanger, solexa, illumina, prb, pssm (default: fasta)
```
(Full dump has ~40 more lines of E-value/score/split options; all left at
default for the `lastal -F15 -P <threads>` invocation used here — a
DNA-vs-protein frameshift alignment is exactly LAST's designed default case
for `-F`, so `-J`/`-p`/`-a` self-select the right values without extra flags.)
