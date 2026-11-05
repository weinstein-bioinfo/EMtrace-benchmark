# Running EPI2ME wf-metagenomics on Roihu — Framework A

`epi2me-labs/wf-metagenomics`: Kraken2 (or minimap2) + Bracken, ONT-native.

## Environment (Roihu)

| what | path / value |
|---|---|
| Nextflow | `module load nextflow` (Roihu-CPU) |
| Container runtime | Apptainer — preinstalled on Roihu compute nodes, no module needed |
| workflow | `epi2me-labs/wf-metagenomics` (pulled automatically by `nextflow run`) |
| raw input FASTQs | `$PROJ/sra/<ACC>.fastq.gz` |
| results | `$PROJ/results/framework_a/epi2me/<dataset>/<ACC>/` |
| Singularity image cache | `$PROJ/singularity_cache` (persistent, shared across array tasks) |

`$PROJ` is `/scratch/project_XXXXXXX/emtrace-benchmark`, same root
every other Framework A script uses.

## Setup

1. **Nextflow** — already in Roihu modules:
   ```bash
   module load nextflow
   ```

2. **Pull the workflow** (optional pre-fetch; `nextflow run` pulls it
   automatically on first use otherwise):
   ```bash
   nextflow pull epi2me-labs/wf-metagenomics
   ```

3. **Container runtime** wf-metagenomics
   ships `standard` (Docker), `singularity`, and `conda` profiles. Roihu
   compute nodes have no Docker daemon, so **`-profile standard` will fail
   there** — use `-profile singularity`. Apptainer (the actively-maintained
   Singularity fork CSC ships) is command-compatible and needs no module
   load; see [docs.csc.fi/apps/nextflow](https://docs.csc.fi/apps/nextflow/).

## Raw reads, not the shared preprocessed reads — and why

`framework_A/metamix/*_preprocess.sbatch` produces a trimmed,
host-filtered FASTQ that both metaMix and MEGAN-LR reuse so the two
EM-analogue comparators see identical input. **EPI2ME deliberately does
NOT reuse that file** — it runs on the raw fetched FASTQ instead. Two
reasons:

1. Per the project notes, EPI2ME's role here is to be the baseline 
   that "demonstrates the FP problem that EM solves."
2. It matches Buddle et al. 2024's own methodology for this exact dataset:
   > "Raw reads were also used as input into either Illumina's Dragen
   > Metagenomics Pipeline or ONT's EPI2ME labs wf-metagenomics as
   > appropriate. The same Kraken2 custom database as used for running
   > Kraken2 through nf-core Taxprofiler was used, and both tools were run
   > using the command-line interface."

   One divergence from Buddle et al.: we run wf-metagenomics' own bundled
   `--database_set Standard-8` rather than their custom Kraken2 DB (not
   available to us).

If you want a run against the host-depleted reads metaMix/MEGAN-LR use instead, 
set `EXCLUDE_HOST=1` in the relevant `framework_A/epi2me/<dataset>.sbatch` 
(uses `--exclude_host` against the same `$hostgenomes/GCF_000001405.40_genomic+rna.fna.gz` 
Homo_sapiens reference `emtrace.config.yaml` points at) — this is a different, non-default run.



## Pipeline

One file per dataset — no RUN_LABEL, no case statement, each file already
has its own `--array` count and accession list baked in:

| dataset | script | produces |
|---|---|---|
| buddle_dna | `framework_A/epi2me/buddle_dna.sbatch` | `.../epi2me/buddle_2024_dna/<ACC>/<ACC>_epi2me.bracken.tsv` |
| buddle_rna | `framework_A/epi2me/buddle_rna.sbatch` | `.../epi2me/buddle_2024_rna/<ACC>/<ACC>_epi2me.bracken.tsv` |
| gauthier | `framework_A/epi2me/gauthier.sbatch` | `.../epi2me/gauthier_2021/<ACC>/<ACC>_epi2me.bracken.tsv` |
| peterson | `framework_A/epi2me/peterson.sbatch` | `.../epi2me/peterson_2026/<ACC>/<ACC>_epi2me.bracken.tsv` |

No DB-build step needed — `--database_set Standard-8` is fetched/managed by
the workflow itself on first use.

```bash
sbatch --array=1 framework_A/epi2me/buddle_dna.sbatch   # smoke-test one task first 
sbatch framework_A/epi2me/buddle_dna.sbatch             # full array
```

Same `assets/*_accs.txt` accession lists as
`framework_A/metamix/<dataset>.sbatch` / `framework_A/megan-lr/<dataset>.sbatch`.

## Output -> ngs-benchmark

The workflow's per-sample Bracken report,
`<out_dir>/bracken/<ACC>.kraken2_bracken.report`, already has the exact
schema `ngs_benchmark.io.adapters.epi2me.load()` expects:

```
name  taxonomy_id  taxonomy_lvl  kraken_assigned_reads  added_reads  new_est_reads  fraction_total_reads
```

The sbatch script copies it to a canonical name,
`<ACC>_epi2me.bracken.tsv`, which is what feeds
`ngs-benchmark --pred ... --adapter epi2me`.






## Reference: `nextflow run epi2me-labs/wf-metagenomics --help`

```
wf-metagenomics v2.14.2-gb060790
--------------------------------------------------------------------------------
Typical pipeline command:

  nextflow run epi2me-labs/wf-metagenomics \ 
	--fastq 'wf-metagenomics-demo/test_data'

Input Options
  --fastq                    [string]  FASTQ files to use in the analysis.
  --bam                      [string]  BAM or unaligned BAM (uBAM) files to use in the analysis.
  --classifier               [choice]  Kraken2 or Minimap2 workflow to be used for classification of reads. [default: kraken2]
                                       * kraken2
                                       * minimap2
  --analyse_unclassified     [boolean] Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified 
                                       directory. 
  --exclude_host             [string]  A FASTA or MMI file of the host reference. Reads that align with this reference will be excluded from the 
                                       analysis. 

Sample Options
  --sample_sheet             [string]  A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory 
                                       containing sub-directories with FASTQ files. 
  --sample                   [string]  A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of 
                                       .fastq(.gz) files. 

Reference Options
  --database_set             [choice]  Sets the reference, databases and taxonomy datasets that will be used for classifying reads. Choices: 
                                       ['ncbi_16s_18s','ncbi_16s_18s_28s_ITS', 'SILVA_138_1', 'Greengenes2_plus', 'Standard-8', 'PlusPF-8', 'PlusPFP-8']. Memory 
                                       requirement will be slightly higher than the size of the database. Standard-8, PlusPF-8 and PlusPFP-8 databases require 
                                       more than 8GB and are only available in the kraken2 approach. [default: Standard-8] 
                                       * ncbi_16s_18s
                                       * ncbi_16s_18s_28s_ITS
                                       * SILVA_138_1
                                       * Greengenes2_plus
                                       * Standard-8
                                       * PlusPF-8
                                       * PlusPFP-8
  --database                 [string]  Not required but can be used to specifically override Kraken2 database [.tar.gz or Directory].
  --taxonomy                 [string]  Not required but can be used to specifically override taxonomy database. Change the default to use a different taxonomy 
                                       file  [.tar.gz or directory]. 
  --reference                [string]  Override the FASTA reference file selected by the database_set parameter. It can be a FASTA format reference sequence 
                                       collection or a minimap2 MMI format index. 
  --ref2taxid                [string]  Not required but can be used to specify a  ref2taxid mapping. Format is .tsv (refname  taxid), no header row.
  --taxonomic_rank           [choice]  Returns results at the taxonomic rank chosen. In the Kraken2 pipeline: set the level that Bracken will estimate abundance 
                                       at. Default: S (species). Other possible options are P (phylum), C (class), O (order), F (family), and G (genus). 
                                       [default: S] 
                                       * S
                                       * G
                                       * F
                                       * O
                                       * C
                                       * P

Kraken2 Options
  --bracken_length           [integer] Set the length value Bracken will use
  --bracken_threshold        [integer] Set the minimum read threshold Bracken will use to consider a taxon [default: 10]
  --kraken2_memory_mapping   [boolean] Avoids loading database into RAM
  --kraken2_confidence       [number]  Kraken2 Confidence score threshold. Default: 0.0. Valid interval: 0-1

Minimap2 Options
  --minimap2filter           [string]  Filter output of minimap2 by taxids inc. child nodes, E.g. "9606,1404"
  --minimap2exclude          [boolean] Invert minimap2filter and exclude the given taxids instead
  --keep_bam                 [boolean] Copy bam files into the output directory.
  --minimap2_by_reference    [boolean] Add a table with the mean sequencing depth per reference, standard deviation and coefficient of variation. It adds a 
                                       scatterplot of the sequencing depth vs. the coverage and a heatmap showing the depth per percentile to the report 
  --min_percent_identity     [number]  Minimum percentage of identity with the matched reference to define a sequence as classified; sequences with a value 
                                       lower than this are defined as unclassified. [default: 90] 
  --min_ref_coverage         [number]  Minimum coverage value to define a sequence as classified; sequences with a coverage value lower than this are defined as 
                                       unclassified. Use this option if you expect reads whose lengths are similar to the references' lengths. 

Antimicrobial Resistance Options
  --amr                      [boolean] Scan reads for antimicrobial resistance or virulence genes
  --amr_db                   [choice]  Database of antimicrobial resistance or virulence genes to use. [default: resfinder]
                                       * resfinder
                                       * ecoli_vf
                                       * plasmidfinder
                                       * card
                                       * argannot
                                       * vfdb
                                       * ncbi
                                       * megares
                                       * ecoh
  --amr_minid                [integer] Threshold of required identity to report a match between a gene in the database and fastq reads. Valid interval: 
                                       0-100 [default: 80] 
  --amr_mincov                [integer] Minimum coverage (breadth-of) threshold required to report a match between a gene in the database and fastq reads. Valid 
                                       interval: 0-100. [default: 80] 

Report Options
  --abundance_threshold      [number]  Remove those taxa whose abundance is equal or lower than the chosen value.
  --n_taxa_barplot           [integer] Number of most abundant taxa to be displayed in the barplot. The rest of taxa will be grouped under the "Other" 
                                       category. [default: 9] 

Output Options
  --out_dir                  [string]  Directory for output of all user-facing files. [default: output]
  --igv                      [boolean] Enable IGV visualisation in the EPI2ME Desktop Application by creating the required files. This will cause the workflow 
                                       to emit the BAM files as well. If using a custom reference, this must be a FASTA file and not a minimap2 MMI format 
                                       index. 
  --include_read_assignments [boolean] A per sample TSV file that indicates the taxonomy assigned to each sequence.
  --output_unclassified      [boolean] Output a FASTQ of the unclassified reads.

Advanced Options
  --min_len                  [integer] Specify read length lower limit.
  --min_read_qual            [number]  Specify read quality lower limit.
  --max_len                  [integer] Specify read length upper limit
  --threads                  [integer] Maximum number of CPU threads to use in each parallel workflow task. [default: 4]

Miscellaneous Options
  --disable_ping             [boolean] Enable to prevent sending a workflow ping.

Other parameters
  --monochrome_logs          [boolean] null
  --validate_params          [boolean] null [default: true]
  --show_hidden_params       [boolean] null

!! Hiding 5 params, use --show_hidden_params to show them !!
--------------------------------------------------------------------------------
If you use epi2me-labs/wf-metagenomics for your analysis please cite:

The nf-core framework
  https://doi.org/10.1038/s41587-020-0439-x
```
