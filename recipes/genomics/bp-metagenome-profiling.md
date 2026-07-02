---
name: bp-metagenome-profiling
description: Profile a shotgun metagenome — read-based taxonomic/functional profiling (Kraken2/Bracken,
  MetaPhlAn via nf-core/taxprofiler) vs de novo assembly + binning into MAGs (nf-core/mag), and when to
  choose which.
when_to_use: Use to decide how to analyze WHOLE-GENOME SHOTGUN metagenomic reads from a microbial community
  (gut/soil/water/host-associated) — pick read-based profiling (who is there / relative abundance / known
  functions) vs assembly+binning to recover MAGs (novel genomes, strain-level, gene clusters). NOT for
  amplicon/16S/ITS/rRNA metabarcoding (that is nf-core/ampliseq), and NOT for single-isolate genome assembly.
keywords:
- shotgun metagenomics
- read-based profiling
- Kraken2
- Bracken
- MetaPhlAn
- taxprofiler
- metagenome assembly
- binning
- MAG
- MetaBAT2
- GTDB-Tk
- host read removal
- taxonomic profiling
- functional profiling
domain: genomics
source: nf-core/mag best-practice pipeline (Krakau et al., NAR Genom Bioinform 2022, PMC8808542) and metagenome
  annotation-approach benchmark (Tamames et al., BMC Genomics 2019, 10.1186/s12864-019-6289-6)
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core
- metagenomics
produces:
- taxonomic_profile_table
- abundance_matrix
- krona_plot
- MAGs_fasta
- bin_quality_report
- gtdbtk_taxonomy
---

# Shotgun metagenome profiling: read-based vs assembly + MAGs (best practice)

Question: given shotgun (whole-genome) sequencing of a microbial community, do I profile reads directly (who/how-much/known-functions) or assemble and bin into metagenome-assembled genomes (MAGs)?

## Recommended approach (with why)

These are **complementary, not mutually exclusive**, but the default entry point is chosen by the scientific question:

- **Read-based profiling is the primary/first-line choice** when the question is *"what taxa are present and at what relative abundance, and how does composition/known-function differ across conditions?"* It is fast, needs no assembly, works at modest depth, and — critically — recovers low-abundance taxa that never assemble. Run **both** a k-mer classifier (Kraken2 + Bracken) and a marker-gene profiler (MetaPhlAn): k-mer tools give high per-read sensitivity but more false positives, while marker-gene tools give higher precision and normalized species-level relative abundance. Benchmarks consistently show these families trade sensitivity for precision, so agreement between them is the robust signal ([Tamames 2019](https://link.springer.com/article/10.1186/s12864-019-6289-6); [sbv IMPROVER, PMC9429340](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9429340/)).

- **Assembly + binning into MAGs is the choice when you need genome-resolved biology**: recovering *novel/uncharacterized* organisms absent from reference databases, full-length genes and biosynthetic gene clusters, strain-level or population genomes, and phylogenomic placement (GTDB-Tk). Confidence in taxonomic assignment is highest for MAGs, then contigs, then reads, and MAGs quantify composition more accurately at species level — but only for organisms abundant/simple enough to assemble ([Integrating MAGs+contigs, PMC11032395](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11032395/)). Deep sequencing is required to assemble genomes from a complex community.

Rule of thumb: **composition / differential abundance / known functional potential → read-based (taxprofiler); discovery of new genomes, strains, or gene content → assembly+binning (mag).** Many studies do read-based first, then assemble deeply-sequenced samples to recover MAGs of interest.

## Alternatives (with caveats — when each is preferred)

- **Amplicon / metabarcoding (16S rRNA, ITS, rpoB)** — a different assay, not a shotgun choice. If the input is PCR-amplified marker reads (not whole-genome), taxonomic resolution is genus/species at best and there is no functional content; use **nf-core/ampliseq**, not this knowhow. Do not run taxprofiler/mag on amplicon data.
- **MetaPhlAn only** — sufficient if you only need well-characterized species-level relative abundance (e.g., human gut) and want low false-positive rates; but it classifies only reads hitting marker genes (often <10% of reads) so it misses under-represented/novel taxa.
- **Kraken2/Bracken only** — fast composition + read-count abundance across a broad database; caveat is false positives, so filter (Bracken threshold, confidence) and confirm with a marker-gene method.
- **Assembly-free functional profiling (HUMAnN, DIAMOND vs gene catalogs)** — for pathway/gene-family abundance without assembly; complements read-based taxonomy. taxprofiler covers DIAMOND-style protein classification.
- **Reference-genome read mapping** — strain-level analysis when a high-quality reference exists; feasible even at shallow depth.

## Decision features (what drives the choice)

- **Goal**: composition / diff-abundance / known functions → read-based. Novel genome, strain, plasmid, or gene-cluster recovery → assembly+MAGs.
- **Sequencing depth**: shallow (~0.5–2M reads, sub-Gbp) is adequate for read-based *composition*; MAG recovery from a complex community needs deep data (commonly ≥5–10 Gbp, tens of millions of reads per sample). Under-sequenced samples yield fragmented, low-quality bins.
- **Community complexity**: low-diversity / high-dominance communities (some clinical, engineered, simple cultures) assemble well at modest depth; highly diverse (soil, sediment) resists assembly → favor read-based, or budget very deep sequencing for MAGs.
- **Reference-database coverage**: well-characterized habitats (human gut/oral) → read-based works because DBs are rich; poorly-characterized environments → assembly to *discover* what references miss.
- **Read length**: long reads (PacBio HiFi, Nanopore) dramatically improve MAG completeness/contiguity and near-complete genomes; both pipelines accept long reads. HiFi is the strongest choice when MAGs are the goal.
- **Host association**: host-associated samples (human/mouse stool, tissue, skin) carry heavy host DNA → **host read removal is mandatory** before both profiling and assembly.
- **Single- vs paired-end / mixed tech**: paired-end short reads are standard; note tech per sample (short vs long) so the pipeline branches correctly.

## Pitfalls

- **Skipping host removal** on host-associated samples → host reads misclassified as microbes, wasted assembly, privacy issues. Always supply a host reference.
- **Attempting MAGs on shallow/complex data** → highly fragmented, contaminated bins; check CheckM/CheckM2 completeness and contamination and GUNC before trusting a bin.
- **Trusting a single classifier** → Kraken2 false positives or MetaPhlAn blind spots; cross-check k-mer vs marker-gene results.
- **Comparing abundances across incompatible outputs** (read counts vs marker-normalized relative abundance) — standardize (taxpasta) before differential analysis.
- **Running shotgun tools on amplicon data** (or vice versa) — wrong assay; use ampliseq for marker genes.
- **Ignoring database choice** — the DB (and its version/date) largely determines what read-based profiling can detect; record it.

## In ABA (which pipeline to run)

- **Read-based taxonomic/functional profiling → nf-core/taxprofiler.** Highly parallel multi-profiler run (Kraken2, Bracken, MetaPhlAn, Centrifuge, Kaiju, mOTUs, DIAMOND, sylph, …) over shotgun short- and long-read FASTQ, with fastp/AdapterRemoval trimming, Bowtie2 (short) / minimap2 (long) host removal, and taxpasta-standardized cross-tool tables + Krona.
- **Assembly + binning into MAGs → nf-core/mag.** MEGAHIT/SPAdes assembly (short + Nanopore/hybrid), binning with MetaBAT2/MaxBin2/CONCOCT/SemiBin2, DAS Tool refinement, CheckM/CheckM2/GUNC/BUSCO QC, GTDB-Tk/CAT taxonomy, geNomad viruses, Prodigal/Prokka annotation.
- **Amplicon data → nf-core/ampliseq** (not covered here; use if reads are 16S/ITS marker amplicons).

Launch with `run_nextflow` on the chosen pipeline (`nf-core/taxprofiler` or `nf-core/mag`). Pin a `revision` to the **latest stable release tag** for reproducibility (do not run `dev`/`master`). Pick an execution `profile` for the compute backend/container engine (e.g. the site Singularity/Apptainer profile on the cluster). **Both pipelines require a database/reference config** (profiler DB CSV for taxprofiler; host-genome + GTDB reference for mag) and a samplesheet declaring per-sample read type — call **`describe_pipeline`** for the exact params (samplesheet schema, `--databases`, host-removal reference, which profilers/binners to enable, depth/QC thresholds) and to resolve the current release tag before launching. Keep read-based and MAG runs as separate invocations; they can share the same trimmed/host-filtered inputs.
