---
name: bp-amplicon-microbiome
description: Best-practice marker-gene amplicon (16S/18S rRNA, ITS) community profiling — primer trimming
  + DADA2 ASV inference + SILVA/UNITE/PR2 taxonomy + diversity, via nf-core/ampliseq.
when_to_use: Use when the input is PCR-amplified marker-gene reads (16S or 18S rRNA, ITS, or COI metabarcoding)
  and the goal is who-is-there community composition, taxonomy, and alpha/beta diversity across samples.
  Requires known forward+reverse PRIMER sequences. NOT for shotgun/whole-metagenome reads or functional/gene/MAG
  profiling — that is bp-metagenomics (nf-core/mag).
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core/ampliseq
- DADA2
- QIIME2
keywords:
- 16S rRNA
- 18S rRNA
- ITS
- amplicon
- metabarcoding
- DADA2
- ASV
- OTU
- primer trimming
- cutadapt
- SILVA
- UNITE
- PR2
- GTDB
- taxonomic assignment
- alpha diversity
- beta diversity
- microbiome community composition
produces:
- ASV_table.tsv
- ASV_seqs.fasta
- taxonomy.tsv
- alpha_diversity.tsv
- beta_diversity/
- barplot/
- qiime2_artifacts/
domain: genomics
source: Callahan, McMurdie & Holmes 2017, ISME J 11:2639 'Exact sequence variants should replace operational
  taxonomic units in marker-gene data analysis' (doi:10.1038/ismej.2017.119); nf-core/ampliseq (DADA2+QIIME2)
---

Question: given PCR-amplified marker-gene reads (16S/18S rRNA, ITS), which microbes are present and how do their abundances and diversity differ across samples?

## Recommended approach (with WHY)
Trim primers, then **infer amplicon sequence variants (ASVs) with DADA2**, assign taxonomy against a marker-appropriate reference, and compute diversity. This is the current field default, implemented end-to-end by **nf-core/ampliseq** (DADA2 + QIIME2).

- **ASVs over OTUs.** Callahan, McMurdie & Holmes (2017, ISME J) argue exact sequence variants "should replace operational taxonomic units": DADA2 models the sequencer error process and resolves biological sequences down to a single nucleotide, giving sensitivity/specificity as good or better than OTU clustering while producing **exact, reusable, dataset-independent sequences** (an ASV means the same thing across studies, unlike a de-novo 97%-OTU whose identity depends on the cluster it landed in). This makes results comparable and re-analysable.
- **DADA2 specifically.** In benchmarks DADA2 gives the best sensitivity; alpha/beta-diversity structure recovered by DADA2 closely matches the intended community. It is the reference denoiser and the pipeline default.
- **Primer trimming is mandatory and non-optional.** Amplicon reads carry the PCR primer at the 5' end; leaving it in corrupts ASV inference and taxonomy. ampliseq removes primers with Cutadapt — you must supply the exact forward/reverse primer sequences (e.g. 515F/806R for 16S V4).
- **Marker-matched reference database.** Taxonomy quality is dominated by using the RIGHT database: **SILVA** for 16S/18S (prokaryotes + eukaryotes), **UNITE** for fungal **ITS**, **PR2** for 18S protists, **GTDB** as a genome-based prokaryotic alternative. Wrong DB = wrong or unassignable taxa.

## Alternatives (with caveats)
- **OTU clustering (VSEARCH/UPARSE, 97%)** — still used for legacy comparability or when over-splitting of ASVs is a concern (strains with intragenomic 16S copy variation can split into several ASVs). ampliseq can optionally post-cluster ASVs into OTUs with VSEARCH. Prefer ASVs by default; add clustering only for cross-study continuity with old OTU tables.
- **Other denoisers (Deblur, UNOISE3)** — comparable to DADA2; DADA2 is the mainstream default and what ampliseq implements. Not a reason to leave the pipeline.
- **Closed-reference / classifier-only tools (e.g. Kraken2, SINTAX)** — faster taxonomy without exact ASVs; useful for quick screens. ampliseq exposes SINTAX/Kraken2/QIIME2 classifiers as alternatives to DADA2's assignTaxonomy.
- **Shotgun / whole-metagenome sequencing** — a DIFFERENT assay, not a DADA2 parameter. Choose shotgun when you need **species/strain resolution, functional potential, or non-bacterial members (viruses, plasmids, genes)**; it samples whole genomes with less amplification bias but costs far more per sample. Amplicon is the workhorse for **large-scale, cost-sensitive who-is-there surveys** of bacteria/archaea (16S) or fungi (ITS). If the reads are NOT primer-amplified marker genes, this is the wrong knowhow — route to **bp-metagenomics (nf-core/mag)**.

## Decision features (what drives the choice)
- **Assay = amplicon vs shotgun.** Primer-amplified single marker region → this knowhow (ampliseq). Random whole-genome fragments → bp-metagenomics/mag.
- **Marker gene → database.** 16S (bacteria/archaea) → SILVA/GTDB; ITS (fungi) → UNITE; 18S (eukaryotes/protists) → SILVA/PR2; COI → dedicated ref. You must know the marker to pick the DB.
- **Primers known?** Forward+reverse primer sequences are REQUIRED inputs. If unknown, they must be recovered before running.
- **Read layout & length / platform.** Paired-end Illumina (must overlap enough to merge — depends on amplicon length vs 2×read length), single-end, or long-read **PacBio/IonTorrent** (full-length 16S; different DADA2 settings). Region choice (V3–V4 vs V4) sets required overlap.
- **Multiple runs / batches.** DADA2 error models are learned per sequencing run; multi-run studies must be denoised per-run then merged — ampliseq handles this when the run is encoded in the samplesheet.
- **Controls.** Negative/blank controls enable decontamination (ampliseq integrates it); mock communities validate ASV accuracy. Low-biomass samples need these.
- **Scale & goal.** Hundreds of samples, budget-limited, community-structure/diversity question → amplicon. Need function or strain ID → shotgun.

## Pitfalls
- **Primers left in reads** → garbage ASVs and misassignment. Always trim; verify primer sequences match the wet-lab protocol.
- **Wrong reference DB for the marker** (e.g. SILVA on ITS) → mostly unassigned or wrong taxa.
- **Merging runs before denoising** → DADA2 error models are per-run; pool ASVs, not raw reads, across runs.
- **Insufficient read overlap** for the amplicon length → paired reads fail to merge and yield is lost; check region vs read length before choosing truncation lengths.
- **Over-aggressive quality truncation** discards large read fractions (DADA2 can retain a minority of reads); tune truncation to the region, don't blindly copy defaults.
- **Treating ASV counts as absolute abundances** — amplicon data are compositional and subject to 16S copy-number and PCR bias; interpret with compositional/diversity methods, and don't compare raw counts across samples without normalization/rarefaction.
- **Over-splitting** of single organisms into multiple ASVs from intragenomic marker copies — consider OTU post-clustering if strain-level splitting confounds the question.

## In ABA (which pipeline to run)
Run **nf-core/ampliseq** via `run_nextflow`. It is the only nf-core pipeline for marker-gene amplicon/metabarcoding (DADA2 ASVs + QIIME2 diversity/taxonomy); the shotgun sibling **nf-core/mag** is for whole-metagenome assembly/binning and belongs to bp-metagenomics.

- **Pin a stable revision** with `-r` (latest stable is 2.18.0 as of 2026-06; prefer the current released tag over `dev`).
- **Profile:** use the platform's institutional/executor profile plus a container profile (`singularity`/`apptainer` on the cluster, per site config); combine with `-profile <site>,singularity`. Use `-profile test` for a smoke run first.
- **Before launching, call `describe_pipeline` for nf-core/ampliseq** to get exact parameter names/defaults. The load-bearing choices to set there: the **primer sequences** (FW/RV), the **marker/database** (SILVA vs UNITE vs PR2/GTDB), read layout (paired/single, or PacBio/IonTorrent flags), per-run/batch grouping in the samplesheet, optional VSEARCH OTU clustering, and decontamination/negative-control settings. Keep those in the pipeline params, not here.
- **Inputs:** a samplesheet of demultiplexed FASTQs (with run/batch column for multi-run studies) plus a metadata sheet for downstream diversity comparisons.
