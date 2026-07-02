---
name: bp-viral-genome-reconstruction
description: Reference-based viral consensus genome reconstruction + intra-host/low-frequency variant
  calling from amplicon (ARTIC tiling) or metagenomic/capture short/long reads, with lineage assignment
  (pangolin/nextclade) — via nf-core/viralrecon.
when_to_use: Use ONLY for VIRAL consensus-genome reconstruction and reference-based variant calling against a
  KNOWN reference virus (e.g. SARS-CoV-2/ARTIC, influenza, RSV, mpox) from amplicon (ARTIC/midnight tiling),
  metagenomic-shotgun, or hybrid-capture reads (Illumina or Nanopore) — i.e. viral genomic surveillance,
  lineage/clade typing, and intra-host viral variant calling. This is virus-specific — do NOT match generic
  "genome" or "variant calling" requests. NOT for host/human or other eukaryotic variant calling, NOT for
  bacterial or eukaryotic genomes, NOT for de novo assembly of a novel/unknown virus (use metagenomic
  assembly), NOT for host RNA/DNA variant calling.
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core
- apptainer
keywords:
- viral consensus genome
- ARTIC amplicon
- primer scheme
- iVar
- bcftools consensus
- intra-host variants
- pangolin lineage
- nextclade
- SARS-CoV-2 surveillance
- viralrecon
- primer trimming
- Freyja wastewater
produces:
- consensus.fa
- variants.vcf
- pangolin_lineage.csv
- nextclade.tsv
- multiqc_report.html
- coverage_plots
domain: genomics
source: nf-core/viralrecon v3.0.0 (https://github.com/nf-core/viralrecon/releases/tag/3.0.0;
  nf-co.re/viralrecon; Patel et al.) + ARTIC Network protocol (artic.network)
  and benchmark literature (PMC12394145; CoVpipe2 lessons PMC11074694)
---

# Viral consensus genome reconstruction & variant calling (best practice)

Question: given reads from a sample of a KNOWN reference virus, what is this sample's
consensus genome, what variants does it carry, and what lineage/clade is it?

The field standard for this (surveillance-scale, known reference) is **reference-based**
mapping → primer-aware trimming → variant calling → masked consensus → lineage assignment,
**not** de novo assembly. **nf-core/viralrecon** implements exactly this and is the community
workhorse for SARS-CoV-2 and other segmented/small viral genomes. The single most important
decision is **library type**, because it dictates the whole downstream tool chain:

- **Amplicon (tiling PCR — ARTIC / midnight / custom):** the dominant surveillance method.
  At low viral load it outperforms shotgun/metagenomic (which fails when the genome is buried in
  host reads), so it is the go-to for low-load clinical samples — **BUT genome completeness still
  degrades at low titre / high Ct**: all amplicon protocols drop off as titre falls (PMC12394145
  found this for every amplicon scheme it tested, which motivated a new ARTIC-Amp method precisely
  to recover completeness at low titre). **You MUST trim primer sequences** (position-based, via the
  scheme BED) before variant calling, or primer-derived bases masquerade as real variants.
  Default caller **iVar** (amplicon-aware). Caveat: PCR/primer bias gives uneven coverage and
  **amplicon dropout** (primer-site mutations or primer-primer dimers → whole amplicons drop →
  runs of `N`), and low-frequency intra-host allele frequencies are distorted by amplification.

- **Metagenomic / shotgun (no primers):** even coverage, faithful allele frequencies (best for
  intra-host diversity), and captures the whole genome without scheme design — but needs **high
  viral load / low host background** or the genome is buried in host reads. Default caller
  **bcftools**. Preferred when titre is high, when studying within-host variation, or when no
  validated primer scheme exists for the target.

- **Hybrid-capture (probe enrichment):** middle ground — tolerates more host background than
  shotgun and is robust to lineage-defining primer-site mutations that break amplicon schemes,
  at higher cost. Good for surveillance across divergent variants. Handled by viralrecon as the
  metagenomic (primer-free) path.

## Alternatives (with caveats)
- **De novo assembly (SPAdes/metaSPAdes, or viralrecon's assembly sub-workflow):** use when the
  target is **novel/highly divergent** or you have **no suitable reference**. viralrecon can run
  de novo assembly on the Illumina path, but for a known reference the mapping/consensus path is
  more accurate and interpretable — reserve assembly for reference-free situations.
- **ARTIC fieldbioinformatics CLI / `artic minion`:** the original Nanopore-amplicon consensus
  tool (medaka/clair3). viralrecon's Nanopore path wraps this; use the standalone CLI only for
  bespoke single-sample runs. For batch/surveillance, viralrecon is preferable (parallel,
  MultiQC, pangolin/nextclade built in).
- **CoVpipe2 / other SARS-CoV-2 pipelines:** equivalent reference-based approach; pick viralrecon
  in ABA because it is the maintained nf-core implementation with the broadest platform/protocol
  matrix.
- **Freyja (relative lineage abundance):** for **mixed samples / wastewater** where you want
  lineage *proportions*, not a single consensus — available inside viralrecon.

## Decision features (what drives the choice)
- **Library prep / assay:** amplicon (tiling) vs metagenomic-shotgun vs capture → sets
  `--protocol` and the whole caller/trim chain. This is the primary switch.
- **Sequencing platform:** Illumina vs Nanopore → `--platform`. On viralrecon, **Nanopore is
  amplicon-only** (ARTIC medaka/clair3 workflow); Illumina supports amplicon, metagenomic, AND
  capture. If reads are Nanopore metagenomic, viralrecon is not the tool.
- **Viral load / Ct value:** low load (high Ct) → amplicon preferred (metagenomic fails outright),
  though even amplicon completeness degrades as titre falls; high load → metagenomic viable and
  gives cleaner intra-host frequencies.
- **Primer scheme identity + version:** the wet-lab scheme (ARTIC V3/V4/V4.1/V5.3.2, midnight
  1200bp, custom) MUST match `--primer_set`/`--primer_set_version` or a BED. Wrong version →
  mis-trimming, spurious variants, and dropout misattribution.
- **Read length / amplicon size:** 400 bp ARTIC vs 1200 bp midnight/long-amplicon changes
  coverage-vs-dropout tradeoff (long amplicons rescue primer-site-mutation dropout).
- **Organism / reference:** SARS-CoV-2 gets pangolin + nextclade + Freyja automatically; other
  viruses reconstruct a consensus but need an appropriate `--genome`/reference (and may have no
  lineage caller — nextclade datasets exist for some, e.g. flu/RSV/mpox).
- **Analysis goal:** single consensus + lineage (surveillance) vs intra-host low-frequency
  variants (needs metagenomic + adequate depth) vs lineage *abundance* in a mixture (Freyja).
- **Host contamination:** clinical/shotgun samples → enable host read removal (Kraken2) so host
  reads don't distort mapping/coverage.

## Pitfalls
- **Skipping / mis-versioning primer trimming (amplicon):** the #1 error — untrimmed or
  wrong-scheme primers create false variants at primer sites. Always confirm the scheme+version
  with the person who did the prep.
- **Amplicon dropout read as deletion:** runs of `N` are usually failed amplicons (primer-site
  mutation, low input), not real deletions — check the per-amplicon coverage/MultiQC before
  interpreting.
- **Trusting low-frequency (iSNV) calls from amplicon data:** PCR amplification distorts allele
  frequencies; for genuine intra-host diversity use metagenomic/capture data.
- **Metagenomic on a low-titre sample:** genome buried in host reads → incomplete consensus;
  amplicon would have recovered more (though even amplicon completeness degrades at very low titre).
- **Consensus N-masking threshold:** low-coverage positions must be masked to `N` (min-depth,
  typically ~10x); an unmasked consensus over-calls bases in poorly covered regions.
- **Stale lineage databases:** pangolin/nextclade calls are only as current as their datasets —
  pin/record the DB version, and note lineage drift over time.
- **Wrong reference/`--genome`:** reference-based reconstruction against the wrong reference or
  variant silently produces garbage for divergent samples — reserve for known-reference cases.

## In ABA (which pipeline to run)
Route to **nf-core/viralrecon** via `run_nextflow`. Pin a released revision — latest is
**`3.0.0`** (Oct 2025); `2.6.0` is the widely-cited prior stable. Use the `singularity`
(Apptainer) profile on this cluster (add the site/institution profile as configured).

Before launching, **call `describe_pipeline("viralrecon", revision=...)`** to pull the exact,
version-correct params — do NOT hardcode param details here. The selection params you steer with
(confirm each against the sample metadata / prep):
- `--platform` (illumina | nanopore),
- `--protocol` (amplicon | metagenomic),
- `--primer_set` + `--primer_set_version` (or a primer `--primer_bed`) for amplicon,
- `--variant_caller` (ivar for amplicon / bcftools for metagenomic — the defaults follow
  `--protocol`) and `--consensus_caller` (bcftools|bedtools default, or ivar),
- `--genome` / reference for the target virus,
- host-removal (Kraken2) and skip-toggles for QC/lineage steps.

Provide a standard nf-core `--input` samplesheet (sample, fastq_1, fastq_2). Outputs to hand
downstream: per-sample `consensus.fa`, `variants.vcf`, pangolin/nextclade lineage tables, and the
MultiQC report (inspect per-amplicon coverage and %N before trusting the consensus). If the target
is a **novel/reference-free** virus, this recipe does not apply — use a metagenomic de novo
assembly route instead.
