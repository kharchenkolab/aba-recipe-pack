---
name: bp-chipseq-peaks
description: Best-practice ChIP-seq occupancy analysis — align, then call TF/punctate binding sites with
  MACS narrow peaks or spread histone domains with MACS broad peaks against an input/IgG control, using
  nf-core/chipseq.
when_to_use: 'Use for crosslinked antibody ChIP-seq (TF or histone-mark) where you want genomic binding
  sites/occupancy: narrow peaks for TFs and punctate marks (H3K4me3/H3K27ac), broad domains for spread
  marks (H3K27me3/H3K36me3/H3K9me3), with an input/IgG control. NOT for ATAC-seq open chromatin (bp-atacseq)
  or CUT&RUN/CUT&Tag targeted-nuclease (bp-cutandrun); NOT for RNA expression.'
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core
keywords:
- ChIP-seq
- peak calling
- MACS2
- MACS3
- narrow peaks
- broad peaks
- transcription factor binding
- histone modification
- input control
- IgG control
- consensus peaks
- IDR
- differential binding
- DiffBind
- phantompeakqualtools
- FRiP
produces:
- narrowPeak/broadPeak files
- consensus peakset
- bigWig coverage tracks
- peak annotation
- MultiQC QC report
- DESeq2 differential-binding results
domain: genomics
source: Landt et al. 2012, ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia, Genome
  Research 22:1813 (PMC3431496); nf-core/chipseq docs (nf-co.re/chipseq)
---

Question: from ChIP-seq reads, where does my TF or histone mark bind the genome, and (optionally) which sites change between conditions?

## Recommended approach (with why)
Run the ENCODE-style pipeline: QC/trim reads, align to the reference, filter duplicates/multimappers/blacklist, then call peaks with MACS against a matched control track, and QC the peaks. This is the community standard implemented by **nf-core/chipseq** (aligner + MACS3 + HOMER annotation + consensus peakset + QC).

The single most consequential decision is **narrow vs broad peak mode**, and it is dictated by the *binding shape* of your target, not by a tuning preference (Landt et al. 2012, ENCODE/modENCODE guidelines):
- **Narrow (point-source) peaks** — for transcription factors and punctate histone marks that localize to focal sites: TFs, H3K4me3, H3K27ac, H3K9ac, and typically H3K4me1. Use MACS narrow mode (`--call-summits`).
- **Broad (spread-source) domains** — for marks that blanket large domains: H3K27me3, H3K36me3, H3K9me3, H3K79me2. Use MACS `--broad`.
- Some marks (e.g. H3K4me1, H3K9me3 near heterochromatin) are "mixed" and are sometimes run both ways.

A **matched control is mandatory**: input (sonicated, non-IP chromatin) is preferred; IgG mock-IP is the alternative. The control models mappability, copy-number, and open-chromatin bias so MACS's local-lambda background does not call artifact peaks. Never call ChIP peaks without a control (this is a core distinction from ATAC-seq, which has no antibody and no input).

**Replicates:** ENCODE requires ≥2 biological replicates. nf-core/chipseq builds a **consensus peakset** (peaks reproducibly present across replicates) and reports FRiP/NSC/RSC. This is the practical reproducibility filter the pipeline ships; the stricter ENCODE **IDR** procedure (rank-consistency of peaks between true replicates) is the gold standard for TF/punctate-mark reproducibility but is NOT implemented in nf-core/chipseq — do it downstream if a submission requires it.

## Alternatives (with when each is preferred)
- **CUT&RUN / CUT&Tag / TIPseq** → use **bp-cutandrun (nf-core/cutandrun)**, NOT this. Targeted-nuclease assays have very low background, need far fewer reads, use IgG (not sonicated input), often use SEACR (or MACS2) for peak calling, and support spike-in normalization. Do not force these through the ChIP-seq pipeline.
- **ATAC-seq** → use **bp-atacseq (nf-core/atacseq)**. No antibody, no input control, Tn5 open-chromatin signal; peak calling and duplicate/Tn5-shift handling differ.
- **Differential binding across conditions** — the built-in consensus→featureCounts→DESeq2 in nf-core/chipseq gives a first-pass quantitative comparison. For a rigorous, publication-grade differential analysis prefer **DiffBind** (MACS peaks → count → DESeq2/edgeR; the most widely used combo for multi-replicate designs) or **csaw** (peak-independent sliding windows + edgeR GLMs, best for complex/ANOVA-like contrasts and diffuse marks). Choose DiffBind for standard two-group TF/mark comparisons; csaw when binding is broad/diffuse or the design has multiple factors.
- **SPP / other callers** — SPP is the ENCODE reference caller for the IDR pipeline; MACS is the default for routine narrow/broad calling and is what nf-core uses.

## Decision features (what drives the choice)
- **Target class** → TF or punctate mark = narrow; spread mark = broad. This is the primary switch.
- **Sequencing depth** → ENCODE target ~≥10M unique reads/replicate for narrow (mammalian), ≥20M for broad marks. Under-sequenced broad marks call poorly.
- **Control present?** → input strongly preferred; IgG acceptable. Absent control = degraded specificity (last resort `--nolambda`, not recommended).
- **Replicates** → ≥2 biological replicates for consensus/reproducibility; more replicates enable differential binding.
- **Read length / single vs paired** → both supported; paired-end improves duplicate detection and fragment-size estimation. Fragment length matters for single-end (cross-correlation/phantompeakqualtools).
- **Organism / genome** → drives reference + blacklist; well-supported for iGenomes references.
- **Goal** → occupancy map only (peaks) vs quantitative comparison between conditions (differential binding) changes whether you stop at consensus peaks or continue to DiffBind/csaw.

## Pitfalls
- Calling narrow peaks on a broad mark (or vice versa) — fragments the domains or misses focal sites. Match mode to target.
- Omitting the input/IgG control — produces false peaks at open/high-copy/artifact regions.
- Not applying a **blacklist** — ENCODE blacklist regions generate spurious high-signal peaks.
- Treating **technical replicates as biological** — inflates reproducibility; consensus/IDR power comes from biological replicates.
- Comparing conditions by overlapping peak lists instead of quantifying reads in a common consensus peakset — do proper count-based differential binding (DESeq2/edgeR via DiffBind/csaw).
- Under-sequencing broad marks and expecting clean domains.
- Assuming nf-core/chipseq did IDR — it did consensus peaks; run IDR separately if required.
- Feeding CUT&RUN/CUT&Tag or ATAC data here — wrong background model and controls.

## In ABA (which pipeline to run)
Route ChIP-seq occupancy/peak-calling to **nf-core/chipseq** via `run_nextflow`.
- Pin a released revision: **`2.1.0`** (latest stable) rather than `dev`/`master`.
- Choose an execution `profile` for the container engine + cluster (e.g. `singularity` on this HPC), composed with the site/institution config.
- Before launching, call **`describe_pipeline`** for nf-core/chipseq to get the exact params (samplesheet columns for sample/replicate/antibody/control, `--narrow_peak` vs default broad, `--macs_gsize`, aligner choice bwa/bowtie2/chromap/star, blacklist, `--read_length`, references). Set narrow-vs-broad and the matched control in the samplesheet per the Decision features above — that is the choice ABA must make for the scientist.
- Sibling routing: **nf-core/cutandrun** for CUT&RUN/CUT&Tag/TIPseq, **nf-core/atacseq** for ATAC. Pick by assay chemistry, not by target.
- For rigorous differential binding beyond the pipeline's DESeq2 pass, follow up with a DiffBind/csaw recipe on the consensus peaks.
