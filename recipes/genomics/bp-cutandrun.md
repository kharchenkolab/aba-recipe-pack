---
name: bp-cutandrun
description: Best-practice CUT&RUN / CUT&Tag occupancy analysis — nf-core/cutandrun (Trim Galore -> Bowtie2
  dual-genome align -> spike-in normalization -> IgG-controlled SEACR/MACS2 peak calling -> consensus
  + fragment/FRiP QC), the low-input in-situ alternative to ChIP-seq.
when_to_use: Use for CUT&RUN, CUT&Tag, or TIPseq FASTQs (antibody + pAG-MNase or pAG-Tn5 in permeabilized
  cells/nuclei, low cell input — hundreds to tens of thousands of cells, e.g. ~50,000 is low-input BULK, NOT
  single-cell — IgG control instead of a sonicated input) to map TF or histone-mark occupancy
  — spike-in-normalized coverage, IgG-gated peaks (SEACR for low-background CUT&Tag, MACS2 as alternative),
  consensus peaksets, and fragment/FRiP QC. NOT crosslink+sonication+input ChIP-seq (use bp-chipseq);
  NOT open-chromatin ATAC-seq (use bp-bulk-atac / bp-atac).
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core
keywords:
- CUT&RUN
- CUT&Tag
- CUTandRUN
- CUTandTAG
- TIPseq
- pAG-MNase
- pAG-Tn5
- nf-core/cutandrun
- spike-in normalization
- E. coli spike-in
- IgG control
- SEACR
- MACS2
- GoPeaks
- low-input chromatin profiling
- low cell number
- limited cell input
- tens of thousands of cells
- 50000 cells
- histone modification
- transcription factor occupancy
- consensus peaks
- FRiP
produces:
- SEACR/MACS2 peaks (bed/narrowPeak)
- consensus_peaks.bed
- spike-in-normalized bigWig tracks
- fragment-length + FRiP + IgG-enrichment QC
- MultiQC report
domain: genomics
source: Meers, Tenenbaum & Henikoff (2019) 'Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin
  profiling', Epigenetics & Chromatin 12:42; Zheng/Henikoff CUT&Tag Data Processing and Analysis Tutorial
  (yezhengstat.github.io/CUTTag_tutorial); Rausch et al. (2025) 'Benchmarking peak calling methods for
  CUT&RUN', Bioinformatics 41(7):btaf375 (PMC12255880); nf-core/cutandrun docs.
---

# CUT&RUN / CUT&Tag occupancy (best practice)

Question: from CUT&RUN or CUT&Tag FASTQs, where does my target (a TF or histone mark) bind, is the library good (spike-in scaling, fragment sizes, FRiP, IgG background), and which regions change between conditions — and when should I have used ChIP-seq instead?

## Recommended approach (and why)
Run the standard **CUT&RUN/CUT&Tag pipeline**: adapter/quality trim (Trim Galore) -> **align with Bowtie2 to BOTH the target genome and a spike-in genome** -> filter/dedup -> **compute a spike-in scale factor from the spike-in read count and apply it to coverage** -> generate normalized bedGraph/bigWig -> **call peaks against the matched IgG control**, using **SEACR** for low-background profiles (the field default for CUT&Tag) or **MACS2** where more background must be modeled -> build a **consensus peakset** and report **fragment-length distribution, FRiP, and IgG-enrichment** QC. This is exactly the flow the Henikoff-lab SEACR paper and CUT&Tag tutorial prescribe, and it is what **nf-core/cutandrun** implements end to end.

Why this shape (and how it differs from ChIP-seq):
- **In-situ tethering, not IP.** An antibody guides pAG-MNase (CUT&RUN) or pAG-Tn5 (CUT&Tag) to cleave/tagment chromatin *in place* in permeabilized cells/nuclei. There is **no crosslinking, sonication, or immunoprecipitation**, so background is far lower and signal-to-noise far higher than ChIP-seq at a tiny fraction of the input (hundreds of cells, down to single cells for CUT&Tag).
- **The control is IgG, not sequenced input.** ChIP-seq subtracts a sonicated *input* track; CUT&RUN/CUT&Tag use a **non-specific IgG** reaction to define background. Peaks are called relative to IgG — the pipeline needs IgG samples flagged as controls.
- **Spike-in normalization is first-order here, not optional.** Because the reaction is so clean, total read depth is a poor normalizer across samples. A defined **spike-in** (E. coli DNA carried with pAG-Tn5 in CUT&Tag, or an exogenous spike-in) gives a per-sample scale factor so tracks and quantitative comparisons are on a common footing. Critically: **once spike-in scaling is applied, tell SEACR NOT to re-normalize** (`non` mode) or you double-normalize.
- **SEACR is purpose-built for low-background data.** Sparse Enrichment Analysis for CUT&RUN thresholds contiguous signal blocks and, given an IgG bedGraph, calls an empirical threshold from the control — well matched to the near-zero background of CUT&Tag. MACS2 (Poisson/local-lambda) is the ChIP-era default and is more appropriate when background is higher or you want ChIP-comparable peak semantics.
- **Assay-specific QC.** Read the **fragment-length distribution** (CUT&Tag of nucleosomal marks shows ~nucleosome-sized and sub-nucleosomal periodicity; TF profiles are shorter), **FRiP**, and **IgG-enrichment ratio** before trusting peaks — the benchmark used IgG-enrichment across all datasets to gate assay quality.

## Alternatives (competing approaches, when preferred)
- **MACS2 instead of SEACR** — prefer for broadly ChIP-comparable calls, when background is elevated (e.g. over-digested CUT&RUN, abundant marks), or to match an existing ChIP-seq analysis. The peak-caller benchmark found MACS2 and SEACR both give balanced precision/recall across marks; SEACR gave the highest SNR.
- **GoPeaks** — newer CUT&Tag-oriented caller with the **highest precision and strongest cross-replicate reproducibility** in the benchmark; prefer when specificity/high-confidence calls matter more than sensitivity, especially for broad repressive marks. (Not in nf-core/cutandrun today — run separately if required.)
- **LanceOtron** — ML caller with maximum sensitivity but lowest precision and weakest reproducibility; use only for discovery with heavy downstream filtering.
- **ChIP-seq (nf-core/chipseq)** — prefer when there is **no CUT&RUN/CUT&Tag-validated antibody**, for **transient/weakly-bound or hard targets that need crosslinking** to capture, when you must match legacy/ENCODE ChIP datasets, or when input material is ample and an established ChIP workflow already exists. This is the main "should I have used ChIP?" fork.
- **Differential binding via DESeq2/edgeR/csaw on the consensus peakset** — for condition contrasts; power comes from **biological replicates**, not depth. nf-core/cutandrun focuses on peaks + QC, so run the DE stage on the consensus count matrix downstream.

## Decision features (what drives the choice)
- **Assay chemistry (the primary fork).** Antibody + pAG-MNase/pAG-Tn5 *in situ*, IgG control, low input -> this recipe. Crosslink + sonicate + IP + sequenced input -> ChIP-seq (bp-chipseq). Tn5 with **no** targeting antibody -> ATAC (bp-bulk-atac/bp-atac).
- **CUT&RUN vs CUT&Tag vs TIPseq** — all handled by the same pipeline; CUT&Tag tolerates the lowest input and highest multiplexing, CUT&RUN is gentler for some antibodies. TIPseq is supported too.
- **Spike-in present?** — is there an E. coli / exogenous spike-in genome to scale by? If yes, enable spike-in normalization and set SEACR to `non`. If no, fall back to read-count/library-size normalization (weaker for cross-condition quantitation).
- **IgG control present?** — matched IgG enables threshold/enrichment-based peak calling and background gating; without it, peak calling is less reliable.
- **Target type / peak shape** — punctate **TF or narrow histone marks** (H3K4me3, H3K27ac) vs **broad domains** (H3K27me3, H3K9me3). Drives narrow vs broad peak mode and caller choice (GoPeaks/SEACR strong on broad).
- **Cell input & depth** — CUT&Tag works from ~hundreds of cells up; a modest count like **tens of thousands (e.g. 50,000) is a LOW-INPUT BULK population and squarely CUT&RUN/CUT&Tag's sweet spot — NOT a reason to reach for a single-cell assay** (single-cell CUT&Tag exists but is a niche, specialized method). Low-complexity libraries need duplicate-aware QC. Typical targets: a few million usable reads for sharp marks/TFs, more for broad domains.
- **Replicates** — >=2 biological replicates per condition for any differential-binding claim (pseudoreplication kills it).
- **Paired-end** — CUT&RUN/CUT&Tag are paired-end; PE fragment sizes are essential for fragment-length QC and for MNase/Tn5 fragment-based peak calling.
- **Organism / genome** — target genome plus a spike-in genome index (commonly E. coli) for the dual alignment.

## Pitfalls
- **Double normalization** — applying spike-in scaling AND leaving SEACR normalization on. Set SEACR to `non` once spike-in scaling is used.
- **Treating IgG like ChIP input, or omitting it** — IgG defines background here; without a matched IgG control, thresholds and enrichment QC degrade.
- **Normalizing by total read depth across samples** — misleading because background is so low; use the spike-in scale factor for cross-sample/condition comparisons.
- **Wrong caller for the mark** — forcing narrow-peak MACS2 settings on broad repressive domains (or vice versa); pick narrow vs broad and consider SEACR/GoPeaks for broad marks.
- **Ignoring fragment-length / FRiP / IgG-enrichment QC** — a flat fragment profile, low FRiP, or poor IgG discrimination means the reaction failed; don't call peaks on a failed library.
- **Using an unvalidated antibody** — CUT&RUN/CUT&Tag are more antibody-sensitive than ChIP; a poor antibody yields clean-looking but wrong occupancy. Validate before trusting.
- **Pseudoreplication in differential binding** — lanes/technical replicates are not biological replicates.
- **Reaching for ChIP-seq tooling with an input track when you have IgG** — modality mismatch; route to the CUT&RUN pipeline, not chipseq.

## In ABA (which pipeline to run)
Route CUT&RUN / CUT&Tag / TIPseq to **nf-core/cutandrun** via `run_nextflow`.
- **Pipeline:** `nf-core/cutandrun`. **Revision:** pin the latest stable release **`3.2.2`** (do not run `dev`/`master`).
- **Profile:** use the site container profile (`singularity`/`apptainer` on this cluster; `docker` elsewhere) plus the cluster executor profile; do NOT use `-profile test` for real data.
- **Inputs:** a samplesheet (group, replicate, fastq_1, fastq_2, **control**) that pairs each target sample with its **IgG control**; the target genome; and a **spike-in genome** (E. coli by default) to enable spike-in normalization. Set the **peak caller** (SEACR default for low-background CUT&Tag; MACS2 for ChIP-comparable/higher-background) and narrow vs broad mode for the target.
- **Params:** call `describe_pipeline` for `nf-core/cutandrun` for the exact samplesheet schema, `--peakcaller` (seacr/macs2), spike-in genome and normalization options, IgG/control handling, `--igg_control`, dedup and consensus settings, and threshold params — keep per-param detail there, not here.
- **Selection:** stay on nf-core/cutandrun for antibody-targeted in-situ profiling. If the assay is actually **crosslink+IP ChIP-seq (sequenced input control)** -> `nf-core/chipseq` (bp-chipseq). If it's **ATAC / open chromatin (no targeting antibody)** -> `bp-bulk-atac` (bulk) or `bp-atac` (single-cell). For **differential binding across conditions**, run peaks here first, then DESeq2/edgeR/csaw on the consensus count matrix.
