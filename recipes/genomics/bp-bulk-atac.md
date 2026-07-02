---
name: bp-bulk-atac
description: Best-practice BULK ATAC-seq — nf-core/atacseq (Trim Galore -> BWA/Bowtie2 -> mito/blacklist/dup
  filter -> Tn5 +4/-5 shift -> MACS2 peaks -> consensus + DESeq2), with ataqv TSS/FRiP QC.
when_to_use: Use for BULK / population / tissue ATAC-seq FASTQs (whole-tissue or sorted bulk
  populations, no barcodes) to get accessible-chromatin peaks, per-sample and consensus peaksets,
  TSS-enrichment/FRiP/fragment-size QC, and DESeq2 differential accessibility across conditions.
  NOT for single-cell / scATAC-seq — see bp-atac.
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core
keywords:
- bulk ATAC-seq
- population ATAC-seq
- tissue ATAC-seq
- bulk chromatin accessibility
- nf-core/atacseq
- Tn5 shift
- nucleosome-free
- MACS2 peak calling
- consensus peaks
- FRiP
- ataqv
- blacklist
- differential accessibility
- DESeq2
produces:
- macs2 narrowPeak/broadPeak
- consensus_peaks.bed
- featureCounts peak x sample matrix
- DESeq2 results
- bigWig tracks
- ataqv/MultiQC QC report
domain: genomics
source: 'Yan, Powell, Curtis & Wu (2020) ''From reads to insight: a hitchhiker''s guide to ATAC-seq data
  analysis'', Genome Biology 21:22 (PMC6996192); nf-core/atacseq docs; ENCODE ATAC-seq standards.'
---

# Bulk ATAC-seq: accessible chromatin (best practice)

Question: from bulk ATAC-seq FASTQs, where is chromatin open, is the library good (TSS enrichment, FRiP, nucleosome laddering), and which regions change accessibility between conditions?

## Recommended approach (and why)
Run the standard **bulk ATAC-seq pipeline**: adapter trim (Trim Galore) -> align (BWA-MEM or Bowtie2, target >80% unique) -> **remove mitochondrial reads, ENCODE blacklist regions, and duplicates** -> apply the **Tn5 shift (+4 bp on +strand, -5 bp on -strand)** to recenter cut sites -> **call peaks with MACS2** -> build a **consensus peakset**, count reads per peak, and test differential accessibility with **DESeq2**. This is exactly the sequence the hitchhiker's-guide review and the ENCODE ATAC standards prescribe, and it is what **nf-core/atacseq** implements end-to-end.

Why this shape:
- **Mito + duplicate + blacklist removal is first-order.** ATAC libraries are heavily contaminated with mitochondrial reads; blacklist/dup removal is what makes replicates reproducible. Skipping it inflates background and wrecks FRiP.
- **Tn5 shift matters** because the transposase inserts as a dimer leaving a 9 bp duplication; +4/-5 recenters reads on the true cut site (critical for peak boundaries and any footprinting). ENCODE-style MACS2 for ATAC uses `--nomodel --shift -37 --extsize 73`.
- **MACS2 is the field default** for count-based peak calling; the review notes there is still no comprehensive ATAC peak-caller benchmark, so MACS2 remains the safe, comparable choice.
- **QC is non-negotiable and specific to ATAC:** fragment-size distribution must show the nucleosome-free (<100 bp) peak plus mono-/di-/tri-nucleosome laddering (~200 bp periodicity); **TSS enrichment** (>6 acceptable, >10 excellent) and **FRiP** gate signal-to-noise. nf-core/atacseq runs **ataqv** + MultiQC to report all of these.

## Alternatives (competing approaches, when preferred)
- **HMMRATAC** — ATAC-specific HMM peak caller that jointly segments nucleosome-free / nucleosomal / background. Prefer when you care about nucleosome positioning or have deep single-sample data; costlier and less standardized than MACS2.
- **Genrich (ATAC mode)** — popular for calling peaks across replicates jointly with built-in mito/multimapper handling; a reasonable MACS2 alternative for replicate pooling.
- **HOMER / SICER-epic2** — broad-domain callers; only if you expect broad accessibility domains rather than punctate peaks.
- **csaw (edgeR, sliding-window)** — for differential accessibility without committing to a fixed peakset; good when peak boundaries shift between conditions. DESeq2-on-consensus-peaks (nf-core default) is simpler and usually sufficient.
- **IDR across replicates** — ENCODE reproducibility metric on peak ranks; use to define high-confidence peaks when you have exactly 2 replicates and need ENCODE-grade peak sets.
- **Genome bins / no peak calling** (as in scATAC snapATAC2) — do NOT use for bulk; peak calling is the right feature definition when signal is not barcode-sparse.

## Decision features (what drives the choice)
- **Single-cell vs bulk** — the primary fork. Per-cell fragments file / cell barcodes -> scATAC (bp-atac). One FASTQ pair per sample/condition -> this recipe.
- **Paired-end vs single-end** — bulk ATAC should be paired-end; PE gives true fragment sizes for nucleosome-free filtering and shift. SE loses the fragment-size QC.
- **Read depth** — ~50M usable reads/sample for open-chromatin + differential analysis; ~200M if you want TF footprinting. Under-sequenced libraries -> weak peaks, unreliable FRiP.
- **Replicates** — >=2 biological replicates per condition are required for DESeq2 differential accessibility; power comes from replicates, not read depth.
- **Peak shape** — punctate regulatory elements (default, narrow/MACS2) vs broad domains (broadPeak / broad callers).
- **Organism / genome** — need a genome with a curated **blacklist** and TSS annotation for QC; iGenomes keys or a custom FASTA+GTF+blacklist.
- **Comparison vs description** — just mapping accessible regions + QC, or contrasting conditions? The latter needs the consensus + DESeq2 stage and a design.

## Pitfalls
- **Not filtering mitochondrial reads / blacklist** — the single most common ATAC failure; destroys FRiP and reproducibility.
- **Forgetting the Tn5 +4/-5 shift** — mis-centers cut sites; matters for boundaries and any footprinting.
- **No fragment-size QC** — a library lacking the nucleosome ladder / with weak TSS enrichment is failing and should not be trusted; always read the ataqv/MultiQC report first.
- **Pseudoreplication in DE** — technical replicates or lanes are not biological replicates; DESeq2 needs true biological replicates.
- **Treating this like scRNA/log-PCA or scATAC/TF-IDF** — wrong modality; bulk ATAC is count-per-peak, analyzed like ChIP-seq/DESeq2, not like a cell x feature matrix.
- **Comparing peaksets across conditions without a consensus set** — always count against one merged consensus peakset before DE.

## In ABA (which pipeline to run)
Route bulk ATAC-seq to **nf-core/atacseq** via `run_nextflow`.
- **Pipeline:** `nf-core/atacseq`. **Revision:** pin the latest stable release **`2.1.2`** (do not run `dev`/`master`).
- **Profile:** use the site container profile (`singularity`/`apptainer` on this cluster; `docker` elsewhere) plus the cluster executor profile; do NOT use `-profile test` for real data.
- **Inputs:** a samplesheet (sample, replicate, fastq_1, fastq_2) and a genome — an iGenomes key or explicit `--fasta`/`--gtf` plus a **blacklist** BED; set narrow vs broad peaks and the DE design as needed.
- **Params:** call `describe_pipeline` for `nf-core/atacseq` to get the exact samplesheet schema, aligner choice (bwa/bowtie2/chromap/star), `--narrow_peak`/`--broad_cutoff`, blacklist/genome params, and DE/consensus options — keep per-param detail there, not here.
- **Selection:** stay on nf-core/atacseq for ATAC. If the assay is actually **ChIP-seq** or **CUT&RUN/CUT&Tag**, switch to `nf-core/chipseq` or `nf-core/cutandrun`. If it's **single-cell** ATAC, this is the wrong recipe -> `bp-atac` (snapATAC2).
