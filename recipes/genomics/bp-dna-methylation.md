---
name: bp-dna-methylation
description: Best-practice DNA methylation sequencing (WGBS / EM-seq / RRBS) — align + per-cytosine methylation
  calling with nf-core/methylseq (Bismark or bwa-meth+MethylDackel), then differential methylation (DMC/DMR)
  downstream in R.
when_to_use: Use for bisulfite or enzymatic (EM-seq) DNA methylation sequencing reads (WGBS, EM-seq/NEBNext,
  RRBS, PBAT, TAPS) where you need C->T-conversion-aware alignment and per-cytosine (CpG/CHG/CHH) methylation-level
  calls, and optionally differentially methylated cytosines/regions between groups. NOT for 5mC from long-read/Nanopore
  basecalling, methylation ARRAYS (EPIC/450K), MeDIP/MBD enrichment-seq, or RNA/other-modality data.
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core
- bismark
- bwa-meth
- methyldackel
keywords:
- WGBS
- EM-seq
- enzymatic methyl-seq
- RRBS
- bisulfite sequencing
- PBAT
- TAPS
- Bismark
- bwa-meth
- MethylDackel
- per-cytosine methylation
- CpG methylation
- differentially methylated region
- DMR
- DMC
- dmrseq
- DSS
- methylKit
produces:
- bismark_methylation_calls
- cytosine_report.txt.gz
- bedGraph.gz
- coverage_files
- multiqc_report.html
- dmr_results.csv
domain: genomics
source: nf-core/methylseq (Bismark + bwa-meth/MethylDackel); Nunn et al. 2021 'Comprehensive benchmarking
  of software for mapping whole genome bisulfite data' (Briefings in Bioinformatics, bbab021) — peer-reviewed
  aligner benchmark; Kerns & Weber 2025 'Variable performance of widely used bisulfite sequencing methods and
  read mapping software' (bioRxiv PREPRINT, not peer-reviewed; PMC11957057) — single non-model (stickleback)
  example; Korthauer et al. 2018 dmrseq (PMID 29481604).
---

# DNA methylation sequencing (WGBS / EM-seq / RRBS) — best practice

Question: given bisulfite- or enzymatically-converted reads, how do I get per-cytosine methylation levels and find methylation differences between groups — and which conversion protocol / aligner is right for this data?

## Recommended approach (with why)

**Single-base-resolution, conversion-based methylation → align conversion-aware, then call per-cytosine methylation with nf-core/methylseq.** In converted data, unmethylated C reads as T, so a bisulfite/EM-aware aligner (three-letter or wildcard) is mandatory; a normal DNA aligner mismaps massively. The de-facto standard toolchain is **Bismark** (align + methylation extractor → CpG/CHG/CHH cytosine reports, bedGraph, coverage), which nf-core/methylseq runs as its default and validated path.

**Choosing the protocol for a NEW experiment: prefer EM-seq (enzymatic methyl-seq, e.g. NEBNext) over classic bisulfite WGBS.** Enzymatic conversion (TET2 + APOBEC) is far gentler than sodium-bisulfite, so EM-seq gives more even genome coverage, higher/uniform CpG coverage, larger insert sizes, higher library complexity, and works from low or degraded/FFPE input — while producing methylation calls highly concordant with WGBS. Recent multi-platform benchmarks rate EM-seq the best balance of accuracy and practicality for genome-wide profiling. **EM-seq is analyzed on the SAME WGBS code path** (same C→T readout); the only difference is trimming/clip settings for its end-repair chemistry (the pipeline's `--em_seq` preset).

**Differential methylation is a DOWNSTREAM step, not part of the alignment pipeline.** nf-core/methylseq produces per-cytosine calls; DMC/DMR testing between conditions is done afterward in R. For **region-level DMRs with valid FDR control, prefer dmrseq or DSS** (they model biological variance and spatial correlation across CpGs). **methylKit** is the most widely used and is fine for per-site (DMC) tests and window-based DMRs. Run these via a downstream R step over the pipeline's coverage/cytosine-report outputs.

## Alternatives (competing approaches, when each is preferred)

- **bwa-meth + MethylDackel** (`--aligner bwameth`) instead of Bismark: often faster and lighter on memory, with methylation profiles closely comparable to Bismark. In a peer-reviewed WGBS aligner benchmark (Nunn et al. 2021, *Briefings in Bioinformatics*, bbab021), bwa-meth ranked among the fastest and most accurate aligners (high F1, low peak memory) with Bismark close behind — i.e. broadly comparable mapping accuracy, not a large general efficiency gap. A dramatic gap (e.g. bwa-meth ~99% vs Bismark ~54% mapping efficiency) has been reported, but only in a single non-model example (threespine stickleback; Kerns & Weber 2025, bioRxiv **preprint**, not peer-reviewed, PMC11957057), so it should not be treated as a general result. **Consider for large cohorts, deep WGBS, non-model / divergent genomes, or when Bismark's mapping rate is poor on a given dataset.** Bismark remains the safer, most-validated default and integrates cleanly with MultiQC/downstream conventions.
- **RRBS** instead of WGBS/EM-seq: enzyme-digest enrichment of CpG-rich regions (islands, promoters). Interrogates only ~6–12% of CpGs but at high per-site depth for a fraction of the cost. **Prefer for large sample cohorts / population or clinical screens focused on regulatory CpG islands**, where cost per sample dominates and distal/intergenic coverage is not needed. Requires `--rrbs` (MspI diversity/adapter trimming).
- **PBAT** (post-bisulfite adapter tagging): for very-low-input / single-cell-ish WGBS; directional handling differs (`--pbat`). EM-seq now generally outperforms PBAT for low input.
- **TAPS / bwa-mem + rastair** (`--aligner bwamem --taps`): for TET-assisted pyridine borane (non-destructive, no C→T of the bulk genome) libraries — a different chemistry needing its own caller.
- **NOT this recipe:** long-read (Nanopore/PacBio) 5mC comes from the basecaller; methylation **arrays** (EPIC/450K) use minfi/sesame; **MeDIP/MBD** enrichment gives regional, not per-base, signal.

## Decision features (what drives the choice)

- **Assay/conversion chemistry** — EM-seq vs classic bisulfite WGBS vs RRBS vs PBAT vs TAPS. This is the single biggest driver; it sets protocol flag (`--em_seq` / `--rrbs` / `--pbat` / `--taps`) and trimming.
- **Genome scope needed** — every cytosine / distal & intergenic regulatory elements → WGBS or EM-seq; CpG-island/promoter focus only → RRBS is enough and cheaper.
- **Input amount / quality** — low, precious, or degraded/FFPE DNA → EM-seq (or PBAT) over bisulfite WGBS.
- **Cohort size & budget** — many samples, cost-bound → RRBS (deep per-site, cheap); few samples, comprehensive → EM-seq/WGBS.
- **Organism / reference quality** — non-model or divergent genome, or need for max read recovery → bwa-meth; well-annotated model organism → Bismark default is fine.
- **Sequencing depth** — WGBS/EM-seq need deep sequencing for adequate per-CpG coverage; pilot a couple of samples deep and check where mean-methylation estimates plateau before committing the cohort depth.
- **Biological replicates** — required for any credible DMR/DMC test; dmrseq/DSS need replicate groups to estimate variance. Without replicates you can only describe, not test.
- **Single vs paired-end & read length** — paired-end aids alignment and enables SNP-aware filtering (MethylDackel) when no genotypes are available; also matters for RRBS MspI-fragment handling.
- **Cytosine context** — plants/other non-CpG methylation → keep CHG/CHH reporting (`--comprehensive`); mammals → CpG focus usually suffices.
- **Targeted panel** — capture/targeted methyl-seq → `--run_targeted_sequencing` with a regions BED.

## Pitfalls

- **Using a normal (non-conversion) aligner** on converted reads → catastrophic mismapping. Always the bisulfite/EM path.
- **Skipping the M-bias / end-repair trim.** Bisulfite and especially EM-seq read ends carry conversion/end-repair artifacts; not clipping (e.g. EM-seq `--clip_r1/--clip_r2`) biases methylation estimates. Inspect the M-bias plot in MultiQC.
- **No spike-in / conversion-rate check.** Include unmethylated lambda + fully-methylated pUC19 controls; report conversion efficiency — low conversion inflates apparent methylation.
- **Calling DMRs without biological replicates**, or with per-site tests that ignore CpG spatial correlation and over-dispersion → uncontrolled FDR. Prefer dmrseq/DSS; treat single-replicate results as descriptive only.
- **Ignoring SNPs**: C/T SNPs masquerade as unmethylated C. Filter with genotypes, or use paired-end + MethylDackel SNP-aware filtering.
- **Comparing across protocols naively** — RRBS (island-biased) and WGBS/EM-seq cover different CpG populations; don't merge their methylation distributions as if equivalent.
- **PCR-duplicate handling differs by assay** — dedup for WGBS/EM-seq, but do NOT dedup RRBS (fragments legitimately share MspI start sites).

## In ABA (which pipeline to run)

Route to **`nf-core/methylseq`** via `run_nextflow`. It aligns converted reads and produces per-cytosine methylation calls (cytosine reports, bedGraph, coverage) plus MultiQC — it does the alignment→calling stage, NOT differential methylation.

- **Inspect first:** call `describe_pipeline("nf-core/methylseq")` for the exact samplesheet columns (`--input`) and the full parameter set (aligner options, `--em_seq`, `--rrbs`, `--pbat`, `--taps`, clip/trim params, `--comprehensive`, targeted options, resources). Keep per-param detail there — this knowhow only picks the approach.
- **Revision:** pin the latest stable release (currently `4.2.0`) via `revision`; don't run `master` for real data.
- **Aligner selection** (the main analytical fork): default `--aligner bismark` (validated, best-supported); switch to `--aligner bwameth` (bwa-meth + MethylDackel) for large cohorts, deep WGBS, non-model genomes, or when mapping rate is poor; `--aligner bwamem --taps` for TAPS.
- **Protocol flags:** set `--em_seq` for EM-seq, `--rrbs` for RRBS (and do not dedup), `--pbat` for PBAT — these are the correctness-critical switches; get them from the user's library prep.
- **Profile / execution:** use `-profile test` for a smoke run; add your site profile (e.g. `test,cbe`) for real runs. Keep `background=True` (default) for real WGBS/EM-seq — head runs as a long Slurm job and fans tasks out; use `execution="local"` only for tiny `-profile test` runs.
- **Reference:** provide the genome (iGenomes key or FASTA); the run builds/uses the bisulfite index.
- **Differential methylation (downstream):** after the run, feed the coverage/cytosine-report files to an R step (`run_r`/`run_python`) using **dmrseq or DSS** for replicate-based DMRs with FDR control, or **methylKit** for per-site/window tests. This is a separate stage, not a methylseq parameter.
