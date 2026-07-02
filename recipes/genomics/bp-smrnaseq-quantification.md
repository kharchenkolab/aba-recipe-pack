---
name: bp-smrnaseq-quantification
description: Best-practice small-RNA-seq / miRNA-seq quantification — protocol-aware adapter+UMI trimming,
  contamination QC, miRBase mapping, isomiR + known/novel miRNA counts via nf-core/smrnaseq.
when_to_use: Use ONLY for dedicated small-RNA / miRNA-seq libraries (size-selected, single-end, ~18-30 nt
  inserts) whose goal is a miRBase-based miRNA / isomiR count matrix. Triggers specifically on miRBase
  mature/hairpin quantification, isomiR annotation (mirtop/mirGFF3), novel-miRNA discovery (miRDeep2), and
  small-RNA-specific profiling (piRNA/tRF). Do NOT use for generic, bulk, poly-A, or total RNA-seq / mRNA
  expression queries (those need an mRNA quant recipe), and NOT for the downstream miRNA differential-expression
  stage (that is bp-differential-expression).
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core
- singularity
keywords:
- small RNA-seq
- miRNA-seq
- miRBase
- mature hairpin
- isomiR
- mirtop
- miRDeep2
- novel miRNA
- miRTrace
- adapter trimming
- UMI
- QIAseq
- NEXTflex
- piRNA
- tRF
produces:
- mirna_counts.tsv
- mirtop.gff3
- mirtrace_qc
- multiqc_report.html
domain: genomics
source: 'Benešová, Kubista & Vališrach 2021, ''Small RNA-Sequencing: Approaches and Considerations
  for miRNA Analysis'' (Diagnostics 11(6):964; PMC8229417); Giraldez et al. 2018 Nat Biotechnol 36:746
  multi-center small RNA-seq benchmark (nbt.4183); nf-core/smrnaseq docs (nf-co.re/smrnaseq).'
---

# Small-RNA-seq / miRNA quantification (best practice)

Question: given small-RNA-seq FASTQs, how do I get accurate miRNA (and isomiR / novel-miRNA / other small-RNA) counts?

## Recommended approach (with why)
Run the **nf-core/smrnaseq** workflow. Small-RNA-seq is NOT a short mRNA-seq run and
must not be quantified like one: the insert (~18-30 nt) is *shorter than the read*, so
every read runs into the 3' adapter and **protocol-specific adapter trimming is mandatory
and load-bearing** — skip or mis-specify it and you quantify adapter dimers, not miRNAs.
Quantification is against the **curated miRBase mature+hairpin reference**, not a
gene/transcript model, because mature miRNAs are 22-nt products of hairpin precursors and
map to many near-identical loci. The recommended flow, all wired in smrnaseq:

1. **Adapter + (optional) UMI trimming** — fastp for the 3' adapter; UMI-tools/UMICollapse
   when the kit carries UMIs (QIAseq). This is the single most consequential step.
2. **Contamination / QC** — **miRTrace** (RNA-type composition, contamination, complexity)
   plus Bowtie2 filtering against rRNA/tRNA/cDNA/ncRNA/piRNA. Small-RNA libraries, especially
   low-input plasma/EV/cfRNA, are dominated by adapter-dimer and rRNA/tRNA fragments.
3. **miRBase mapping + known-miRNA counts** — Bowtie1 alignment to mature + hairpin; per-sample
   count matrix (edgeR-style summarization).
4. **isomiR annotation** — **mirtop** collapses reads and emits the community-standard
   **mirGFF3**, capturing 5'/3' trimming, templated/non-templated additions and editing.
5. **Novel + known miRNA discovery** — optional **miRDeep2** against the genome for candidates
   not in miRBase.

WHY nf-core over a hand-rolled script: library-prep bias dominates small-RNA-seq. The
multi-center benchmark (Giraldez 2018) showed every kit imposes strong protocol- and
sequence-specific ligation bias, so a *reproducible, protocol-parameterised* pipeline with
built-in contamination QC is what makes results interpretable and comparable within a study.

## Alternatives (with caveats)
- **miRge3.0** — fast all-in-one (miRNA + **tRF** + other ncRNA, UMI handling, A-to-I). Prefer
  when you need tRNA-fragment / broad-ncRNA quantification or a lightweight single-tool run.
- **miRDeep2 alone** — the reference for *novel* miRNA prediction, but **no isomiR support**.
  Use when discovery of unannotated miRNAs is the primary goal; smrnaseq already wraps it.
- **sRNAbench / seqcluster / isomiR-SEA / Prost!** — specialist isomiR/clustering tools; all
  can export mirGFF3 via mirtop. Reach for these only for method-specific isomiR studies.
- **Targeted RT-qPCR / miRNA panels (e.g. TaqMan, QIAGEN GeneGlobe)** — preferred when the
  candidate set is small and clinical-grade validation (not discovery) is the goal.
- **Bulk mRNA-seq pipelines (STAR/salmon → gene counts)** — NEVER appropriate here; they
  assume long spliced transcripts, no adapter read-through, and a gene model.

## Decision features (what drives the choice)
- **Library-prep protocol** — the pivotal parameter. Two-adapter ligation (Illumina TruSeq),
  randomized/degenerate adapters (**NEXTflex**), UMI-based (**QIAseq**), poly-A/template-switch
  (**CATS/SMARTer**), circularization (RealSeq). The pipeline's `protocol` (adapter + UMI
  handling) MUST match the kit or counts are garbage.
- **UMIs present?** → enable UMI dedup (QIAseq); otherwise PCR duplicates inflate counts.
- **Read layout** — small-RNA-seq is **single-end, short**; paired-end or long reads signal
  it is NOT a small-RNA library.
- **Organism** — must be in miRBase; set the miRTrace species (e.g. `hsa`, `mmu`).
- **Goal** — known-miRNA counts only (mature/hairpin) vs **isomiRs** (needs mirtop) vs
  **novel miRNAs** (needs a genome + miRDeep2) vs **other small RNAs** (piRNA/tRF → miRge3.0).
- **Input material** — low-input plasma/EV/cfRNA raises adapter-dimer + contamination; lean on
  miRTrace QC and expect low complexity.
- **Replicates / design** — biological replicates per condition are required for any downstream
  differential expression; quantification here just produces the matrix.

## Pitfalls
- Treating it as mRNA-seq (no adapter trim / transcriptome aligner) — the #1 fatal error.
- Wrong `protocol`/adapter → mass loss of reads or quantifying adapter dimers.
- **Cross-protocol comparison** — Giraldez 2018: kits disagree strongly; only compare samples
  prepared with the *same* protocol, and treat absolute abundances as relative-within-study.
- Ignoring contamination — high rRNA/tRNA/dimer fractions (flagged by miRTrace) invalidate counts.
- **miRBase version** drift changes miRNA IDs/counts; pin and record it.
- Multi-mapping across paralogous miRNA loci; forgetting UMI dedup on UMI kits.

## In ABA (which pipeline to run)
Route to **nf-core/smrnaseq** via `run_nextflow`.
- **Pin a revision**: `-r 2.4.1` (current stable) for reproducibility.
- **Profile**: `singularity` on the cluster (Docker not available); combine with an institutional
  config if present.
- **Input**: single-end samplesheet (`sample,fastq_1`). Set the **`protocol`** to match the kit
  (illumina / nextflex / qiaseq / cats), the reference **genome** (for genome-map + miRDeep2),
  and the **miRTrace species**. Toggle miRDeep2 for novel discovery and mirtop for isomiRs.
- Call **`describe_pipeline`** for nf-core/smrnaseq to get the exact param names/defaults
  (protocol, three_prime_adapter, mirtrace_species, mirGeneDB, skip flags) before launching —
  per-param detail lives there, not here.
- **Downstream**: the resulting miRNA count matrix feeds **`bp-differential-expression`**
  (edgeR/DESeq2 on raw counts, needs biological replicates) — quantification alone does not
  test for change between conditions.
