---
name: bp-bulk-rnaseq-quantification
description: Best-practice bulk RNA-seq FASTQ→gene/transcript count matrix — STAR genome alignment + Salmon
  quantification (decoy-aware) via nf-core/rnaseq, producing tximport gene counts ready for DESeq2/edgeR.
when_to_use: Use for the bulk RNA-seq QUANTIFICATION stage only — you have bulk (not single-cell) short-read
  RNA-seq FASTQs from multiple samples and need a raw gene/transcript count matrix ready for differential
  expression. Picks aligner-based (STAR+Salmon/RSEM) vs alignment-free (Salmon/kallisto) and handles strandedness.
  NOT for the DE stats stage (bp-differential-expression), NOT for droplet/plate single-cell FASTQs (bp-raw-data-processing,
  quantify-fastq-to-counts-kb), NOT for small-RNA/miRNA (nf-core/smrnaseq) or splicing/DTU (nf-core/rnasplice).
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core/rnaseq
keywords:
- bulk RNA-seq quantification
- FASTQ to counts
- STAR Salmon
- star_salmon
- alignment-free pseudoalignment
- kallisto
- RSEM
- tximport gene counts
- strandedness auto-detection
- nf-core rnaseq
- count matrix for DESeq2
produces:
- salmon.merged.gene_counts.tsv
- salmon.merged.gene_counts_length_scaled.tsv
- salmon.merged.gene.SummarizedExperiment.rds (SummarizedExperiment)
- transcript counts
- multiqc_report.html
domain: genomics
source: Conesa et al. 2016, A survey of best practices for RNA-seq data analysis (Genome Biology 17:13);
  Salmon (Patro 2017), selective/decoy-aware alignment (Srivastava 2020), tximport (Soneson 2015); nf-core/rnaseq
  docs (nf-co.re/rnaseq).
---

# Bulk RNA-seq quantification: FASTQ → count matrix (best practice)

Question: I have bulk RNA-seq FASTQs from several samples — how do I get a **raw
gene (and transcript) count matrix** I can hand to differential expression?

> This is the **quantification** stage. The DE model/stats (DESeq2/edgeR contrasts,
> replicates, pseudobulk) live in **`bp-differential-expression`** — do not duplicate them here.
> This knowhow feeds that one.

## Recommended: STAR alignment + Salmon quantification (default `star_salmon`)

Map reads to the **genome** with **STAR** (splice-aware), then quantify against the
transcriptome with **Salmon** in alignment-based mode. This is the nf-core/rnaseq
**default** and the field consensus. Why:

- **Salmon corrects the biases that distort raw counting** — fragment GC content,
  positional/sequence bias, and transcript-length differences from isoform switching
  (Patro 2017). Plain `featureCounts`/HTSeq counting does none of this. STAR+Salmon has
  always been nf-core's recommended primary quantification route for exactly this reason —
  Salmon corrects GC/length/positional bias. `featureCounts` was never the default gene
  quantifier; nf-core still runs it, but for **biotype-level QC** (summarising how reads
  distribute across feature types and flagging rRNA), not as the count matrix you hand to DE.
- **Genome alignment gives you QC and reusable BAMs** that pure pseudoalignment cannot:
  rRNA/intronic/intergenic rates, duplication, gene-body coverage, junction saturation
  (RSeQC/Qualimap/MultiQC), plus BAMs for coverage tracks, variant calling, or fusion work.
  Conesa 2016 lists these alignment-QC metrics as core to a trustworthy RNA-seq analysis.
- **Gene-level counts are built correctly** via **tximport** (Soneson 2015): transcript
  estimates are summarised to genes with a length offset, so downstream DESeq2/edgeR see
  proper counts, not TPM. nf-core emits both the raw and **length-scaled** gene count
  matrices plus a `SummarizedExperiment` `.rds`.
- Salmon uses **decoy-aware selective alignment** (the genome as decoy, Srivastava 2020),
  which suppresses spurious assignment of intronic/intergenic reads to transcripts — the
  main accuracy failure of naive transcriptome mapping, especially for total RNA.

**Set strandedness to `auto`.** nf-core sub-samples reads and infers library strandedness
(forward/reverse/unstranded) by pseudoalignment; getting this wrong silently loses most of
your counts, so let it detect and then verify against the RSeQC report.

## Alternatives (and when each is preferred)

- **Alignment-free / pseudoalignment (`--pseudo_aligner salmon` or `kallisto`, no genome
  alignment):** 5–20× faster and far lighter on compute/storage. Prefer when you have a
  **well-annotated genome**, standard **poly-A mRNA**, **many samples** on a tight compute
  budget, and you do **not** need genome BAMs. Accuracy on annotated protein-coding genes is
  comparable to STAR+Salmon; it is weakest for total/rRNA-depleted RNA and poorly annotated
  organisms (alignment-free tools mis-handle reads from unannotated regions). In nf-core you
  can run it **in addition** to `star_salmon` as a fast cross-check.
- **STAR + RSEM (`--aligner star_rsem`):** choose for legacy comparability (e.g. TCGA/GTEx
  pipelines) or when a collaborator expects RSEM TPM/expected-counts. Slower than Salmon,
  similar accuracy.
- **HISAT2 (`--aligner hisat2`):** genome alignment **only — it produces no counts**. Use
  when the goal is a splice-aware BAM for variant calling or novel-transcript work and you'll
  quantify separately; not a route to a count matrix on its own.
- **Different assay entirely → different pipeline:** small RNA / miRNA → **nf-core/smrnaseq**;
  differential transcript usage / splicing → **nf-core/rnasplice**; single-cell/droplet →
  **bp-raw-data-processing** or **quantify-fastq-to-counts-kb** (NOT this pipeline).

## Decision features (what about the data drives the choice)

- **Bulk vs single-cell:** this knowhow is bulk only. Cell/droplet barcodes + UMIs → use the
  single-cell path. (3′ **QuantSeq/Lexogen tag** bulk is still bulk, but see Pitfalls.)
- **Library type:** poly-A mRNA vs total/rRNA-depleted. Total RNA-seq favors **alignment-based**
  STAR+Salmon (better handling of intronic/unannotated reads); poly-A is safe for either.
- **Genome/annotation quality:** well-annotated (human/mouse) → alignment-free is fine;
  non-model or draft annotation → alignment-based (STAR/HISAT2) to genome.
- **Need for genome BAMs:** coverage tracks, variant calling, fusion detection, or novel
  transcript discovery all require STAR/HISAT2 alignment, not pseudoalignment.
- **Sample count & compute:** dozens–hundreds of samples on limited compute → pseudoalignment;
  a handful with full QC needs → STAR+Salmon.
- **Read length:** Salmon's default **k=31** index needs reads **≥75 bp**; short reads
  (≤50 bp, older data) need a smaller k-mer or a genome aligner.
- **Single- vs paired-end:** both supported; declare per-sample in the samplesheet.
- **Strandedness:** unknown → set `auto`; only hard-code if you truly know the kit.
- **UMI-tagged bulk** (e.g. QuantSeq FWD-UMI): enable UMI deduplication.
- **Replicates/depth:** the count matrix must serve DE — aim for **≥3 biological replicates
  per condition**; that requirement is enforced downstream in `bp-differential-expression`.

## Pitfalls

- **Feeding TPM/normalized values into DESeq2/edgeR.** DE needs **raw counts** — use the
  `salmon.merged.gene_counts.tsv` (or the length-scaled matrix / tximport offsets), never
  the TPM/scaled tables.
- **Wrong strandedness** silently zeroes most counts. Use `auto` and confirm in the RSeQC/
  MultiQC strandedness check before trusting the matrix.
- **Skipping the decoy.** A transcriptome-only Salmon/kallisto index inflates counts from
  intronic/intergenic reads; nf-core builds a **decoy-aware** index automatically — keep it.
- **Reference/annotation drift.** The GTF used to build the index must match the one used to
  summarise to genes and the one used downstream — pin a single genome+GTF version.
- **rRNA / DNA contamination** in total RNA-seq — inspect MultiQC rRNA and duplication metrics
  before DE; don't quantify blindly.
- **3′ tag-based assays (QuantSeq):** reads pile at the 3′ end, so gene-length correction is
  wrong — disable Salmon length correction and treat as gene-level 3′ counts.
- **Using this for single-cell FASTQs** — route to the single-cell recipes instead.

## In ABA (which pipeline to run)

Run **nf-core/rnaseq** via **`run_nextflow`** (verified live: latest release **3.26.0**,
nf-co.re/rnaseq). Pin the revision and let strandedness auto-detect:

- **Samplesheet** columns: `sample,fastq_1,fastq_2,strandedness` — set `strandedness: auto`
  per row; leave `fastq_2` empty for single-end.
- **Reference:** `--genome GRCh38` (iGenomes) or explicit `--fasta` + `--gtf`.
- **Aligner:** default `--aligner star_salmon` (recommended). Switch to `star_rsem` for legacy
  comparability, or add `--pseudo_aligner salmon` for a fast alignment-free cross-check. Use
  `--aligner hisat2` only when you want a BAM and will quantify elsewhere.
- **Revision/profile:** pass `-r 3.26.0` (reproducibility) and the container profile for this
  cluster (`singularity`/`apptainer`). Validate the setup first with `-profile test`.
- **UMI bulk:** enable UMI dedup for QuantSeq-UMI libraries.

Call **`describe_pipeline`** for the full, current parameter list (index k-mer, strandedness
thresholds, UMI, trimming, skip flags) — per-param detail lives there, not here.

**Hand-off:** take `salmon.merged.gene_counts.tsv` (raw) / the length-scaled matrix / the
`SummarizedExperiment` `.rds` into **`bp-differential-expression`** (DESeq2/edgeR contrasts),
and the QC into review before DE. For small RNA use **nf-core/smrnaseq**; for splicing/DTU use
**nf-core/rnasplice**.
