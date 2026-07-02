---
name: bp-scrnaseq-quantification
description: Best-practice droplet scRNA-seq FASTQ→cell×gene count matrix — barcode/UMI-aware mapping
  (alevin-fry/simpleaf default, or STARsolo/Cell Ranger/kallisto|bustools) with empty-droplet calling,
  via nf-core/scrnaseq.
when_to_use: Use for the droplet single-cell QUANTIFICATION stage only — you have 10x/Chromium (or other
  droplet) scRNA-seq FASTQs with cell barcodes + UMIs and need a raw cell×gene count matrix (h5ad/mtx)
  to feed cell QC and the single-cell chain. Picks among alevin-fry/simpleaf, STARsolo, Cell Ranger, and
  kallisto|bustools, and routes feature-barcode/multiplexed (CITE-seq, hashing, cellranger multi) and
  multiome (GEX+ATAC) variants. NOT the conceptual FASTQ→counts walkthrough (bp-raw-data-processing) or
  the Python kb-python path (quantify-fastq-to-counts-kb); NOT bulk RNA-seq (bp-bulk-rnaseq-quantification);
  NOT cell QC/filtering (bp-quality-control) or downstream clustering/DE.
keywords:
- droplet scRNA-seq quantification
- 10x Chromium FASTQ to counts
- cell barcode UMI
- alevin-fry simpleaf
- STARsolo
- Cell Ranger
- kallisto bustools
- empty droplet calling
- emptyDrops
- single-nucleus GeneFull
- cellranger multi feature barcode
- nf-core scrnaseq
domain: genomics
source: Single-cell Best Practices (Heumos et al. 2023, Nat Rev Genet) — sc-best-practices.org/introduction/raw_data_processing.html;
  Brüning et al. 2022 'Comparative analysis of common alignment tools for single-cell RNA sequencing'
  (GigaScience 11:giac001); You et al. 2021 'Benchmarking UMI-based single-cell RNA-seq preprocessing
  workflows' (Genome Biology 22:339); nf-core/scrnaseq docs (nf-co.re/scrnaseq).
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core/scrnaseq
produces:
- count_matrix.h5ad
- filtered_feature_bc_matrix
- raw_feature_bc_matrix
- multiqc_report.html
- alevinqc report
---

# Droplet scRNA-seq quantification: FASTQ → cell×gene count matrix (best practice)

Question: I have droplet (10x/Chromium) scRNA-seq FASTQs with cell barcodes + UMIs —
how do I get a raw **cell × gene count matrix** to hand to cell QC and the single-cell chain?

> This is the **upstream preprocessing / quantification** stage: map reads → correct cell
> barcodes → resolve UMIs → quantify → call cells vs empty droplets. Cell QC/filtering lives in
> **`bp-quality-control`**; normalization, clustering, annotation, DE are the downstream bp-* chain.
> For the CONCEPTUAL flow (what each step means) see **`bp-raw-data-processing`**; for a hands-on
> Python kb-python run see **`quantify-fastq-to-counts-kb`**. This knowhow is the **pipelined,
> reproducible** route and picks the mapper for the agent.

## Recommended: alevin-fry / simpleaf (nf-core/scrnaseq default)

Map with **simpleaf** (piscem selective/pseudo-alignment → **alevin-fry** for barcode
correction, UMI resolution, and quantification). This is the nf-core/scrnaseq **default** and a
sound field choice for a scripted, cluster-friendly run. Why:

- **Concordant with the reference standard.** Benchmarks show alevin-fry (and STARsolo)
  produce count matrices and gene sets closely matching **Cell Ranger**, the de-facto 10x
  standard (You et al. 2021; Brüning et al. 2022). You do not sacrifice accuracy to gain speed.
- **Fast and memory-frugal.** Lightweight mapping to a transcriptome (or augmented
  transcriptome) is several-fold faster than full genome alignment and runs in bounded memory,
  which matters at scale on shared compute (sc-best-practices, Heumos et al.).
- **Open-source and reproducible.** No commercial licence, pinned containers, deterministic
  across runs — the right default for an automated pipeline. Cell Ranger's nf-core container is
  **community-built and not supported by 10x** (see Pitfalls), so it is a deliberate opt-in.
- **Correct barcode/UMI handling built in.** 10x whitelist barcode correction (Hamming ≤1),
  UMI deduplication with EM resolution of multimapping reads, and an unfiltered + knee-based cell
  call — exactly the four steps the best-practice book prescribes — with an **AlevinQC** report.

## Alternatives (and when each is preferred)

- **STARsolo (`--aligner star`):** spliced alignment to the **full genome**. Prefer when you
  need **genome BAMs** (coverage, variant calling, novel-transcript work), a **non-model / poorly
  annotated** organism, or **single-nucleus** data where intronic reads matter (`--star_feature
  GeneFull` counts pre-mRNA). Benchmarks rate STARsolo as the closest substitute for Cell Ranger
  — nearly identical cells/genes, ~5× faster (Brüning et al. 2022).
- **Cell Ranger (`--aligner cellranger`):** the 10x commercial standard. Choose when a
  collaborator/journal/regulator **requires exact Cell Ranger output**, or for downstream tools
  that expect its precise matrix. In nf-core the container is community-built (not 10x-supported).
- **Cell Ranger multi (`--aligner cellrangermulti`):** the route for **Feature Barcode**
  experiments — **CITE-seq** (antibody-derived tags), **cell hashing / multiplexing**, and
  combined GEX + antibody/VDJ libraries. Plain single-modality aligners will not demultiplex
  these; route here (and see **`bp-cite-seq`** downstream). Multiome **GEX + ATAC** →
  `--aligner cellrangerarc`.
- **kallisto | bustools (`--aligner kallisto`):** fastest, constant-memory pseudoalignment.
  Usable, but benchmarks flag that it reports **more low-content barcodes** and **spurious genes**
  (e.g. Olfr/Vmn families as likely mapping artifacts) (Brüning et al. 2022). Prefer only when
  speed dominates and you QC the cell/gene calls hard, or for a fast cross-check. ABA's Python
  equivalent is **`quantify-fastq-to-counts-kb`**.

## Decision features (what about the data drives the choice)

- **Droplet vs plate/other:** this knowhow is barcoded **droplet** scRNA (10x/Chromium and
  droplet variants). Plate-based (Smart-seq2/3) has no UMIs on the same footing → different route.
- **Assay chemistry / protocol:** 10x **v2 / v3 / v4** (and others) differ in barcode/UMI
  geometry — set `--protocol` (or `auto`); a mismatch silently wrecks barcode calling (Pitfalls).
- **Single-cell vs single-nucleus:** snRNA-seq is intron-rich → use an **augmented
  transcriptome** (simpleaf includes intronic sequence) or STARsolo **GeneFull**; a spliced-only
  reference undercounts nuclei.
- **Feature Barcode / multiplexing:** CITE-seq, hashing, or pooled/multiplexed libraries →
  **cellranger multi** (not a plain aligner). Multiome GEX+ATAC → **cellranger arc**.
- **RNA velocity downstream:** need spliced **and** unspliced counts → choose an augmented
  reference now (simpleaf / STARsolo Velocyto mode) and route to **`bp-rna-velocity`**.
- **Organism / annotation quality:** well-annotated (human/mouse) → transcriptome-based
  alevin-fry is fine; non-model / draft annotation or need for intergenic capture → genome-based
  STARsolo.
- **Need for genome BAMs:** coverage tracks, variant calling, novel transcripts → STARsolo /
  Cell Ranger, not pseudoalignment.
- **Reproducibility / compliance:** exact-Cell-Ranger requirement → Cell Ranger; otherwise the
  open-source simpleaf default is more portable and deterministic.
- **Compute budget & sample count:** many samples on limited compute → simpleaf/kallisto
  (lightweight); a few with full genome QC → STARsolo.

## Pitfalls

- **Wrong `--protocol` / chemistry.** Barcode+UMI length and whitelist are chemistry-specific;
  choosing v2 for a v3 run zeroes valid barcodes. Use `auto` or set the true chemistry, and check
  the barcode-rank/knee plot and AlevinQC.
- **Single-nucleus with a spliced-only reference.** Nuclei are dominated by unspliced/intronic
  reads — use augmented (simpleaf) or STARsolo GeneFull, or you throw away most signal.
- **kallisto's low-content cells + phantom genes.** Its default whitelist and multimapping
  behaviour over-call barcodes and inflate certain gene families — verify cell/gene calls or use it
  only as a cross-check.
- **Treating Feature-Barcode data as plain GEX.** CITE-seq/hashing/multiplexed libraries need
  **cellranger multi** to demultiplex; a single-modality aligner silently drops the tag/hash reads.
- **Cell Ranger container is community-built.** The nf-core `cellranger` container is not
  supported by 10x; for a licensed/validated Cell Ranger run use it deliberately (and confirm
  reference build), not as a casual default.
- **Empty-droplet over-filtering.** Aggressive knee/UMI thresholds delete real small/low-RNA
  cells; keep the **unfiltered** (`raw_feature_bc_matrix`) output and let cell calling
  (knee / **emptyDrops**, `--skip_emptydrops` to disable) plus biological QC in
  **`bp-quality-control`** decide — do not filter on UMI count alone here.
- **Reference/GTF drift.** The FASTA+GTF (or `--genome` iGenomes / `--cellranger_index`) used to
  build the index must match downstream annotation — pin one genome+GTF version.
- **Barcode (R1) vs biological (R2) reads.** For 10x, R1 carries barcode+UMI and R2 the cDNA;
  the samplesheet order matters, and FastQC on R1 will look "abnormal" by design.

## In ABA (which pipeline to run)

Run **nf-core/scrnaseq** via **`run_nextflow`** (verified live: latest release **4.1.0**,
nf-co.re/scrnaseq). Pick the aligner from the decision features above:

- **Aligner:** `--aligner simpleaf` (default, recommended) | `star` (STARsolo, genome/BAMs/snRNA)
  | `cellranger` (exact 10x output) | `cellrangermulti` (CITE-seq / hashing / multiplex /
  Feature Barcode) | `cellrangerarc` (multiome GEX+ATAC) | `kallisto` (fastest / cross-check).
- **Protocol:** `--protocol auto` (or set `10XV2` / `10XV3` / `10XV4` explicitly to match chemistry).
- **Reference:** `--genome` (iGenomes) or explicit `--fasta` + `--gtf`; `--cellranger_index` /
  `--transcript_fasta` where a prebuilt index applies. For snRNA / velocity choose an augmented
  reference (simpleaf) or STARsolo `--star_feature GeneFull`.
- **Samplesheet:** `sample,fastq_1,fastq_2` with the barcode read as `fastq_1` and cDNA as
  `fastq_2` per 10x convention.
- **Revision/profile:** pass `-r 4.1.0` for reproducibility and the cluster container profile
  (`singularity`/`apptainer`). Validate first with `-profile test`.

Call **`describe_pipeline`** for the full current parameter list (protocol strings, `--star_feature`,
`--skip_emptydrops`, index/reference options, multi/feature-barcode config) — per-param detail lives
there, not here.

**Hand-off:** take the resulting count matrix (`filtered_feature_bc_matrix` / h5ad; keep the
unfiltered `raw_feature_bc_matrix` for cell calling) into **`bp-quality-control`**, then the
single-cell bp-* chain (normalization → feature selection → dimensionality reduction → clustering →
annotation). CITE-seq/multimodal → **`bp-cite-seq`**; velocity → **`bp-rna-velocity`**.
