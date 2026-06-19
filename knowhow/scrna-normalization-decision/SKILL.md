---
name: scrna-normalization-decision
description: Decision guide for picking a scRNA-seq normalization — LogNormalize (shifted-log), SCTransform v2, analytic Pearson residuals, and where CLR fits (ADT only). Use when the user asks "which normalization should I use" or is choosing between scanpy's `normalize_total`+`log1p`, Seurat's `NormalizeData`, `SCTransform`, or `sc.experimental.pp.normalize_pearson_residuals`.
when_to_use: User is preparing scRNA-seq data for clustering/UMAP/DE/integration and has to commit to a normalization method, or is debugging a downstream artifact (over-mixed clusters, depth-correlated PCs, drop-out of rare cells) that points back to the normalization step. Trigger phrases include "which normalization", "LogNormalize vs SCTransform", "should I use Pearson residuals", "do I need SCTransform".
avoid_when: User wants the end-to-end runnable pipeline (use scrna-qc-clustering-v2 or seurat-scrna-v2); user is on CITE-seq / ADT (use bp-cite-seq — different statistics, CLR/DSB, not the methods here); user has bulk RNA-seq (size-factor / VST / rlog territory, not this doc).
invocation: interactive
kind: knowhow_draft
requires_tools: [WebFetch, Read]
keywords: [normalization, LogNormalize, SCTransform, SCTransform v2, Pearson residuals, analytic Pearson residuals, sctransform, log1p, normalize_total, scran, size factor, CLR, shifted log, depth normalization, glmGamPoi, NB regression, scRNA-seq]
domain: genomics
source: "Hafemeister & Satija 2019 (SCT v1); Choudhary & Satija 2022 (SCT v2 / Pearson residuals); Lause, Berens & Kobak 2021 (analytic Pearson residuals); Heumos et al. 2023 (Best Practices); Booeshaghi et al. 2022 (depth normalization dominance); Cole et al. 2019 (SCnorm comparison); recipe behavior of scrna-qc-clustering-v2, seurat-scrna-v2, bp-normalization."
audience: both
last_reviewed: null
reviewed_by: null
---

# scRNA-seq normalization — which method for which dataset?

The question this knowhow answers: **"Which normalization should I use for my
scRNA-seq data — LogNormalize, SCTransform, or Pearson residuals?"**

Audience: both naive (first-time biologist on a clustering pipeline) and
experienced (deciding whether to leave the default for a particular dataset).
The decision is multi-axis (median UMI per cell × cell-type-depth heterogeneity
× downstream task × speed budget). The executable recipes this knowhow points
at are `scrna-qc-clustering-v2` (scanpy + shifted-log), `seurat-scrna-v2`
(Seurat + LogNormalize, with an SCTransform alternative), and `bp-normalization`
(method-choice stage in the best-practices flow — implements both shifted-log
and analytic Pearson residuals; routes scran size factors too).

The most defensible empirical claim in the field is **not** "method X wins" but:
*depth normalization is the single largest determinant of downstream geometry;
the choice among shifted-log, SCT, and Pearson residuals is a second-order
refinement of that* (Booeshaghi et al. 2022). Read §2 with that calibration.

## Quick decision

For a typical 10x-style droplet dataset (3k–500k cells, median UMI 1k–10k,
human/mouse tissue):

| Your situation | Use | Recipe |
|---|---|---|
| First pass / unsure / fast iteration | **Shifted-log** (`normalize_total` + `log1p`, or Seurat `LogNormalize`) | `scrna-qc-clustering-v2` (scanpy) / `seurat-scrna-v2` (R) |
| Low median UMI (≲2k) **and** strongly heterogeneous cell-type depth (e.g. tissue with neurons + immune cells) | **SCTransform v2** | `seurat-scrna-v2` (uncomment the SCT alternative) |
| Feature selection / rare-cell detection is the *goal* | **Analytic Pearson residuals** | `bp-normalization` (path 3) |
| Downstream is **scVI / scANVI / totalVI / generative model** | **Skip normalization** — feed raw counts | `scvi-integration`, `scvi-totalvi-citeseq`, `scvi-de` |
| Downstream is **pseudobulk DE** (DESeq2/edgeR/limma) | **Skip normalization** — aggregate raw counts per sample, let the DE tool model depth | `bulk_rnaseq_de`, `deseq2_r`, `limma_voom` (after pseudobulking) |
| ADT / CITE-seq surface-protein modality | **CLR or DSB** (NOT the RNA methods on this row) | `bp-cite-seq`, `scvi-totalvi-citeseq` |

If you're not in the table: default to **shifted-log**. It is the workhorse,
underdetermined choices rarely hurt clustering at typical depths, and you can
swap to SCT or Pearson residuals later without re-running QC.

**Strong default rationale.** Shifted-log is what `scrna-qc-clustering-v2`,
`seurat-scrna-v2`, and the Heumos best-practices doc reach for first. It is
fast, deterministic, well-understood, and downstream tools (scanpy/Seurat
plotting, marker DE on log-normalized data) all expect it. SCT and Pearson
residuals are *better in known regimes* (see §4); they are not blanket
upgrades.

# REVIEW(scrna-normalization): is "shifted-log is the safe default at typical
depths" still the field consensus in 2025, or has the analytic-Pearson-residuals
path (Lause 2021 + Heumos 2023) overtaken it for first-pass clustering on
modern droplet data? Cite the most recent benchmark known.

## When NOT to normalize (or: when this decision doesn't apply)

- **You're feeding a count-model downstream.** scVI/scANVI/totalVI, pydeseq2
  on pseudobulk, glmGamPoi for marker DE — these consume *raw counts* and
  model depth internally. Normalizing first throws away the count statistic
  and can degrade the model. Keep raw counts in `adata.layers["counts"]` /
  the `RNA` assay's `counts` slot, and hand those to the count-model tool.
  (Heumos et al. 2023 — keep-raw-counts is a section heading, not a tip.)
- **You're doing pseudobulk DE.** Aggregate raw counts per sample × cell-type
  first; the bulk DE method (DESeq2/edgeR/limma-voom) handles depth. Running
  scRNA normalization and *then* pseudobulking double-corrects and can
  invert the variance structure.
- **You only have ADT / surface-protein counts.** Use `bp-cite-seq`. The
  methods in §4 below were developed against negative-binomial-like RNA
  count distributions; ADT counts are bimodal (background + specific) and
  need CLR or DSB. CLR on RNA is mathematically defined but loses
  information and is not the right tool — do not transplant it.
- **You're combining RNA across multiple samples for joint analysis.**
  Normalize per-sample (or jointly inside the integration recipe), then run
  Harmony / scVI / RPCA — the integration step is what aligns samples, not
  the normalization. See `scrna_pipeline.md` for the don't-`sc.concat`
  guardrail.

## Alternatives matrix

| Method | Characterization | When it wins | When it fails | Cost | Recipe |
|---|---|---|---|---|---|
| **Shifted-log (LogNormalize / `normalize_total`+`log1p`)** | Divide each cell by its size factor (default: median library size, or fixed 10k for "CP10k"), add 1, take log. The single-cell field's default for ~7 years. | Typical droplet data at moderate-to-high depth (≳3k UMI/cell), homogeneous-enough cell-type depth, downstream is PCA/UMAP/clustering or per-cluster marker DE. Fast, deterministic, every tool consumes log-norm input. (Heumos 2023) | Strongly depth-heterogeneous data (e.g. neurons at 50k UMI mixed with PBMCs at 2k) — the constant pseudo-count of 1 distorts low-count genes more in low-depth cells than high (Lause 2021); residual depth-correlation in PC1 is a common smell. | Seconds even on 1M cells; in-place on `.X`. | `scrna-qc-clustering-v2`, `seurat-scrna-v2`, `bp-normalization` path 1 |
| **SCTransform v2 (regularized NB + glmGamPoi)** | Per-gene regularized negative-binomial GLM on raw counts, with sequencing depth as the predictor; returns Pearson residuals as the "normalized" values and an adjusted count matrix. v2 (Choudhary & Satija 2022) re-parameterizes for stability and uses `glmGamPoi` for speed. | Heterogeneous depth across cell types **and** low median UMI (≲2k). Choudhary & Satija 2022 specifically report improved clustering separability + better marker stability vs v1 on low-depth + heterogeneous data; their headline benchmark vs LogNormalize is sharpest in that regime. | Atlas-scale data where the per-gene fit cost matters; data where SCT's regularization "smooths away" rare-cell signal you wanted to keep (less common since v2's reparameterization, but documented for v1). The "residual" output is unfamiliar to people expecting log-normalized counts; downstream marker DE has to be done correctly (see anti-patterns). | ~minutes on 50k cells with glmGamPoi; substantially slower without. R-only as a turnkey package (`Seurat::SCTransform`). | `seurat-scrna-v2` (the SCTransform alternative block) |
| **Analytic Pearson residuals (`sc.experimental.pp.normalize_pearson_residuals`)** | Closed-form Pearson residuals against a null NB(μ_ij, θ) where μ_ij = (gene_total × cell_total) / grand_total. No GLM fit per gene — analytic. (Lause, Berens & Kobak 2021) | Feature selection (HVG identification) — Lause 2021 shows analytic Pearson residuals select more informative genes than mean-variance-based HVG on shifted-log data; rare-cell-type detection benefits. Heumos 2023 recommends it specifically for HVG selection. Fast. | As the *input to downstream*, residuals can have negative values and are not log-counts — many tools and plots assume log-normalized expression and behave oddly on residuals. Best used in the HVG/PCA step, with shifted-log retained as the expression layer for plotting and DE. | Seconds even on 1M cells; cheaper than SCT because no per-gene optimization. | `bp-normalization` path 3 (scanpy `sc.experimental.pp.normalize_pearson_residuals`) |
| **CLR (centered log ratio)** | Compositional normalization: `log(x / geometric_mean(x))`. Mathematically appropriate for compositional data where only ratios are meaningful. | **ADT / surface-protein counts** (and similar bimodal bounded counts). Seurat uses CLR as the default ADT normalization. | **RNA counts.** Do not use on RNA. The compositional assumption is weaker for the RNA modality (you care about absolute expression for marker genes), and the field does not use CLR-on-RNA in practice. | Seconds. | `bp-cite-seq` (CLR / DSB); `scvi-totalvi-citeseq` (totalVI handles both modalities) |

# REVIEW(scrna-normalization): the "SCT v2 wins on low-depth + heterogeneous"
claim is sourced to Choudhary & Satija 2022 — confirm the specific median-UMI
threshold and the heterogeneity definition the paper uses; we're paraphrasing
as "≲2k median UMI and visibly different per-cell-type depth distributions",
which may understate or overstate. Pull the paper's Figure 2 / Table 1.

# REVIEW(scrna-normalization): for atlas-scale data (>1M cells, many tissues,
deep heterogeneity), is there a community consensus on normalization as of
2024-2025? Heumos 2023 names trade-offs without a winner; some atlas papers
have moved to feeding raw counts into scVI and skipping normalization
entirely for the embedding step. Confirm whether that pattern has hardened
into consensus.

## Anti-patterns and common mistakes

- **Using shifted-log output as input to a count model.** scVI/scANVI,
  pydeseq2 (after pseudobulking), glmGamPoi marker DE — all expect raw
  counts. Passing log-normalized values silently degrades the model.
  *Fix:* keep raw counts in `adata.layers["counts"]` (scanpy) or the
  `counts` slot of the `RNA` assay (Seurat) and pass those.

- **Running SCTransform and then `NormalizeData` again** on the same assay,
  or running both with DefaultAssay set to "SCT" while doing marker DE on
  the SCT residuals. The seurat-scrna-v2 recipe's SCT block is explicit:
  for marker DE switch `DefaultAssay` back to `"RNA"` (with `NormalizeData`
  run once), or use `PrepSCTFindMarkers`. Doing DE on raw SCT residuals
  without that step produces uninterpretable fold changes.

- **Normalizing per-sample and then `sc.concat`-ing, calling it
  "integration".** Normalization is not batch correction. Concatenation
  after per-sample shifted-log still gives you a joint embedding dominated
  by batch. Use a real integration recipe (`harmony-integration-scanpy`,
  `scvi-integration`, `seurat-integration`); see `scrna_pipeline.md`.

- **Switching to SCT because the dataset is "big".** SCT's advantages are
  about depth + heterogeneity, not size. On a large but
  uniformly-deep-and-typed dataset, shifted-log is the right choice; SCT
  just adds runtime without improving the clustering. The 1M-cell argument
  for SCT is not in the literature; the 1M-cell argument for *scVI* (skip
  normalization, model counts directly) is.

- **Pearson residuals as the expression layer for plotting / DE.** Residuals
  contain negative values and don't read like log-expression. Use Pearson
  residuals for HVG selection and as the PCA input; keep shifted-log values
  as the layer you display (`sc.pl.umap(adata, color="GENE")`) and as the
  input to marker DE. The `bp-normalization` recipe enforces this by
  computing both and storing them in named layers.

- **Picking SCTransform because "newer is better".** Choudhary & Satija 2022
  is real evidence, not a recency claim. The evidence is *regime-specific*
  (§4 row 2). On a typical PBMC-like dataset at ~5k median UMI, shifted-log
  and SCT v2 give similar clusters; SCT's cost only pays off when the
  preconditions hold.

- **Forgetting to re-run normalization after filtering.** If you filter
  cells (QC) *after* calling `NormalizeData` / `normalize_total`, the size
  factors were computed on the pre-filter set. Some downstream tools care;
  most don't. The Seurat / scanpy recipes both filter first, normalize
  second, which is the right order.

# REVIEW(scrna-normalization): "Pearson residuals are unfamiliar for
plotting" is a UX claim, not a methods claim — is there a benchmark where
plotting on residuals vs log-norm leads to different *interpretive* outcomes,
or is it purely an idiom mismatch? If purely idiom, soften the anti-pattern.

## Sanity checks — how to know your choice was right

After normalizing, before you commit downstream:

- **PC1 vs total counts per cell.** Plot total UMI (`adata.obs["n_counts"]`
  or `nCount_RNA`) against PC1 (and PC2). Strong correlation
  (|r| > ~0.4) is a residual depth artifact — your normalization didn't
  fully remove sequencing depth. With shifted-log on heterogeneous-depth
  data this is the canonical warning; switching to SCT v2 or analytic
  Pearson residuals usually drops the correlation.

- **HVG list reasonableness.** Look at the top 20–50 highly-variable genes.
  They should be biology (cell-type markers, immediate-early genes,
  stress-response, cell-cycle if you didn't regress it) — *not*
  dominated by ribosomal RPS/RPL genes or mitochondrial MT-* genes. If
  RPS/RPL/MT dominate, the normalization isn't pulling out biological
  variance well; this is often a Pearson-residuals-fixable failure for
  the HVG step (Lause 2021).

- **UMAP coloured by total counts per cell.** If cells gradient smoothly
  by depth across the UMAP rather than clustering by biology, the
  embedding is still depth-driven. Most often a QC problem (too-wide
  count distribution let through), but normalization choice contributes.

- **Marker genes match canonical biology.** Top marker per cluster should
  be a known marker for the cell type you'd guess from the UMAP geography
  (or at least a sensible candidate). If all clusters' top markers are
  ribosomal, you have a normalization + DE-method mismatch.

- **Switching method changes < ~10% of cluster assignments.** Run
  shifted-log; then run SCT v2 or Pearson residuals on the same upstream;
  cross-tab the cluster assignments (`pd.crosstab` /
  `table(seurat$shifted, seurat$sct)`). On a dataset where the
  preconditions for SCT/Pearson don't hold, the clusters should be ~the
  same — that's evidence the default was fine. Large reshuffles mean the
  dataset is in a regime where normalization matters; pick the method
  whose marker biology you trust more.

# REVIEW(scrna-normalization): the |r| > 0.4 threshold for PC1-vs-depth is
a rule-of-thumb I've seen used in practice but I don't have a cite for the
specific cutoff. Either find a citation or soften to "any visibly strong
correlation".

## See also

**Method papers (cited above; cite formats simplified):**

- Hafemeister & Satija 2019, Genome Biology — *Normalization and variance
  stabilization of single-cell RNA-seq data using regularized negative
  binomial regression* (SCT v1)
  (https://doi.org/10.1186/s13059-019-1874-1)
- Choudhary & Satija 2022, Genome Biology — *Comparison and evaluation of
  statistical error models for scRNA-seq* (SCT v2 + Pearson residuals
  discussion) (https://doi.org/10.1186/s13059-021-02584-9)
  # REVIEW(citations): confirm this is the SCT v2 paper and not a sibling
  paper from the Satija lab; the title above is paraphrased.
- Lause, Berens & Kobak 2021, Genome Biology — *Analytic Pearson residuals
  for normalization of single-cell RNA-seq UMI data*
  (https://doi.org/10.1186/s13059-021-02451-7)
- Booeshaghi et al. 2022 — *Depth normalization for single-cell genomics
  count data* (preprint / paper — the "depth normalization is the most
  important step" line). # REVIEW(citations): confirm the venue (bioRxiv vs
  journal) and final DOI.

**Best-practices / consensus:**

- Heumos et al. 2023, Nat Rev Genet — *Best practices for single-cell
  analysis across modalities* — the sc-best-practices.org book
  (https://doi.org/10.1038/s41576-023-00586-w;
  https://www.sc-best-practices.org/preprocessing_visualization/normalization.html)
- Cole et al. 2019, Cell Systems — *Performance assessment and selection
  of normalization procedures for single-cell RNA-seq* (SCnorm benchmark
  era; older but still cited)
  (https://doi.org/10.1016/j.cels.2019.03.010)

**Recipes that execute the paths above:**

- `scrna-qc-clustering-v2` — scanpy single-sample QC + clustering; uses
  shifted-log (`sc.pp.normalize_total` + `sc.pp.log1p`).
- `seurat-scrna-v2` — Seurat single-sample; default path uses
  `NormalizeData(method="LogNormalize")`; carries an explicit SCTransform
  alternative block.
- `bp-normalization` — best-practices normalization stage; implements
  shifted-log, scran pooling size factors, and analytic Pearson residuals
  with a downstream-task → method map.
- `bp-cite-seq` — CITE-seq / ADT normalization (CLR / DSB) — the CLR
  pointer from §2.
- `scvi-integration`, `scvi-de`, `scvi-totalvi-citeseq` — recipes that
  consume **raw counts** directly; do not pre-normalize.

**Adjacent knowhow:**

- `scrna_pipeline.md` — upstream framing (one sample vs many; do not
  `sc.concat`).
- `bulk_rnaseq_de.md` — pseudobulk DE path (skip scRNA normalization,
  aggregate raw counts, let DESeq2 handle depth).

**Recipe-coverage gaps (for the expert reviewer):**

- `# TODO(recipe): seurat-sctransform` — the assignment named a standalone
  `seurat-sctransform` recipe but the catalogue carries SCT only as the
  alternative block inside `seurat-scrna-v2`. Decide whether to split SCT
  into its own recipe (cleaner pointer) or keep it inlined (current
  state). Both are defensible; current state is what the matrix points
  at.
- `# TODO(recipe): scanpy-pearson-residuals-clustering` — analytic
  Pearson residuals are accessible only through `bp-normalization` (a
  *stage*, not a full pipeline). A user who wants "scrna-qc-clustering
  but with Pearson residuals for HVG" has to compose two recipes. A
  combined recipe would close this gap.
