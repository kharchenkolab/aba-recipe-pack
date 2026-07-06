---
name: scrna-de-methodology
description: Decision guide for differential expression in scRNA-seq — picks between per-cell tests (Wilcoxon / MAST / LR / ROC) and pseudobulk DESeq2/edgeR depending on whether the question is cluster markers, within-sample two-group, or cross-sample condition effects. Settles the per-cell-vs-pseudobulk question via Squair 2021. Use when the user asks "Wilcoxon or MAST or DESeq2 for my scRNA DE?" or "per-cell or pseudobulk?".
when_to_use: 'User has an scRNA-seq dataset (one sample or many) and wants to pick a DE test. Triggers: "DE in scRNA-seq", "per-cell vs pseudobulk", "Wilcoxon vs MAST vs DESeq2", "FindMarkers test.use", "condition DE on my Seurat/anndata object", "stim vs ctrl across donors".'
avoid_when: Bulk RNA-seq DE (use bulk-rnaseq-de / deseq2-r directly — no decision needed). Spatial / multiplexed-imaging DE (different statistical regime). Variant/allele-specific expression. The user wants the runnable recipe — point them at the recipes named in §4 + §7.
invocation: interactive
kind: knowhow_draft
requires_tools: [WebFetch, Read]
keywords: [scRNA-seq, differential expression, DE, Wilcoxon, MAST, DESeq2, edgeR, pseudobulk, pydeseq2, FindMarkers, FindAllMarkers, rank_genes_groups, scVI DE, cluster markers, condition DE, pseudoreplication, Squair, decision]
domain: genomics
source: "Squair et al. 2021, Nat Commun (https://doi.org/10.1038/s41467-021-25960-2); Crowell et al. 2020, Nat Commun (muscat); Soneson & Robinson 2018, Nat Methods; Finak et al. 2015 (MAST); Love et al. 2014 (DESeq2); Heumos et al. 2023 single-cell best-practices (sc-best-practices.org/conditions/differential_gene_expression.html)."
audience: both
produces: []
capabilities_needed: []
---

# scRNA-seq differential expression — per-cell or pseudobulk? Wilcoxon, MAST, or DESeq2?

The question this knowhow answers: **"I want to test differential expression
in scRNA-seq — should I do it per-cell or pseudobulk? Wilcoxon, MAST, or
DESeq2?"**

Audience: both naive (first DE analysis on a single-cell object — usually
asking implicitly "which `test.use` do I pass to `FindMarkers`?") and
experienced (knows the alternatives exist, wants to defend a choice in a
methods section).

The decision shape is **a tree, not a continuum**: the right method depends
on WHICH QUESTION you are actually asking, not on a dataset feature. The
three questions get three different answers, and the most common mistake is
applying the answer to one question against another.

The executable recipes this knowhow connects to: `scrna_qc_clustering`
(`rank_genes_groups`, Python/scanpy), `seurat-scrna`
(`FindAllMarkers` / `FindMarkers`, R/Seurat), `bp-differential-expression`
(pseudobulk aggregation + DESeq2/edgeR/pydeseq2, scanpy + decoupler),
`deseq2-r` (the R DESeq2 backend on pseudobulk counts), `bulk-rnaseq-de`
(pydeseq2 backend on pseudobulk counts), and `scvi-de` (scVI Bayesian DE
on a trained model).

## Quick decision

Before picking a method, **answer this question first**: what are you
comparing?

| Your question | Use | Recipe |
|---|---|---|
| **Q1. Cluster markers (one cluster vs the rest) within one dataset / object** | Per-cell **Wilcoxon** (Seurat `FindAllMarkers` default; scanpy `sc.tl.rank_genes_groups(method='wilcoxon')`) | `seurat-scrna` (Seurat) / `scrna_qc_clustering` (scanpy) |
| **Q2. Two groups WITHIN one sample, with per-cell covariates** (e.g. cells expressing gene X vs not, controlling for %MT / nCount) | **MAST** or **LR** with `latent.vars` (Seurat `FindMarkers(test.use='MAST', latent.vars=…)`) | `seurat-scrna` (test.use=MAST/LR) — # REVIEW(recipes): no dedicated Seurat-DE recipe yet |
| **Q3. Condition effect across MULTIPLE samples** (stim vs ctrl, disease vs healthy — multiple donors per group) | **Pseudobulk** (sum raw counts per sample × cell type) then **DESeq2 / edgeR** | `bp-differential-expression` then `deseq2-r` or `bulk-rnaseq-de` |
| Q3 but you already have a trained scVI/scANVI | scVI Bayesian DE (`model.differential_expression`) — second-line | `scvi-de` |

**The hard rule (load-bearing):** for **Q3 (cross-sample condition DE)**,
bulk DE methods (DESeq2, edgeR, limma, pydeseq2) applied **directly to a
per-cell count matrix are wrong** — they inflate the false positive rate by
roughly **10–100×** because they treat each cell as an independent biological
replicate (pseudoreplication). The unit of replication is the **subject /
sample**, not the cell. Aggregate to pseudobulk first; then run the bulk
method on the pseudobulk matrix. This is settled methodology — Squair et al.
2021, Nat Commun (https://doi.org/10.1038/s41467-021-25960-2) demonstrated it
on multiple ground-truth datasets, Crowell 2020 (muscat) reached the same
conclusion, and Heumos 2023 (single-cell best-practices) codifies it.
# REVIEW(scrna-de): is "10-100×" the right magnitude to quote? Squair's
figures show inflation depends on the cells-per-subject ratio — give the
range with the regime that produces each end.

If your question isn't in the table, OR your dataset shape is unusual
(one sample with no replicates and you want condition DE, or paired
within-subject design, or repeated measures), **read §4 (alternatives
matrix) before picking** — the choice changes.

## When NOT to do this analysis

- **One sample, no replicates, but you want a condition effect.** You
  cannot do condition DE without replication — there is no statistical
  basis for separating biology from sample variability. The honest answer
  is "the experiment doesn't support this comparison." If forced, run
  cluster markers (Q1, Wilcoxon) and call them "differential between
  populations within this one sample", not "differential between conditions".
  # REVIEW(scrna-de): is this the right framing for the n=1-per-arm case?
  Some labs use MAST with a subject random effect even at n=1; Squair
  argues this still inflates FPR.
- **You only have a gene list, no count matrix.** DE needs counts. Use
  enrichment (GSEA / ORA) on the gene list — see `bp-gsea-pathway`.
- **You want to call cell types, not test genes.** "Differential expression"
  is sometimes used loosely to mean "find markers for annotation". For
  annotation, use the cluster-marker workflow (Q1) + a reference-mapping
  tool (`SingleR`, `Azimuth`, `scvi-reference-mapping`), NOT a condition
  DE pipeline.
- **You haven't run QC + clustering yet.** DE on un-QC'd, un-clustered
  cells will be dominated by quality artifacts (dying cells, doublets).
  Run `scrna_qc_clustering` / `seurat-scrna` first.
- **You want to "test the integration".** Integration is corrected at
  the embedding/latent level; DE should be done on the **raw counts**
  (with the integrated labels as the grouping), NOT on the integrated
  expression values. Don't run DE on `X_scVI` / Harmony-corrected counts.

## Alternatives — full matrix

Every method below is a valid answer to **at least one** of Q1/Q2/Q3, but
they answer **different questions**. The "When it wins" column names the
question; the "When it fails" column names what happens if you apply it to
the wrong question.

| Method | Characterization | When it wins | When it fails | Cost | Recipe |
|---|---|---|---|---|---|
| **Wilcoxon (per-cell)** | Rank-sum test, per gene, on log-normalized counts; non-parametric | **Q1 cluster markers** — the default. Seurat `FindAllMarkers` and scanpy `rank_genes_groups` defaults. Fast (with `presto`, ~seconds on 100k cells). | **Q3 condition DE** — every cell treated as independent → 10-100× inflated FPR (Squair 2021). Also weak when per-cell covariates matter. | Linear in cells × genes; fast with `presto` | `seurat-scrna`, `scrna_qc_clustering` |
| **MAST (per-cell)** | Two-part hurdle model — logistic for detection + Gaussian for log-expression — accommodates zero-inflation, supports covariates as `latent.vars` and a subject random effect | **Q2 within-sample two-group with per-cell covariates** (%MT, nCount, cell-cycle score). Also acceptable for **Q3 IF a subject random effect is included** (Heumos 2023). | Plain MAST without a subject term → same pseudoreplication failure as Wilcoxon for Q3 (Squair 2021). Slower than Wilcoxon. | Minutes on a typical object; slower with random effects | `seurat-scrna` (test.use='MAST', latent.vars=…); # REVIEW(recipes): no dedicated MAST + random-effect recipe |
| **LR / Logistic regression (per-cell)** | Per-gene logistic regression of group membership on expression, with `latent.vars` as covariates | **Q2 within-sample two-group with covariates** — when MAST is too slow / fragile and a simpler covariate-aware test suffices | Same Q3 failure as Wilcoxon (pseudoreplication). LR returns a Wald p-value but no fold change directly. | Minutes; comparable to MAST | `seurat-scrna` (test.use='LR', latent.vars=…) |
| **ROC (per-cell)** | Per-gene AUC of the binary group label using gene expression as a classifier — descriptive, not inferential | **Q1 cluster markers when you want a "cleanness" score per marker**, not a p-value. Good for picking the visually crispest marker. | Returns no p-value → cannot be used for significance reporting. Not for Q3. | Fast | `seurat-scrna` (test.use='roc') |
| **Pseudobulk + DESeq2 (R)** | Sum raw counts per (sample × cell type) → bulk count matrix → DESeq2 negative-binomial GLM with sample-level design. Supports LRT, multi-factor, interactions, custom contrasts. | **Q3 cross-sample condition DE** — the canonical correct path (Squair 2021, Crowell 2020, Heumos 2023). Mandatory when you need LRT / interactions / arbitrary contrasts. | Needs ≥2-3 biological replicates per condition. Fails on n=1 per arm (no dispersion estimate). NEVER on per-cell counts directly. | Seconds-to-minutes on a pseudobulk matrix | `bp-differential-expression` (aggregation) → `deseq2-r` |
| **Pseudobulk + edgeR (R)** | Pseudobulk matrix → edgeR negative-binomial GLM. Comparable to DESeq2; some benchmarks (Crowell 2020) put it marginally ahead for complex designs / few replicates. | **Q3 with complex designs** (multiple factors, low replicate count) where edgeR's empirical Bayes shrinkage is robust. | Same replicate / pseudobulk requirements as DESeq2. NEVER on per-cell counts. | Comparable to DESeq2 | `bp-differential-expression` then edgeR — # REVIEW(recipes): no dedicated edgeR recipe in catalogue; deseq2-r is the closest, edgeR call is a 5-liner inside `run_r` |
| **Pseudobulk + pydeseq2 (Python)** | Owkin's Python re-implementation of DESeq2; Wald-only (no LRT). Operates on a pseudobulk matrix in Python. | **Q3 in a Python-native session**, when LRT/multi-factor isn't needed | NEVER on per-cell `AnnData.X` directly (the most common anti-pattern — see §5). No LRT, fewer features than `deseq2-r`. | Seconds on pseudobulk | `bp-differential-expression` → `bulk-rnaseq-de` |
| **scVI DE** | Posterior-based Bayesian DE on a trained scVI/scANVI model; reports `proba_de`, `bayes_factor`, `lfc_mean`; can marginalize over batch | **Q3 when you already have a trained scVI model and a heavily batched dataset.** Atlas-scale, batch structure baked in. | Requires a trained model (training is expensive). No explicit subject random effect — Squair 2021 argues scVI DE alone is still prone to pseudoreplication; combine with per-sample resampling or use as a sanity check against pseudobulk DESeq2. | GPU-hours for training; DE itself is fast | `scvi-de` (requires upstream `scvi-integration`) |

For per-method internals (the MAST hurdle parameterization, DESeq2's
size-factor estimation, scVI's posterior sampling), see the method papers
in §7. Per-recipe operational detail lives in the recipes themselves.

## Anti-patterns and common mistakes

Five common mistakes — **the first four all reduce to the same root
cause**: using a bulk method on per-cell counts. Squair 2021 is explicit
that this is the dominant failure mode in published scRNA DE.

- **`DESeqDataSetFromMatrix(per_cell_counts, design = ~condition)` directly
  on a per-cell matrix.** The canonical wrong thing. Each cell becomes a
  pseudo-sample, dispersion is estimated on cell-level variance, and the
  Wald test sees ~10,000s "samples" → every gene is significant. Treats
  cells as biological replicates → pseudoreplication. **Aggregate to
  pseudobulk per sample × cell type first**; then call
  `DESeqDataSetFromMatrix` on the pseudobulk matrix (rows = pseudo-samples
  = per-donor sums; design = `~ condition` with donor-level variability).
  (Squair 2021, Crowell 2020, Heumos 2023). See `bp-differential-expression`
  / `deseq2-r`.
- **`FindMarkers(obj, test.use = "DESeq2", group.by = "condition")` directly
  on a per-cell Seurat object for a CONDITION effect.** Seurat exposes
  `test.use = "DESeq2"` and it runs, but it runs DESeq2 with cells-as-samples
  → exactly the failure above. The `test.use="DESeq2"` flag is documented but
  is **not** the right tool for cross-sample condition DE. For Q3, leave
  `FindMarkers` entirely — aggregate (e.g. Seurat v5 `AggregateExpression(obj,
  group.by=c("sample","cell_type"))`) and run DESeq2 on the resulting bulk
  matrix.
- **`limma::lmFit` on the log-normalized per-cell matrix for a condition
  test.** limma was designed for bulk + per-sample replicates. Applied to a
  per-cell `logcounts` matrix with cells as columns, every cell is a "sample"
  → same pseudoreplication failure. Use limma-voom on pseudobulk counts
  (see `limma_voom` recipe) for Q3.
- **`pydeseq2.DeseqDataSet(counts=adata.X.T, metadata=adata.obs,
  design_factors="condition")` on a per-cell AnnData.** Same failure in
  Python. pydeseq2 will happily fit; the result is wrong. Pseudobulk first
  (see `bp-differential-expression`), then `bulk-rnaseq-de`.
- **Reporting Wilcoxon p-values for a CONDITION contrast across donors as if
  they were "the DE result"**, then complaining the FDR is "too good to be
  true". The p-values are correct under H0 of "this cell is exchangeable",
  but cells from one donor are NOT exchangeable with cells from another — the
  null is wrong. Re-do as pseudobulk DESeq2 and expect ~10× fewer hits.
  (Squair 2021).

Other common mistakes that aren't variants of the above:

- **Running DE on integrated / batch-corrected expression values
  (Harmony-corrected counts, `X_scVI`, ComBat-adjusted).** Integration adjusts
  the EMBEDDING, not the counts. DE on adjusted expression destroys the
  count-variance relationship that DESeq2/edgeR/Wilcoxon assume. **Use the
  raw counts; carry the integrated cluster labels and batch covariate
  separately into the design / model.**
- **Pseudobulk with too few cells per (sample × cell type) bin.** A
  pseudobulk row built from 3 cells is unreliable. Filter to ≥30 cells /
  pseudo-sample (Heumos 2023); drop pseudo-samples below the threshold;
  drop cell types where >half the samples are under-filled. Filter rule:
  `min_cells=30, min_counts=1000` in `decoupler.get_pseudobulk` is a sane
  default. # REVIEW(scrna-de): is 30 cells the current consensus threshold,
  or has it tightened?
- **Cluster-vs-cluster DE within a single sample, then calling DESeq2 with
  n=1 per group.** Pseudobulking two clusters within ONE sample makes
  n=1 per arm — DESeq2 has no dispersion estimate and the statistics are
  meaningless. Use Wilcoxon (Q1) for within-sample cluster comparisons;
  reach for pseudobulk DESeq2 only when there are real replicates.
- **Using `padj < 0.05` alone as the cutoff on very large objects.**
  With 100k+ cells and a per-cell test, the standard FDR threshold leaves
  thousands of "significant" hits driven by tiny effect sizes. Combine FDR
  with a fold-change threshold (`|log2FC| ≥ 0.25` for markers, `≥ 0.5–1.0`
  for condition DE).
- **Using SCT-normalized data for DE without `PrepSCTFindMarkers`.** Seurat
  v5's SCTransform stores per-cell residuals that aren't directly
  comparable across cells. Either `PrepSCTFindMarkers` first, or do DE
  on the `RNA` assay's log-normalized counts (the older default), not the
  `SCT` assay.

## Sanity checks — how to know your choice was right

After running the chosen DE, verify:

- **The number of hits is plausible.** A per-cell condition test (wrong)
  typically returns thousands of "significant" genes with weak fold
  changes; a pseudobulk DESeq2 of the same comparison returns dozens to
  low hundreds. If you got 5000+ significant genes from a 2-condition
  comparison, you are almost certainly looking at pseudoreplication
  inflation — re-do as pseudobulk. (Squair 2021).
- **Volcano plot has a visible "cap" or banana shape.** Pseudobulk DESeq2
  on real biology produces a clear volcano; per-cell DESeq2 produces a
  near-vertical wall (huge `-log10(padj)`) that should be a red flag.
- **Per-sample PCA / MDS of the pseudobulk matrix** (Q3): samples within a
  condition should cluster together, samples across conditions should
  separate. If one donor is far from the others within its condition,
  inspect — may be a batch / quality outlier and a candidate for
  exclusion or for entry into the design as a covariate.
- **Stability across resampling.** Drop one sample at a time, re-run DE,
  count overlap of top-100 hits. Pseudobulk results are usually stable
  (>70% overlap); per-cell results on the wrong question are unstable
  (often <20% overlap). Cheap diagnostic for inflation. # REVIEW(scrna-de):
  do you have a citation for the 70%/20% rule of thumb, or is this
  community lore?
- **Top hits are biologically interpretable.** Cross-condition pseudobulk
  hits in a stim vs ctrl PBMC dataset should be dominated by interferon
  response genes (ISG15, IFIT1/2/3, OAS1, IFI6). Cluster markers in a
  PBMC analysis should be the canonical lineage markers
  (CD3D/CD8A/MS4A1/LYZ/…). If neither pattern shows, suspect the wrong
  comparison or upstream QC issues.
- **You can defend the choice in a methods section.** A reviewer asking
  "why DESeq2 on pseudobulk rather than per-cell Wilcoxon?" should be
  answerable in one sentence ("cells from one donor are not independent
  replicates; per-cell DE inflates FDR ~10× — Squair 2021").

## See also

**The benchmark / canonical paper (the position-setter):**

- Squair et al. 2021, Nat Commun — "Confronting false discoveries in
  single-cell differential expression"
  (https://doi.org/10.1038/s41467-021-25960-2). Demonstrates per-cell DE
  methods inflate FPR by an order of magnitude on cross-condition tests
  with biological replicates; pseudobulk methods (DESeq2, edgeR, limma)
  on aggregated counts give well-calibrated FDR. THE canonical paper for
  this decision.

**Supporting benchmarks + best-practices:**

- Crowell et al. 2020, Nat Commun — muscat — pseudobulk method comparison
  across DESeq2, edgeR, limma, and mixed models
  (https://doi.org/10.1038/s41467-020-19894-4). Reaches the same
  conclusion: pseudobulk DESeq2/edgeR are the practical defaults for
  cross-condition DE. # REVIEW(scrna-de): verify the DOI; reading from
  memory.
- Soneson & Robinson 2018, Nat Methods — earlier scRNA DE benchmark, before
  the pseudobulk consensus
  (https://doi.org/10.1038/nmeth.4612). # REVIEW(scrna-de): confirm DOI;
  per-cell methods looked better here than under Squair 2021, but the
  ground truth set was different.
- Heumos et al. 2023, Nat Rev Genet — single-cell best-practices,
  condition-DE chapter
  (https://sc-best-practices.org/conditions/differential_gene_expression.html).
  The community consensus statement: pseudobulk + DESeq2/edgeR/limma is the
  default; MAST with a subject random effect is acceptable; per-cell tests
  without a subject term are pseudoreplication.

**Method papers:**

- Finak et al. 2015, Genome Biology — MAST
  (https://doi.org/10.1186/s13059-015-0844-5). Two-part hurdle model;
  the canonical per-cell test for Q2.
- Love et al. 2014, Genome Biology — DESeq2
  (https://doi.org/10.1186/s13059-014-0550-8). The bulk DE workhorse; on
  pseudobulk it is the recommended cross-condition method.
- Robinson et al. 2010, Bioinformatics — edgeR
  (https://doi.org/10.1093/bioinformatics/btp616). Pseudobulk DE alternative
  to DESeq2.
- Ritchie et al. 2015, NAR — limma / voom
  (https://doi.org/10.1093/nar/gkv007). Pseudobulk DE alternative
  optimized for small sample counts.
- Hafemeister & Satija 2019, Genome Biology — SCTransform / SCT residuals
  (https://doi.org/10.1186/s13059-019-1874-1). Relevant to the
  SCT-assay-DE caveat in §5.
- Lopez et al. 2018, Nat Methods — scVI
  (https://doi.org/10.1038/s41592-018-0229-2). Trained-model DE
  (`model.differential_expression`).

**Recipes that execute paths from §2 / §4:**

- `seurat-scrna` — Seurat single-sample flow, includes `FindAllMarkers`
  (Q1, Wilcoxon default; `test.use='MAST'/'LR'/'roc'/'DESeq2'/'negbinom'`
  variants for Q2). Path A (Wilcoxon) and Path B (MAST/LR with
  `latent.vars`) live here.
  # REVIEW(recipes): should we factor out a dedicated `seurat-de-testing`
  recipe that organizes Paths A/B/C as the prompt's structure suggests?
  Currently they're embedded inside `seurat-scrna` Step 7.
- `scrna_qc_clustering` — scanpy single-sample flow, includes
  `sc.tl.rank_genes_groups` (Q1, Wilcoxon default; t-test / logreg
  alternatives).
- `bp-differential-expression` — the canonical Q3 pipeline: pseudobulk
  aggregation (decoupler) → pydeseq2 / DESeq2 / edgeR. Handoff step
  before `deseq2-r` or `bulk-rnaseq-de`.
- `deseq2-r` — R DESeq2 on a bulk or pseudobulk count matrix. The
  Q3 backend when LRT / multi-factor / arbitrary contrasts are needed.
- `bulk-rnaseq-de` — pydeseq2 on a bulk or pseudobulk count matrix.
  The Q3 backend in a Python-native session (Wald-only).
- `limma_voom` — limma + voom on a bulk or pseudobulk count matrix.
  Q3 alternative, especially for small sample counts.
- `scvi-de` — Bayesian DE on a trained scVI/scANVI model. Q3 sanity check
  / atlas-scale DE.

**Coverage gaps surfaced during authoring:**

- `# TODO(recipe): seurat-de-testing` — a dedicated Seurat DE recipe
  factoring Path A (Wilcoxon cluster markers), Path B (MAST/LR with
  `latent.vars` for Q2), and Path C (Seurat v5 `AggregateExpression`
  → DESeq2 on pseudobulk for Q3) would make this knowhow's pointers
  resolve cleanly. Currently §2 and §4 hedge with "embedded in
  `seurat-scrna` Step 7".
- `# TODO(recipe): edgeR-r` (or fold into `deseq2-r`) — the alternatives
  matrix lists edgeR but the catalogue has DESeq2 + pydeseq2 + limma_voom
  only; edgeR is a one-liner in `run_r` but no dedicated recipe.
- `# TODO(recipe): scrna-pseudobulk-de-r` — the R equivalent of
  `bp-differential-expression` (the Python pipeline). Aggregate inside R
  (`Seurat::AggregateExpression`) → DESeq2/edgeR — a natural companion
  to `seurat-scrna`.

**Adjacent knowhow:**

- `scrna_pipeline.md` — upstream framing (one sample → keep separate,
  multiple samples → integrate first OR do per-sample then aggregate).
  The "keep samples separate" principle there is what enables pseudobulk
  Q3 here.
- `bulk_rnaseq_de.md` — sibling reference for the bulk DE pipeline that
  pseudobulk feeds into.
- `scrna-integration-knowhow` (if drafted) — choice of integration method
  upstream of Q3 DE; integration is at the embedding layer, DE is on raw
  counts. # REVIEW(coverage): is the integration knowhow drafted yet?
