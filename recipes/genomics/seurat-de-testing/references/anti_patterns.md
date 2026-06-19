# Anti-patterns — what NOT to do for scRNA-seq DE

The bulk-DE-tools-on-per-cell-counts mistake — why it's wrong, why it's
common, what to do instead. This is the load-bearing reference behind the
"hard rule" in SKILL.md's framing.

Load this reference when the user (or an upstream code path) is about to
run DESeq2 / edgeR / limma / pydeseq2 on a per-cell matrix, or when MAST is
being used with `latent.vars = "donor"` as a condition test. The failure
mode in both cases is the same: anti-conservative p-values, inflated false
positives, hits that don't replicate. The fix in both cases is the same:
pseudobulk to the sample level first.

## The headline rule

> **Bulk DE tools require BIOLOGICAL REPLICATES as input columns.** A
> "biological replicate" is an independent sample (donor, animal, mouse,
> well) — NOT a single cell. Cells from the same donor are correlated:
> they share genotype, environment, technical batch, condition history.
> Treating thousands of cells from one donor as thousands of independent
> replicates **inflates the false-positive rate 10–100×** (Squair et al.
> 2021, Nat Commun, DOI 10.1038/s41467-021-25960-2).
>
> The correct uses of bulk DE tools (DESeq2, edgeR, limma, pydeseq2) on
> single-cell data:
>
> 1. **Pseudobulk** — `AggregateExpression(return.seurat = T,
>    group.by = c("donor_id", "condition", "celltype"))` first; THEN the
>    bulk tool on the sample-level matrix. The aggregated columns are
>    valid biological replicates.
> 2. **Bulk RNA-seq input** (not single-cell at all) — one column per
>    sample, raw counts. Use `deseq2-r` / `bulk_rnaseq_de` recipes.
>
> Anything else with these tools is wrong.

## Why this is rule zero — the Squair 2021 result

Squair et al. 2021 (Nat Commun, "Confronting false discoveries in
single-cell differential expression") benchmarked 14 DE methods on 18
ground-truth scRNA-seq datasets where the "answer" was known from matched
bulk RNA-seq. The finding:

- Per-cell DE methods (Wilcoxon, MAST, edgeR-LRT on cells, etc.) reported
  **2–100× more "significant" genes than the matched bulk DE**, even
  after FDR adjustment. The excess hits do not replicate when more donors
  are added.
- Pseudobulk methods (DESeq2, edgeR, limma-trend on aggregated counts)
  reported hit counts and FDRs **consistent with the matched bulk
  ground truth**.
- The FPR inflation in per-cell methods is driven by treating cells as
  replicates when they are not. Adding more cells from the same donor
  does NOT add statistical power for cross-donor inference — it only
  shrinks the p-value of the per-cell null, which is the wrong null.

Soneson & Robinson 2018 (Nat Methods, "Bias, robustness and scalability in
single-cell differential expression analysis") is the earlier complementary
benchmark on 36 methods that reached a similar conclusion for the
WITHIN-sample comparison (per-cell methods are calibrated within sample,
not across donors).

## Anti-pattern 1 — DESeq2 directly on per-cell counts

```r
# WRONG
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")   # cells × genes
coldata <- obj@meta.data[, c("stim", "donor_id")]
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                       colData = coldata,
                                       design  = ~ stim)
dds <- DESeq2::DESeq(dds)
```

This call WILL run without erroring. DESeq2's machinery treats each cell
as a column (replicate), fits dispersion across cells, computes Wald
statistics. The output looks legitimate — there's an `padj`, there are
"significant" hits. Almost all of them are false positives by the
bulk-RNA-seq ground truth.

**Why it's wrong:** cells from one donor are not independent replicates of
that donor's biology. DESeq2's dispersion model assumes biological
replicates contribute independent variance estimates. With thousands of
cells per donor, the dispersion shrinks toward the technical (intra-donor)
variance, not the biological (cross-donor) variance — so the test is
effectively asking "are there any differences between conditions in any
cell?" rather than "is the donor-level mean shifted by condition?".

**Right answer:** Path C in SKILL.md.

```r
# RIGHT
pseudo <- AggregateExpression(obj, assays = "RNA", return.seurat = TRUE,
                              group.by = c("stim", "donor_id", "seurat_annotations"))
# Now build DESeqDataSet from the pseudobulk matrix; see pseudobulk_workflow.md
# for the full multi-factor design.
```

## Anti-pattern 2 — `FindMarkers(obj, test.use = "DESeq2")` on a per-cell object for a condition effect

Same failure mode, dressed in a friendlier API.

```r
# WRONG — `obj` here is the per-cell Seurat object
de <- FindMarkers(obj, ident.1 = "STIM", ident.2 = "CTRL",
                  test.use = "DESeq2", group.by = "stim")
```

Seurat's `DESeq2DETest` (verified via source) builds a
`DESeqDataSetFromMatrix(countData = data.use, design = ~group)` where
`data.use` is the per-cell counts of the cells in ident.1 + ident.2. Same
as Anti-pattern 1, just via the wrapper.

**Right answer:** invoke `FindMarkers(test.use = "DESeq2")` on the
**pseudobulk** Seurat object (after `AggregateExpression(return.seurat =
TRUE)`), where each column is a donor × condition × cell-type pseudobulk.

## Anti-pattern 3 — limma / limma-voom on log-normalized per-cell expression

```r
# WRONG
expr <- GetAssayData(obj, assay = "RNA", layer = "data")       # cells × genes, log-normalized
design <- model.matrix(~ stim, data = obj@meta.data)
fit <- limma::lmFit(expr, design)
fit <- limma::eBayes(fit)
top <- limma::topTable(fit, coef = "stimSTIM", number = Inf)
```

limma assumes each column is a biological replicate of a sample. Same
inflation as Anti-pattern 1.

**Right answer:** pseudobulk first, then `limma::voom` + `lmFit` on the
sample-level matrix (limma-trend is also fine on log-CPM pseudobulk).
limma's `voom` is well-tuned for the sample-level mean-variance shape of
pseudobulk counts.

## Anti-pattern 4 — pydeseq2 on `.X` of an AnnData (cells × genes)

```python
# WRONG
import pydeseq2
from pydeseq2.dds import DeseqDataSet
ds = DeseqDataSet(counts = adata.X.T.toarray(),     # cells become "samples"
                  metadata = adata.obs[["stim"]],
                  design_factors = "stim")
ds.deseq2()
```

pydeseq2 is the Python port of DESeq2 — same model, same assumptions,
same wrong-when-per-cell behavior. AnnData `.X` is the cells × genes
matrix; transposing to feed pydeseq2 with cells as "samples" reproduces
Anti-pattern 1 in Python.

**Right answer:** pseudobulk in Python first.

```python
# RIGHT (sketch)
import scanpy as sc
pb = sc.get.aggregate(adata, by = ["donor_id", "stim", "cell_type"], func = "sum")
# pb is now an AnnData whose .X is the pseudobulk matrix; columns are
# donor × stim × cell_type. Then run pydeseq2 on pb.X.T.
```

(For the canonical Python pseudobulk recipe, see ABA's scanpy DE recipes.)

## Anti-pattern 5 — MAST with `latent.vars = "donor"` as a condition test

```r
# WRONG (when the question is a cross-donor condition effect)
de_mast <- FindMarkers(obj, ident.1 = "STIM", ident.2 = "CTRL",
                       group.by = "stim",
                       test.use = "MAST",
                       latent.vars = "donor_id")
```

MAST's `zlm` accepts `latent.vars` as DESIGN-MATRIX columns — they are
treated as per-cell covariates of the cell's expression, NOT as a random
effect that bumps the effective sample size from N_cells to N_donors. The
test is still per-cell; the p-values are still anti-conservative for
cross-donor inference.

This is one of the most common live-session mistakes because the API
surface ("just add latent.vars = 'donor'") looks like it should account
for donor-level variance. It doesn't. Squair 2021 explicitly tested MAST
with `latent.vars = donor` and showed it remained anti-conservative
vs the bulk ground truth.

**Right answer:** Path C. If the question is "stim vs ctrl across
donors", pseudobulk first. MAST is the right test for WITHIN-DONOR
comparisons with per-cell covariates (cell-cycle phase, percent.mt) — not
a substitute for sample-level aggregation.

A more advanced alternative is a **mixed model** (`glmer.nb` or
`muscat::mmDS`) that puts donor as a random effect — that does
properly account for cross-donor variance, but it's substantially more
complex than pseudobulk and Squair 2021 shows the two approaches reach
similar calibration. Pseudobulk wins on simplicity.

## Anti-pattern 6 — n = 1 or n = 2 donors per condition, then "robust" DE

```r
# WRONG — only 2 donors per condition
# (Even with correct pseudobulk methodology)
table(obj$donor_id, obj$stim)
#         STIM CTRL
#   d1    1000  0
#   d2    1200  0
#   d3       0 1100
#   d4       0  900
de <- FindMarkers(pseudo, ident.1 = "..._STIM", ident.2 = "..._CTRL",
                  test.use = "DESeq2")
```

DESeq2 will run with 2 donors per side, but its dispersion estimate is
unstable and its p-values are unreliable. With n = 1 per condition, there
is NO replication for DE — only descriptive comparison.

**Right answer:**
- If `min(n_per_condition) < 3`: **fall back to description**, not a
  test. Report the direction of change in mean expression with a CAVEAT
  that replication is insufficient for inference.
- If `min(n_per_condition) ≥ 3` but borderline (3–4): run the test, but
  flag the replication count in the report so the reader sees the
  fragility.

## Anti-pattern 7 — running per-cell DE THEN pseudobulk and combining

```r
# WRONG
de_cell <- FindMarkers(obj, ...,  test.use = "wilcox")    # per-cell
de_pb   <- FindMarkers(pseudo, ..., test.use = "DESeq2")  # pseudobulk
common  <- intersect(de_cell$gene, de_pb$gene)            # report both?
```

The "intersection of per-cell and pseudobulk" reasoning is appealing
("they agree, so it's robust") but it is **not** a statistical
combination — it's just two tests with different nulls. The pseudobulk
result is the calibrated one; per-cell adds nothing to the cross-donor
inference. Report pseudobulk; ignore per-cell for the condition question.

Per-cell DE is fine to report SEPARATELY as a within-sample question
("of the genes pseudobulk-significant for stim vs ctrl, which are also
heterogeneously expressed across cells within the stim condition?") —
but that's a different question, not a robustness check.

## What to do when you can't pseudobulk

Some legitimate scenarios where Path C doesn't apply:

| Scenario | What's available | What to do |
|---|---|---|
| One sample, two cell-states within it | Per-cell only | Path B (Wilcoxon or MAST) — WITHIN-sample is fine |
| Atlas-scale (200+ donors) but cell-type-imbalanced | Pseudobulk possible but variable | Pseudobulk + mixed model (muscat::mmDS) |
| Cell-state continuum, no discrete groups | No groups to compare | Trajectory DE (tradeSeq, scLANE) — different recipe |
| <3 donors per condition | Pseudobulk possible but underpowered | Description only; flag insufficient replication |
| Differential abundance, not expression | Counts of cell-types per sample | scProportionTest / propellor (different recipe) |

The thing that is ALMOST NEVER right: running a bulk DE tool against a
per-cell counts matrix and reporting the result as a condition effect.
That's what this entire reference exists to prevent.

## Sources

- Squair et al. 2021 — Confronting false discoveries in single-cell DE, Nat Commun (DOI 10.1038/s41467-021-25960-2). **The canonical reference.**
- Soneson & Robinson 2018 — Bias, robustness and scalability in single-cell DE analysis, Nat Methods (DOI 10.1038/nmeth.4612)
- Crowell et al. 2020 — muscat (mixed-model DS analysis), Nat Commun (DOI 10.1038/s41467-020-19894-4) — for the muscat::mmDS alternative
- Murphy & Skene 2022 — A balanced measure shows superior performance of pseudobulk methods, Nat Commun (DOI 10.1038/s41467-022-32604-6) — confirms the Squair finding with a different metric
- `Seurat:::DESeq2DETest`, `Seurat:::MASTDETest` (Seurat 5.5.0 source) — for what the wrappers actually build
- Seurat v5 DE vignette — satijalab.org/seurat/articles/de_vignette
