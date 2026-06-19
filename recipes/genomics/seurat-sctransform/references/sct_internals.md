# SCTransform internals — the model, the residuals, the slots

What SCTransform actually does (per-gene regularized negative-binomial fit,
Pearson residuals), how the SCT assay differs from RNA, where the
per-gene residual variance table lives (and the v5 gotcha that `meta.features`
is empty), and the `vst.flavor` versions.

Load this reference when:
- Inspecting residual-variance ranks (the recipe's Step 4 visualization).
- A pretrained agent guess about the SCT assay's structure was wrong
  (looking for `@meta.features` and finding it empty).
- You need to explain "Why SCT?" vs LogNormalize to the user.
- You need to choose `vst.flavor` between `"v2"` and `"v1"`.
- Considering `vars.to.regress` covariates.

## The model — regularized negative binomial per gene

SCTransform fits, per gene g:

    y_g ~ NegativeBinomial(μ_g(c), θ_g)
    log μ_g(c) = β_0,g + β_1,g · log(s_c) + (Xβ)_g

where `s_c` is the total UMI count of cell c (the sequencing depth covariate),
and `Xβ` is the linear contribution of any cell-level covariates passed via
`vars.to.regress` (`percent.mt`, cell-cycle scores, batch, …). Each gene
gets its own intercept (`β_0`), depth slope (`β_1`), and overdispersion
(`θ`).

**The regularization step** (the "regularized" in regularized NB): the
per-gene θ estimates are very noisy at low expression. SCTransform smooths
them by kernel regression on the mean-expression dimension, then refits each
gene with its smoothed θ. This is the load-bearing trick from Hafemeister &
Satija 2019 — the per-gene fit is noisy but the smoothed surface is reliable.

**Pearson residuals** — the output of the model is, per (cell, gene),

    r_{c,g} = (y_{c,g} - μ_g(c)) / sqrt(μ_g(c) + μ_g(c)^2 / θ_g)

Cells with more counts of g than the model predicts have positive residuals;
fewer counts than predicted, negative. These residuals are stored in
`scale.data` of the SCT assay and are what PCA reads.

The "v2" version (Lause et al. 2021) tightens the regularization with
several technical improvements over v1; the math is mostly the same.

## Why SCT vs LogNormalize

`NormalizeData(LogNormalize)` does:

    x_{c,g} = log1p(y_{c,g} / s_c * 10000)

— per-cell scaling by total counts (size-factor normalization) followed by
log1p. Two problems:

1. **Variance after log normalization is depth-dependent.** Genes with
   higher mean expression have higher variance even after log. Lause et al.
   2021 showed this distorts HVG selection — high-expression genes are
   preferentially called variable.
2. **Heterogeneous depth between cell types is uncorrected.** If T cells
   sequence at 5k UMI/cell and plasma cells at 50k, the same biological
   gene-expression difference produces different normalized values per
   cell type.

SCTransform corrects both: the per-gene NB fit predicts the expected count
GIVEN cell depth, so the residual is depth-corrected by construction.
Pearson residuals have ~unit variance for non-DE genes.

**Practical consequence:** SCT recovers finer cluster structure than
LogNormalize on heterogeneous-depth or low-depth data. The PBMC3k vignette
shows ~12–14 SCT clusters at `resolution = 0.8` vs ~8 LogNormalize clusters
at the same resolution.

**When SCT does NOT help:** clean 10x v3 PBMC at typical depth (5–15k UMI/
cell) is the regime LogNormalize was tuned on. SCT is identical-to-slightly-
better on those. Reach for SCT when:

- Median UMI per cell well below 5k.
- Heterogeneous depth across cell types (e.g. T vs plasma vs erythrocyte
  contamination).
- The user explicitly names SCTransform / Pearson residuals / regularized
  NB.

## Source-verified `SCTransform` signature (Seurat 5.5.0)

```r
SCTransform(object,
            assay = "RNA",
            new.assay.name = "SCT",
            reference.SCT.model = NULL,
            do.correct.umi = TRUE,
            ncells = 5000,
            residual.features = NULL,
            variable.features.n = 3000,        # NOTE: 3000, not LogNormalize's 2000
            variable.features.rv.th = 1.3,     # residual-variance threshold (overridden if .n given)
            vars.to.regress = NULL,
            do.scale = FALSE,
            do.center = TRUE,
            clip.range = c(-sqrt(ncol(object)/30), sqrt(ncol(object)/30)),
            vst.flavor = "v2",                  # NOTE: v2 is the default in v5
            conserve.memory = FALSE,
            return.only.var.genes = TRUE,
            seed.use = 1448145,
            verbose = TRUE,
            ...)
```

Argument notes (from source + experience):

- **`assay = "RNA"`** — reads from the RNA assay. The new SCT data goes
  into `new.assay.name = "SCT"` and the default assay flips to SCT.
- **`do.correct.umi = TRUE`** — produces a corrected UMI count matrix in
  the SCT `counts` layer. Used by `data` (`log1p` of these). Set FALSE to
  keep raw counts; rarely useful.
- **`ncells = 5000`** — random subsample of cells used to fit the
  regularized NB model. For >5000 cells the model fits a subsample then
  scores all cells; for ≤5000, it uses all. Increase for very heterogeneous
  data; the default is fine for most PBMC-scale samples.
- **`variable.features.n = 3000`** — the SCT default returns 3000 variable
  features (vs LogNormalize's `nfeatures = 2000`). The "variable" set is
  ranked by residual variance. PCA reads only these via `scale.data`.
- **`variable.features.rv.th = 1.3`** — alternative selection by residual
  variance threshold. `variable.features.n` takes priority if both given.
- **`vars.to.regress`** — covariates to remove during the model fit. See
  next section.
- **`do.scale = FALSE`, `do.center = TRUE`** — the Pearson residuals are
  approximately unit-variance by construction, so `do.scale` is unneeded.
  `do.center = TRUE` shifts each gene to mean 0.
- **`clip.range`** — residuals are clipped to ±sqrt(N/30) before storing.
  Prevents extreme cells (sparse high-count genes in rare cells) from
  dominating PCA. The default usually works; ignore unless residuals show
  pathological tails.
- **`vst.flavor = "v2"`** — see below.
- **`return.only.var.genes = TRUE`** — `scale.data` is restricted to the
  variable features. Set FALSE if you need residuals for a non-variable
  gene later (FeaturePlot, marker validation). Doubles memory for large
  objects.

## `vst.flavor` — v1 vs v2

| Value | Method | When |
|---|---|---|
| `"v2"` (default in Seurat v5) | Lause et al. 2021 — corrected per-gene NB with theta regularization via `glmGamPoi` if available, else built-in optimizer | Default. Recommended for all new analyses |
| `"v1"` | Hafemeister & Satija 2019 — original implementation, slower | ONLY for reproducing a pre-v5 analysis that depended on v1 results |

The v1-vs-v2 numerical differences are small but real — HVG ranks shift
~5–15%, PCA loadings shift accordingly, cluster boundaries jiggle. If you're
extending a published analysis that used v1, pass `vst.flavor = "v1"`
explicitly.

`glmGamPoi` (Bioconductor) accelerates the per-gene NB fit by 5–20×.
SCTransform auto-detects it; the recipe's Install block ensures it's present.

## `vars.to.regress` — covariates the model removes

Pass a character vector of metadata columns to remove their linear
contribution during normalization. The vignette uses `vars.to.regress =
"percent.mt"`.

| Covariate | Why regress | When NOT |
|---|---|---|
| `"percent.mt"` | Dying-cell artifact bleeds into the scaled data and dominates a low PC | If MT% is biologically relevant (mitochondrial disease study, energy metabolism comparison) |
| `"S.Score"`, `"G2M.Score"` | Cell cycle dominates a PC and obscures lineage structure | If cell-cycle dynamics ARE the biology of interest (cancer proliferation, stem-cell differentiation) |
| `"nCount_RNA"` | Depth is already in the model — DON'T add it | Always; this is the variable SCT already covers |
| `"orig.ident"` / batch | Removes inter-sample batch effects within a single SCTransform fit | Use proper integration (`seurat-integration` recipe) for multi-sample; regressing here is the wrong tool |
| `"percent.ribo"` | Rarely the dominant confound; defer | Most analyses |

Caveats:
- **Every additional covariate shrinks degrees of freedom** in the per-gene
  NB fit. The model has ~3 parameters per gene already (intercept, depth
  slope, theta); adding 3 more covariates means 6 per gene, and small
  populations don't have enough cells to estimate that many parameters
  cleanly.
- **Add covariates only when there's EVIDENCE** they're a confound: the PC
  loadings (Step 5) show a PC dominated by `MT-*` genes, ribosomal genes,
  or cell-cycle genes. Pre-emptive regression of nothing-much wastes
  degrees of freedom.

## SCT assay structure — where things live

Source-verified slotnames (Seurat 5.5.0):

```
SCTAssay slots:
  SCTModel.list, counts, data, scale.data, assay.orig, var.features,
  meta.features, misc, key

SCTModel (each element of SCTModel.list) slots:
  feature.attributes, cell.attributes, clips, umi.assay, model, arguments,
  median_umi
```

| Slot | Contents |
|---|---|
| `counts` | Corrected UMI counts (cells × genes, sparse). Used as the basis for `data` |
| `data` | `log1p` of corrected counts. Used by `FeaturePlot`, `DotPlot` (default), DE under `slot = "data"` |
| `scale.data` | Pearson residuals (cells × variable_features, dense). Used by `RunPCA`, `DimHeatmap`. Restricted to variable features if `return.only.var.genes = TRUE` |
| `var.features` | Character vector of selected variable feature names |
| `meta.features` | **EMPTY** for a single-sample SCT object in Seurat v5. Used to hold per-gene metadata; the per-gene residual variance lives in `SCTModel.list[[i]]@feature.attributes`, NOT here |
| `SCTModel.list` | One `SCTModel` object per sample. Single-sample: length 1. Merged: length N |

### The `meta.features` gotcha — KEY for residual-variance retrieval

Pretrained intuition says "the SCT assay's variable-feature ranking lives
in `obj[["SCT"]]@meta.features$residual_variance`". This is WRONG for
Seurat v5 — `@meta.features` is an empty data frame for a single-sample
SCT object. The actual table lives one level deeper:

```r
# CORRECT — single-sample SCT
sct_var <- slot(obj[["SCT"]]@SCTModel.list[[1]], "feature.attributes")
# Or equivalently:
sct_var <- obj[["SCT"]]@SCTModel.list[[1]]@feature.attributes

# WRONG — returns empty data frame
sct_var <- obj[["SCT"]]@meta.features
```

`feature.attributes` carries columns: `mean`, `detection_rate`, `gmean`,
`variance`, `residual_mean`, `residual_variance`, `theta`, `(Intercept)`,
`log_umi` (the depth slope coefficient), and any `vars.to.regress`
coefficient. The recipe uses `residual_variance` for the rank plot.

For a MERGED multi-sample SCT object, every element of `SCTModel.list`
carries its own `feature.attributes`. The recipe handles single-sample;
multi-sample workflows (integration / merging) should aggregate
appropriately and are out of scope here — see the integration recipe.

This bug was caught live on the pilot recipe validation (`test-seurat-
sctransform/TEST_LOG.md` Round 1 patch). Without the fix, `top20` prints
`character(0)` and `sct_residuals.png` is a blank plot titled "top 0 HVGs
labeled".

## Default assay flip — implications

After `SCTransform()`, `DefaultAssay(obj) == "SCT"`. Every downstream call
that reads `DefaultAssay(obj)` implicitly — `RunPCA`, `FindNeighbors`,
`FindClusters`, `RunUMAP`, `FindAllMarkers`, `DimPlot`, `DotPlot`,
`FeaturePlot` — uses the SCT assay.

**You can still access the RNA assay** for plots / DE:
- `DefaultAssay(obj) <- "RNA"` switches the default.
- `FeaturePlot(obj, features = ..., assay = "RNA")` overrides per call.
- `DotPlot(obj, assay = "RNA", ...)` overrides per call.

This matters because:
- The RNA `data` layer (log-normalized counts) may need re-running
  (`NormalizeData(obj, assay = "RNA")`) if you ran SCT on a fresh object —
  SCT doesn't normalize the RNA assay.
- DE on the RNA assay is sometimes preferred over SCT residuals — see
  `sct_de_gotchas.md`.

**DO NOT call `NormalizeData` / `FindVariableFeatures` / `ScaleData` on the
SCT assay.** SCTransform already produced the equivalents. Running them
again overwrites the residuals and breaks the model.

## Single-sample SCT vs multi-sample integration

This recipe covers SINGLE-sample SCT. For multi-sample integration with
SCT:

1. SCTransform each sample independently.
2. `SelectIntegrationFeatures()` + `PrepSCTIntegration()` align the SCT
   models across samples.
3. `FindIntegrationAnchors(normalization.method = "SCT", ...)` + 
   `IntegrateData(normalization.method = "SCT", ...)` (Seurat v4 API) or
   `IntegrateLayers(method = CCAIntegration, normalization.method = "SCT",
   ...)` (Seurat v5 API).
4. `PrepSCTFindMarkers()` is REQUIRED before DE on the merged object — it
   re-scales the per-model count corrections to a shared median UMI.

The full integration workflow lives in the `seurat-integration` recipe.
