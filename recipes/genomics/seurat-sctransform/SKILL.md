---
name: seurat-sctransform
description: Single-sample scRNA-seq workflow with SCTransform regularized-NB normalization (R/Seurat v5) — replaces NormalizeData + FindVariableFeatures + ScaleData with one SCTransform() call, downstream PCA/UMAP/clustering/markers run on the SCT assay, and PrepSCTFindMarkers gates DE on merged multi-model objects. Preferred when sequencing depth is low or varies a lot across cell types.
when_to_use: ONE single-cell RNA-seq sample (10x CellRanger dir, .h5, or counts matrix) where the user asks for SCTransform / sctransform / regularized negative binomial normalization, OR the data is low-depth or has heterogeneous depth across populations (the regime where log-normalize over-corrects high-count cells). Use this when the session is R and Seurat v5 is in play. For log-normalize defaults, use seurat-scrna-v2. For multi-sample integration use seurat-integration. For a Python-native equivalent see the scanpy / pearson-residuals recipe.
invocation: interactive+batch
requires_tools: [run_r]
capabilities_needed: [Seurat, glmGamPoi]
keywords: [Seurat, Seurat v5, SCTransform, sctransform, regularized negative binomial, Pearson residuals, scRNA-seq, single cell, single-cell, normalization, percent.mt, glmGamPoi, vst.flavor, v2, SCT assay, PrepSCTFindMarkers, FindAllMarkers, RunPCA, RunUMAP, FindClusters, low depth, heterogeneous depth, R]
produces: [qc_violins_pre.png, qc_violins_post.png, sct_residuals.png, pca_elbow.png, umap_clusters.png, markers_dotplot.png, markers_featureplot.png, cluster_markers.csv, seurat_sct_processed.rds, seurat_sct_processed.lstar.zarr]
domain: genomics
source: "Seurat SCTransform vignette (Satija Lab) — https://satijalab.org/seurat/articles/sctransform_vignette — plus SCTransform / PrepSCTFindMarkers v5 reference pages for argument defaults. Method paper: Hafemeister & Satija 2019, Genome Biology (regularized NB); Lause et al. 2021 (v2 flavor)."
---

# Single-sample scRNA-seq with SCTransform (R/Seurat v5)

Same downstream shape as `seurat-scrna-v2` — load, QC, dim-reduce, cluster,
UMAP, find markers — but Step 4 (normalize + HVG + scale) collapses to **one
`SCTransform()` call** that fits a regularized negative-binomial model per
gene and returns Pearson residuals. The vignette's framing: a more effective
normalization that better separates technical from biological variance,
especially when sequencing depth varies a lot across cell populations or the
sample is shallow overall.

**When to pick THIS over `seurat-scrna-v2`:** low-depth datasets (median UMI
per cell well below 5k), heterogeneous-depth datasets (some clusters
sequenced 10× deeper than others — e.g. plasma cells vs. T cells), or any
case where the user explicitly names SCTransform / sctransform / regularized
NB / Pearson residuals.

**Verified against Seurat v5.x.** `SCTransform()` defaults to `vst.flavor =
"v2"` in v5, which uses `glmGamPoi` for the per-gene NB fits and is
substantially faster than v1 — install `glmGamPoi` (Bioconductor) once. The
SCT assay becomes the default downstream; for DE on the SCT assay,
`PrepSCTFindMarkers()` must be called first.

## Bundled references — load on demand

This SKILL.md is self-contained for the standard workflow. Load these only
when the task crosses into a topic the body summarizes but doesn't unpack:

- `references/sct_internals.md` — the regularized-NB model math, what
  Pearson residuals are, the `vst.flavor` v1 vs v2 difference, full
  `SCTransform()` arg semantics (`ncells`, `variable.features.n`,
  `clip.range`, `return.only.var.genes`), the SCT assay slot structure,
  and the **key v5 gotcha** that `obj[["SCT"]]@meta.features` is empty
  for single-sample objects (the residual variance lives in
  `SCTModel.list[[1]]@feature.attributes`).
- `references/sct_de_gotchas.md` — `PrepSCTFindMarkers` semantics
  (mandatory, even single-sample), SCT vs RNA assay DE tradeoffs, why
  `slot = "data"` (not `"scale.data"`) for Wilcoxon, what changes when
  objects are merged, and pseudobulk DE pointers.
- `references/installation_and_io.md` — reader functions (`Read10X` /
  `Read10X_h5` / `ReadMtx`), 10x v2 vs v3 detection, gz handling,
  organism-specific MT prefixes, dependency table including
  `glmGamPoi`/`sctransform`.
- `references/qc_and_thresholds.md` — QC metric definitions,
  tissue-specific cutoff conventions, the `subset()` filter, doublet
  detection, cell-cycle scoring.
- `references/figure_style.md` — Seurat-collection figure conventions
  (palette, theme ordering, alpha-poke pattern, dpi/bg/width
  conventions).

## Install

Idempotent — re-running is a no-op.

```r
suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(patchwork)
  library(tidyr)
})

# glmGamPoi accelerates vst.flavor="v2"; SCTransform auto-detects it.
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("glmGamPoi", update = FALSE, ask = FALSE)
}

# Optional Wilcoxon accelerator for FindAllMarkers — auto-detected by Seurat v5.
if (!requireNamespace("presto", quietly = TRUE)) {
  # devtools::install_github("immunogenomics/presto")
}

stopifnot(packageVersion("Seurat") >= "5.0.0")  # this recipe is v5-only
```

`library(Seurat)` does NOT attach `ggplot2`, `dplyr`, or `cowplot` — load
them by name. The `sctransform` package is the underlying engine; Seurat's
`SCTransform()` wraps it.

For the full dependency table (what each package is for + alternatives) see
`references/installation_and_io.md`.

## Decisions to surface up front

Tell the user these are the analysis-defining choices:

1. **Species / MT prefix** — `^MT-` (human, uppercase), `^mt-` (mouse),
   organism-specific for others. `percent.mt` is computed BEFORE
   `SCTransform()` because the recipe regresses it out in the model.
   Per-organism table in `references/installation_and_io.md`.
2. **QC thresholds** — `nFeature_RNA` floor / ceiling and `percent.mt`
   ceiling. SCT is less sensitive to depth than log-normalize, but
   empty-droplet / dying-cell removal still matters. Read thresholds off
   the per-metric quantile tables in Step 2; tissue-specific defaults in
   `references/qc_and_thresholds.md`.
3. **`vars.to.regress`** — covariates the SCT model removes during
   normalization. The vignette regresses `"percent.mt"`. Add `"S.Score"` /
   `"G2M.Score"` if a Step-5 PC is dominated by cell cycle. Don't pile
   unmotivated covariates here — every additional one shrinks degrees of
   freedom in the per-gene fit. Full guidance in `references/sct_internals.md`.
4. **`vst.flavor`** — `"v2"` (Seurat v5 default; uses `glmGamPoi`,
   recommended). Drop to `"v1"` ONLY to reproduce a pre-v5 analysis. The
   v1 vs v2 numerical differences are documented in
   `references/sct_internals.md`.
5. **Number of PCs (`dims`)** — the vignette uses `1:30`; SCT's improved
   normalization tends to spread structure across more PCs than
   log-normalize. Default `DIMS_CHOSEN = 30`.
6. **Clustering resolution** — `FindClusters(resolution = ...)`; default
   `0.8` per the vignette (slightly higher than log-normalize's `0.5`
   because SCT recovers finer structure).
7. **Mode** — INTERACTIVE (default) when this is the primary analysis;
   BATCH when invoked as one of many sub-runs over multiple samples.

Figures the user will see as the analysis proceeds:
`qc_violins_pre.png`, `qc_violins_post.png`, `sct_residuals.png`,
`pca_elbow.png`, `umap_clusters.png`, `markers_dotplot.png`,
`markers_featureplot.png`.

---

## Step 1 — Load the counts matrix and create the Seurat object

Same three input shapes as `seurat-scrna-v2`; pick by what's on disk.

```r
# 10x CellRanger directory.
counts <- Read10X(data.dir = "/path/to/filtered_feature_bc_matrix")
```

```r
# 10x .h5
counts <- Read10X_h5("/path/to/filtered_feature_bc_matrix.h5")
```

```r
# Loose, GSM-prefixed GEO triplets (Read10X cannot find them — non-standard
# names). feature.column = 2 → gene symbols (column 1 is Ensembl).
counts <- ReadMtx(
  mtx      = "/path/to/GSMxxxxxxx_matrix.mtx.gz",
  cells    = "/path/to/GSMxxxxxxx_barcodes.tsv.gz",
  features = "/path/to/GSMxxxxxxx_features.tsv.gz",
  feature.column = 2
)
```

The vignette uses **no** `min.cells` / `min.features` on
`CreateSeuratObject` (those args default to 0 — see `?CreateSeuratObject`),
because SCTransform's variable.features.rv.th + variance filter handle gene
selection downstream. We keep a light empty-droplet filter here for the
same reason `seurat-scrna-v2` does — to prune obviously-empty barcodes
before QC plotting.

```r
obj <- CreateSeuratObject(
  counts       = counts,
  project      = "sample",
  min.cells    = 3,         # drop genes seen in <3 cells (sparsity, not biology)
  min.features = 200        # drop barcodes with <200 genes (empty droplets)
)

n_genes_raw <- nrow(counts); n_cells_raw <- ncol(counts)
cat(sprintf("Pre-filter: %d genes x %d cells\n", n_genes_raw, n_cells_raw))
cat(sprintf("After CreateSeuratObject (min.cells=3, min.features=200): %d genes x %d cells\n",
            nrow(obj), ncol(obj)))

# Sanity-check the loaded object — symbol vs Ensembl, MT prefix coverage.
stopifnot(inherits(obj, "Seurat"))
head(rownames(obj))
n_mt <- sum(grepl("^MT-", rownames(obj)))
cat(sprintf("MT genes matched by '^MT-': %d  (0 means wrong prefix; mouse uses '^mt-')\n", n_mt))
```

**Report:** input format, raw cells/genes, post-load cells/genes, whether MT
prefix matched. If `n_mt == 0`, stop and fix the prefix.

For reader signatures, gz auto-detection, AnnData / h5Seurat bridges, and
per-organism MT prefix table, read `references/installation_and_io.md`.

---

## Step 2 — Compute QC metrics and pick thresholds

Run BEFORE `SCTransform()` so `percent.mt` exists for the model's
`vars.to.regress`, and so QC plots reflect the raw input (not residuals).

```r
obj <- PercentageFeatureSet(obj, pattern = "^MT-",   col.name = "percent.mt")
obj <- PercentageFeatureSet(obj, pattern = "^RP[SL]", col.name = "percent.ribo")

# Quantile tables — pick floors from qs_lo, ceilings from qs_hi.
qs_hi <- c(0.50, 0.75, 0.90, 0.95, 0.975, 0.99)
qs_lo <- c(0.005, 0.01, 0.025, 0.05, 0.10)

qtab_hi <- data.frame(
  quantile     = qs_hi,
  nFeature_RNA = quantile(obj$nFeature_RNA, qs_hi),
  nCount_RNA   = quantile(obj$nCount_RNA,   qs_hi),
  percent.mt   = quantile(obj$percent.mt,   qs_hi),
  percent.ribo = quantile(obj$percent.ribo, qs_hi)
)
print(qtab_hi, row.names = FALSE)

qtab_lo <- data.frame(
  quantile     = qs_lo,
  nFeature_RNA = quantile(obj$nFeature_RNA, qs_lo),
  nCount_RNA   = quantile(obj$nCount_RNA,   qs_lo)
)
print(qtab_lo, row.names = FALSE)

THRESH <- list(
  nFeature_low  = 200,    # REPLACE — read off the quantile table
  nFeature_high = 5000,   # REPLACE
  mt_high       = 15      # REPLACE
)
```

SCT is more robust to depth heterogeneity than log-normalize, but it
**doesn't rescue dying cells or empty droplets** — apply the same QC
thresholds you would otherwise.

**Report:** thresholds chosen and the corresponding quantile.

For metric definitions, tissue-specific cutoffs, and doublet detection,
read `references/qc_and_thresholds.md`.

---

## Step 3 — Apply QC: pre-filter plot, filter, post-filter plot

### 3a. Pre-filter violin coloured by qc_kept

```r
obj$qc_kept <- with(obj@meta.data,
  ifelse(nFeature_RNA > THRESH$nFeature_low  &
         nFeature_RNA < THRESH$nFeature_high &
         percent.mt   < THRESH$mt_high,
         "kept", "filtered"))

qc_long <- as.data.frame(obj@meta.data[, c("nFeature_RNA","nCount_RNA",
                                           "percent.mt","percent.ribo","qc_kept")])
qc_long$cell <- rownames(qc_long)
qc_long <- pivot_longer(qc_long, c(-cell, -qc_kept),
                        names_to = "metric", values_to = "value")
qc_long$metric <- factor(qc_long$metric,
                         levels = c("percent.ribo","percent.mt",
                                    "nCount_RNA","nFeature_RNA"))

thresholds <- data.frame(
  metric = factor(c("percent.ribo","percent.mt","nCount_RNA",
                    "nFeature_RNA","nFeature_RNA"),
                  levels = levels(qc_long$metric)),
  value  = c(NA_real_, THRESH$mt_high, NA_real_,
             THRESH$nFeature_low, THRESH$nFeature_high)
)

p_vln <- ggplot(qc_long, aes(x = "", y = value)) +
  geom_violin(aes(fill = metric),
              width = 0.85, colour = "black", linewidth = 0.4,
              scale = "width", trim = FALSE, alpha = 0.85) +
  geom_jitter(aes(colour = qc_kept),
              width = 0.32, height = 0, size = 0.20, alpha = 0.10) +
  geom_hline(data = na.omit(thresholds), aes(yintercept = value),
             colour = "red", linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ metric, scales = "free_x", ncol = 1, strip.position = "left") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_colour_manual(values = c(kept = "black", filtered = "red"),
                      guide = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
  labs(x = NULL, y = NULL, colour = NULL,
       title = sprintf("QC metrics, pre-filter (n = %d cells)", ncol(obj))) +
  theme_cowplot() +
  theme(strip.placement   = "outside",
        strip.background  = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
        axis.text.y       = element_blank(),
        axis.ticks.y      = element_blank(),
        axis.line.y       = element_blank(),
        panel.spacing.y   = unit(2, "pt"),
        plot.title        = element_text(size = 12, face = "bold"),
        legend.position   = "bottom")

ggsave("qc_violins_pre.png", plot = p_vln,
       bg = "white", dpi = 120, width = 7, height = 4.7)
```

### 3b. Apply the filter

```r
n_before <- ncol(obj)
obj <- subset(obj, subset =
  nFeature_RNA > THRESH$nFeature_low  &
  nFeature_RNA < THRESH$nFeature_high &
  percent.mt   < THRESH$mt_high)
n_after <- ncol(obj)
cat(sprintf("Cells: %d -> %d  (removed %d, %.1f%%)\n",
            n_before, n_after, n_before - n_after,
            100 * (n_before - n_after) / n_before))
```

### 3c. Post-filter violin

```r
qc_long_post <- as.data.frame(obj@meta.data[, c("nFeature_RNA","nCount_RNA",
                                                "percent.mt","percent.ribo")])
qc_long_post$cell <- rownames(qc_long_post)
qc_long_post <- pivot_longer(qc_long_post, -cell,
                             names_to = "metric", values_to = "value")
qc_long_post$metric <- factor(qc_long_post$metric,
                              levels = c("percent.ribo","percent.mt",
                                         "nCount_RNA","nFeature_RNA"))

p_vln_post <- ggplot(qc_long_post, aes(x = "", y = value, fill = metric)) +
  geom_violin(width = 0.85, colour = "black", linewidth = 0.4,
              scale = "width", trim = FALSE, alpha = 0.85) +
  geom_jitter(width = 0.32, height = 0, size = 0.20,
              alpha = 0.10, colour = "grey15") +
  facet_wrap(~ metric, scales = "free_x", ncol = 1, strip.position = "left") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = NULL,
       title = sprintf("QC metrics, post-filter (n = %d cells)", ncol(obj))) +
  theme_cowplot() +
  theme(strip.placement   = "outside",
        strip.background  = element_blank(),
        strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
        axis.text.y       = element_blank(),
        axis.ticks.y      = element_blank(),
        axis.line.y       = element_blank(),
        panel.spacing.y   = unit(2, "pt"),
        plot.title        = element_text(size = 12, face = "bold"))

ggsave("qc_violins_post.png", plot = p_vln_post,
       bg = "white", dpi = 120, width = 7, height = 4.5)
```

**Report:** cells before/after, percent removed, dominant failure mode, any
truncation visible at the cutoff in the post-filter violin.

For the figure conventions (palette, Set2 fills, theme ordering) read
`references/figure_style.md`.

---

## Step 4 — SCTransform (replaces NormalizeData + HVG + ScaleData)

ONE call. Defaults verified against `?SCTransform` for Seurat v5: assay
`"RNA"` → new assay `"SCT"`, `variable.features.n = 3000`, `vst.flavor =
"v2"` (uses `glmGamPoi`), `do.correct.umi = TRUE`, `seed.use = 1448145`.

```r
# vars.to.regress = "percent.mt" matches the vignette. Add S.Score/G2M.Score
# only if Step 5 PC loadings show cell-cycle domination.
obj <- SCTransform(obj,
                   vars.to.regress = "percent.mt",
                   vst.flavor = "v2",            # uses glmGamPoi automatically
                   verbose = FALSE)

# After SCTransform, the default assay is "SCT". All downstream calls
# (RunPCA, FindNeighbors, FindClusters, RunUMAP) read it implicitly.
cat("Default assay:", DefaultAssay(obj), "\n")
cat("SCT variable features:", length(VariableFeatures(obj, assay = "SCT")), "\n")
```

The SCT assay stores three layers: `counts` (corrected UMI counts),
`data` (log1p of corrected counts), and `scale.data` (Pearson residuals
for the variable features). PCA reads `scale.data`; FeaturePlot reads
`data`.

### Visualize the residuals — sanity-check the model fit

If residuals look ratty (extreme tails, MT/ribo dominating), the model
likely needs additional `vars.to.regress` or a tighter QC cutoff.

**The v5 gotcha (caught on the pilot validation):** the per-gene residual
variance table is NOT in `obj[["SCT"]]@meta.features` (that slot is empty
for a single-sample SCT object in Seurat v5). It lives in
`obj[["SCT"]]@SCTModel.list[[1]]@feature.attributes`. See
`references/sct_internals.md` for the slot-structure map.

```r
# ggrepel is a Seurat dependency; install BEFORE we build the plot that uses it.
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")

# Single-sample SCT — one SCTModel; per-gene residual variance lives here.
sct_var <- slot(obj[["SCT"]]@SCTModel.list[[1]], "feature.attributes")
sct_var <- sct_var[order(-sct_var$residual_variance), ]
# Restrict the labeled set to the SCT-selected variable features so the
# rank axis matches what downstream PCA actually sees.
hvg_sct <- VariableFeatures(obj, assay = "SCT")
sct_var_hvg <- sct_var[rownames(sct_var) %in% hvg_sct, ]
top20 <- head(rownames(sct_var_hvg), 20)
print(top20)

# Rank vs residual variance scatter, top 30 labeled.
df_rv <- data.frame(rank = seq_len(nrow(sct_var_hvg)),
                    residual_variance = sct_var_hvg$residual_variance,
                    gene = rownames(sct_var_hvg))
df_rv$label <- ifelse(df_rv$rank <= 30, df_rv$gene, NA_character_)

p_rv <- ggplot(df_rv, aes(rank, residual_variance)) +
  geom_point(alpha = 0.35, size = 0.9, colour = "grey40") +
  ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 30,
                           na.rm = TRUE) +
  scale_x_log10() +
  labs(title = sprintf("SCT residual variance (top %d HVGs labeled)",
                       min(30, nrow(df_rv))),
       x = "gene rank (log10)", y = "residual variance") +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave("sct_residuals.png", plot = p_rv,
       bg = "white", dpi = 120, width = 8, height = 5.5)

# Flag suspicious dominators.
suspect <- grep("^(MT-|RPS|RPL|MALAT1)", top20, value = TRUE)
cat("QC-warning SCT-variable features:",
    if (length(suspect)) paste(suspect, collapse = ", ") else "none", "\n")
```

**Report:** SCT variable-feature count, top 20 by residual variance, any
suspicious dominators (suggests tightening QC or adding `vars.to.regress`).

> **Pitfall — DON'T call `NormalizeData` / `FindVariableFeatures` /
> `ScaleData` on the SCT assay.** SCTransform already produced the
> equivalents. Running them again on the `"SCT"` assay overwrites the
> residuals and breaks the model. (`NormalizeData` on the RNA assay is
> OK — and is required ONLY if you want to use the RNA assay for DE in
> Step 7, see the pitfall there.)

For the regularized-NB model math, `vars.to.regress` covariates and when
to add each, the full SCTransform arg semantics, the SCT slot structure,
and the multi-sample integration pointer, read
`references/sct_internals.md`.

---

## Step 5 — PCA on the SCT assay, choose `dims`

PCA reads `DefaultAssay(obj) = "SCT"` implicitly. The vignette uses
`dims = 1:30` for downstream; the elbow below confirms.

```r
obj <- RunPCA(obj, verbose = FALSE)

# Print top loadings on PCs 1-5 — sanity-check PC1 isn't MT/ribo/cell-cycle.
print(obj[["pca"]], dims = 1:5, nfeatures = 8)
```

Elbow plot — per-PC variance (solid blue) and cumulative (grey dashed)
on one axis, normalized to the SCT scale.data total variance.

```r
DIMS_CHOSEN <- 30   # vignette default; override if elbow shows a clear cliff

sd_mat    <- GetAssayData(obj, assay = "SCT", layer = "scale.data")
total_var <- sum(apply(sd_mat, 1, var))
var_pc    <- Stdev(obj, reduction = "pca")^2
prop      <- var_pc / total_var
cum       <- cumsum(prop)
df_elb    <- data.frame(PC = seq_along(prop), prop = prop, cum = cum)
df_long   <- pivot_longer(df_elb, c(prop, cum),
                          names_to = "kind", values_to = "value")
df_long$kind <- factor(df_long$kind, levels = c("prop","cum"),
                       labels = c("per-PC","cumulative"))

p_elbow <- ggplot(df_long, aes(x = PC, y = value,
                               colour = kind, linetype = kind)) +
  geom_vline(xintercept = DIMS_CHOSEN, colour = "red",
             linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = DIMS_CHOSEN - 0.7, y = max(df_long$value) * 0.95,
           label = sprintf("dims = 1:%d (chosen)", DIMS_CHOSEN),
           colour = "red", hjust = 1, size = 3.4) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 0.4, 0.05),
                     limits = c(0, max(df_long$value) * 1.05)) +
  scale_x_continuous(breaks = c(1, seq(5, 50, 5))) +
  scale_colour_manual(values = c("per-PC" = "#1f77b4",
                                 "cumulative" = "grey40")) +
  scale_linetype_manual(values = c("per-PC" = "solid",
                                   "cumulative" = "dashed")) +
  labs(title = "PCA variance explained (SCT assay)",
       x = "principal component",
       y = "proportion of total variance",
       colour = NULL, linetype = NULL) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = c(0.98, 0.55),
        legend.justification = c(1, 0.5),
        legend.background = element_rect(fill = alpha("white", 0.5), colour = NA),
        legend.key = element_blank(),
        legend.title = element_blank())

ggsave("pca_elbow.png", plot = p_elbow,
       bg = "white", dpi = 120, width = 7, height = 4.8)
```

**Report:** `DIMS_CHOSEN`, top-loading genes on PCs 1–5, any PC that looks
technical.

> **SCT spreads structure across more PCs than log-normalize.** A
> log-normalize elbow that suggested `dims = 1:10` will typically extend
> to `1:30` under SCT — using too few PCs underclusters.

---

## Step 6 — Neighbor graph, clustering, UMAP (all on SCT)

```r
obj <- FindNeighbors(obj, dims = 1:DIMS_CHOSEN, verbose = FALSE)

# Vignette uses default resolution (0.8 under FindClusters); start there
# under SCT (SCT recovers finer structure than log-normalize at the same
# resolution).
obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)

obj <- RunUMAP(obj, dims = 1:DIMS_CHOSEN, verbose = FALSE)

cat("Clusters:", length(levels(Idents(obj))), "\n")
print(table(Idents(obj)))
```

UMAP coloured by cluster — alpha-poke `GeomPoint` for readable density
(DimPlot doesn't expose `alpha`).

```r
p_umap <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE,
                  pt.size = 0.4) +
  ggtitle(sprintf("SCT . Louvain res=0.8 . dims=1:%d . n=%d cells . %d clusters",
                  DIMS_CHOSEN, ncol(obj), length(levels(Idents(obj))))) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold")) +  # shrink so the full title fits at width=7.5in
  NoLegend() +              # AFTER theme_cowplot — cowplot 1.2.0 overrides legend.position
  coord_fixed()

for (k in seq_along(p_umap$layers)) {
  if (inherits(p_umap$layers[[k]]$geom, "GeomPoint")) {
    p_umap$layers[[k]]$aes_params$alpha <- 0.6
  }
}

ggsave("umap_clusters.png", plot = p_umap,
       bg = "white", dpi = 120, width = 7.5, height = 6.5)
```

**Report:** algorithm + resolution used, cluster count, cluster sizes, any
cluster <1% of cells.

For the algorithm enum (`1`/`2`/`3`/`4` = Louvain/multilevel/SLM/Leiden),
resolution semantics + cluster-count mapping, and the v4 → v5 UMAP default
changes, the same `references/clustering_choices.md` in the `seurat-scrna-v2`
sibling applies — read it there.

---

## Step 7 — Find cluster markers (PrepSCTFindMarkers first)

For DE on the SCT assay, **`PrepSCTFindMarkers()` must run before
`FindAllMarkers` / `FindMarkers`.** It re-corrects counts using the minimum
median UMI across SCT models (no-op for a single-sample object with one
model, but documented as required for general SCT use and harmless here).

```r
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)

# FindAllMarkers reads the SCT assay because DefaultAssay(obj) == "SCT".
# only.pos / min.pct / logfc.threshold match the conventional Wilcoxon
# marker discovery settings.
markers <- FindAllMarkers(obj,
                          assay           = "SCT",
                          only.pos        = TRUE,
                          min.pct         = 0.25,
                          logfc.threshold = 0.25,
                          verbose         = FALSE)

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

write.csv(markers, file = "cluster_markers.csv", row.names = FALSE)

cat("Total significant markers:", nrow(markers), "\n")
print(top5 %>% select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj), n = Inf)
```

> **Pitfall — RNA-assay DE alternative.** Some workflows prefer running DE
> on the un-regularized log-normalized RNA assay even when clustering used
> SCT. To do that:
>
> ```r
> DefaultAssay(obj) <- "RNA"
> obj <- NormalizeData(obj, verbose = FALSE)          # RNA assay only
> markers_rna <- FindAllMarkers(obj, assay = "RNA",
>                               only.pos = TRUE, min.pct = 0.25,
>                               logfc.threshold = 0.25, verbose = FALSE)
> DefaultAssay(obj) <- "SCT"
> ```
>
> Pick ONE (SCT or RNA) per analysis and report which. Mixing in the same
> figure misleads the reader. Full tradeoffs in
> `references/sct_de_gotchas.md`.

Two figures — dotplot and FeaturePlot — same conventions as
`seurat-scrna-v2`.

```r
# Figure 1: dotplot, top 5 markers per cluster.
genes_to_show <- unique(top5$gene)

p_dot <- DotPlot(obj, features = genes_to_show, cluster.idents = FALSE,
                 assay = "SCT") +
  RotatedAxis() +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0, name = "avg expr") +
  ggtitle(sprintf("Top 5 markers per cluster (SCT, Wilcoxon, only.pos, n=%d genes)",
                  length(genes_to_show))) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10))

suppressMessages(
  ggsave("markers_dotplot.png", plot = p_dot, bg = "white", dpi = 120,
         width = max(12, 0.18 * length(genes_to_show)), height = 6.5)
)
```

```r
# Figure 2: FeaturePlot, 4-8 hand-picked canonical lineage markers.
# REPLACE with markers for YOUR tissue/organism. PBMC placeholder below.
canonical <- c("CD3D","CD8A","MS4A1","GNLY","LYZ","PPBP")
canonical <- canonical[canonical %in% rownames(obj)]

# FeaturePlot reads the SCT "data" layer (log1p of corrected counts) by
# default — the right slot for expression on the embedding.
p_feat <- FeaturePlot(obj, features = canonical, order = TRUE,
                      pt.size = 0.3, ncol = 3) &
  scale_colour_gradient(low = "grey85", high = "#b2182b") &
  theme_cowplot() &
  theme(legend.position = "right", legend.key.size = unit(0.4, "cm"))

for (i in seq_along(p_feat)) {
  pl <- p_feat[[i]]
  if (!is.null(pl$layers)) {
    for (k in seq_along(pl$layers)) {
      if (inherits(pl$layers[[k]]$geom, "GeomPoint")) {
        pl$layers[[k]]$aes_params$alpha <- 0.6
      }
    }
    p_feat[[i]] <- pl
  }
}

n_rows <- ceiling(length(canonical) / 3)
suppressMessages(
  ggsave("markers_featureplot.png", plot = p_feat, bg = "white", dpi = 120,
         width = 13, height = max(4, 4 * n_rows))
)
```

**Report:** total significant markers, top 5 per cluster, any cluster with
zero genes meeting threshold.

For the full `test.use` enum on SCT (Wilcoxon, MAST, ROC, LR; why
`negbinom`/`poisson` need the RNA assay; why `slot = "data"` not
`"scale.data"` for Wilcoxon), the `PrepSCTFindMarkers` mechanism on merged
objects, and the SCT vs RNA DE decision, read
`references/sct_de_gotchas.md`.

---

## Step 8 — Save the processed object

```r
saveRDS(obj, file = "seurat_sct_processed.rds")

sz <- file.info("seurat_sct_processed.rds")$size / 1e6
cat(sprintf("Wrote seurat_sct_processed.rds (%.1f MB)\n", sz))
obj_check <- readRDS("seurat_sct_processed.rds")
stopifnot(identical(dim(obj_check), dim(obj)))
rm(obj_check); invisible(gc())
```

The `.rds` carries both the `RNA` assay (raw + log-normalized if you ran
that path) AND the `SCT` assay (corrected counts + residuals + the SCT
model). The `pca`, `umap`, graphs, and `seurat_clusters` column are
attached. A follow-up session resumes from markers, sub-clustering, or
annotation without re-running SCTransform.

> **SCT models survive `merge()` but require `PrepSCTFindMarkers` after
> merging** — that's the canonical motivating use case for the helper.
> Single-sample SCT pickles cleanly; multi-sample merge-then-DE needs the
> re-prep step before `FindAllMarkers`. See
> `references/sct_de_gotchas.md`.

---

## Batch variant — use INSTEAD of Steps 1–8 when invoked with args="batch"

Branch on `$ARGUMENTS == "batch"` at the top of the body. In batch mode:

- No `ggsave` / `png` calls (skip all per-step figures).
- No quantile tables, no residual scatter, no top-marker print.
- Still save the canonical processed object (`seurat_sct_processed.rds`)
  AND the marker table (`cluster_markers.csv`).
- Print ONE final summary line:
  `"batch ok | <N> cells | <K> clusters | <M> markers | <MB> MB"`.

```r
# args="batch": single consolidated SCT pipeline, no figures.
counts <- Read10X(data.dir = "/path/to/filtered_feature_bc_matrix")
obj <- CreateSeuratObject(counts = counts, project = "sample",
                          min.cells = 3, min.features = 200)
obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent.mt")

obj <- subset(obj, subset =
  nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

obj <- SCTransform(obj, vars.to.regress = "percent.mt",
                   vst.flavor = "v2", verbose = FALSE)

DIMS_CHOSEN <- 30
obj <- RunPCA(obj, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:DIMS_CHOSEN, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:DIMS_CHOSEN, verbose = FALSE)

obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)
markers <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25,
                          verbose = FALSE)

saveRDS(obj, "seurat_sct_processed.rds")
write.csv(markers, "cluster_markers.csv", row.names = FALSE)

sz <- file.info("seurat_sct_processed.rds")$size / 1e6
cat(sprintf("batch ok | %d cells | %d clusters | %d markers | %.1f MB\n",
            ncol(obj), length(levels(Idents(obj))), nrow(markers), sz))
```

---

## Final response checklist

Summarize:

- input format, sample ID, MT-prefix sanity
- raw cells/genes and post-load cells/genes
- chosen QC thresholds, percent of cells removed
- `vars.to.regress` passed to SCTransform, `vst.flavor` used
- SCT variable-feature count, top residual-variance genes (subset),
  any suspicious dominators
- `DIMS_CHOSEN` and any PC that looked technical
- clustering algorithm + resolution, cluster count, cluster sizes
- DE assay used (`SCT` after `PrepSCTFindMarkers`, or `RNA` after
  `NormalizeData`), total significant markers, top markers per cluster,
  clusters with zero markers
- figures shown (filenames)
- saved files (`seurat_sct_processed.rds`, `cluster_markers.csv`)
- caveats: SCT model fit failures on very small/sparse populations,
  vst.flavor version pin, single-sample so PrepSCTFindMarkers is a
  no-op (would matter if data are merged)

---

## See also

- `seurat-scrna-v2` — same workflow shape with the classic LogNormalize +
  FindVariableFeatures + ScaleData path. Switch to it for typical
  10x-depth PBMC-style samples, or when downstream tools require the
  log-normalized RNA assay as the primary slot.
- `seurat-integration` — multi-sample integration. SCT integration uses
  `SCTransform()` per sample, then `IntegrateLayers(method =
  CCAIntegration, normalization.method = "SCT", ...)`; the integration
  recipe covers this.
- `scrna-qc-clustering` — scanpy / Python equivalent. For the
  Pearson-residuals analog in scanpy, see `pp.experimental.normalize_pearson_residuals`.

## Offer an interactive view

Write a viewer-optimized store DIRECTLY from the live Seurat object with lstar (pure R,
highest fidelity — do NOT route through `.h5ad`) and **proactively offer to open it** (a
required part of delivering the result):
```r
DefaultAssay(obj) <- "RNA"        # view raw-count RNA, not SCT residuals
d <- lstar::read_seurat(obj)
lstar::lstar_write_viewer(d, "seurat_sct_processed.lstar.zarr")   # precomputes DE / HVGs /
                                                                  # cell-major counts (optimized)
```
Then call `open_viewer(file_path="seurat_sct_processed.lstar.zarr")` and present the returned
link so the user can explore the UMAP, clusters, and markers in pagoda3 — it opens instantly
(pre-optimized, no on-launch conversion). If `open_viewer` returns `ok:false`, relay the
error rather than a dead link. Format / sharing → **`scrna-viewing-and-interchange`**.
