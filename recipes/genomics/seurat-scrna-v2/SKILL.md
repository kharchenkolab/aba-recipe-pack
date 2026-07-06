---
name: seurat-scrna-v2
description: Single-sample scRNA-seq QC + clustering + markers with R/Seurat v5 — Read10X/Read10X_h5 → CreateSeuratObject → percent.mt + nFeature/nCount QC subset → NormalizeData/FindVariableFeatures/ScaleData → RunPCA + ElbowPlot → FindNeighbors/FindClusters (Louvain) → RunUMAP → FindAllMarkers (Wilcoxon). Default canonical Seurat path; SCTransform sibling recipe handles low-depth / regularized-NB normalization.
when_to_use: ONE single-cell RNA-seq sample (10x CellRanger output dir, .h5 file, or an in-memory counts matrix) for ANY tissue/organism, and the user wants the standard Seurat workflow — QC, normalization, PCA, clustering, UMAP, marker discovery — culminating in a labeled UMAP + a per-cluster marker table. Use this when the session is R-based or the user names Seurat. For SCTransform-based normalization instead of log-normalize, use seurat-sctransform. For multi-sample integration use seurat-integration. For a Python session prefer the scanpy single-sample recipe.
invocation: interactive+batch
requires_tools: [run_r]
capabilities_needed: [Seurat]
keywords: [Seurat, Seurat v5, scRNA-seq, single cell, single-cell, QC, percent.mt, NormalizeData, FindVariableFeatures, ScaleData, RunPCA, ElbowPlot, FindNeighbors, FindClusters, Louvain, RunUMAP, FindAllMarkers, Wilcoxon, marker genes, DimPlot, FeaturePlot, DotPlot, PBMC3k, R]
produces: [qc_violins_pre.png, qc_scatters_pre.png, qc_violins_post.png, hvg_plot.png, pca_elbow.png, pca_heatmap.png, umap_clusters.png, markers_dotplot.png, markers_featureplot.png, cluster_markers.csv, seurat_processed.rds, seurat_processed.lstar.zarr]
domain: genomics
source: "Seurat PBMC3k guided clustering tutorial (Satija Lab) — https://satijalab.org/seurat/articles/pbmc3k_tutorial — generalized to a tissue/species-agnostic recipe with data-driven QC thresholds and the v5 Louvain/uwot defaults."
---

# Single-sample scRNA-seq QC + clustering with R/Seurat (v5)

Generic single-sample recipe — works for **any tissue, any organism** as long as
the input is a raw integer counts matrix. The PBMC3k vignette is the reference
shape; every concrete threshold below (`nFeature` floors, `percent.mt` ceilings,
canonical marker panels) is **dataset-specific** and the recipe says, at each
step, what to look at to set sensible values for *your* data.

Prefer this (R / Seurat v5) when the session is already in R, the user names
Seurat, or downstream tools are Bioconductor. For SCTransform normalization in
the same workflow shape, see `seurat-sctransform`. For multi-sample
integration, see `seurat-integration`. For a Python-native equivalent
(scanpy), see `scrna-qc-clustering`.

**Verified against Seurat v5.x** — `FindClusters` defaults to Louvain
(`algorithm = 1`); `RunUMAP` defaults to `uwot`/cosine. Pin Seurat ≥ 5.0 if
reproducibility across sessions matters.

## Bundled references — load on demand

This SKILL.md is self-contained for the standard workflow. Load these only
when the task crosses into a topic the body summarizes but doesn't unpack:

- `references/installation_and_io.md` — reader function signatures
  (Read10X / Read10X_h5 / ReadMtx), 10x v2 vs v3 detection, gz handling,
  AnnData/h5Seurat bridges, organism-specific MT prefixes, gene name
  normalization gotchas, dependency table.
- `references/qc_and_thresholds.md` — QC metric definitions, how to read
  the quantile tables, tissue-specific cutoff conventions, the
  `subset()` filter, post-filter diagnostic, doublet detection (when to
  add scDblFinder/DoubletFinder), cell-cycle scoring.
- `references/figure_style.md` — the Seurat-collection figure
  conventions: palette (diverging blue/grey/red for signed, grey→red for
  sequential), `theme_cowplot()` ordering trap, alpha-poke pattern for
  DimPlot/FeaturePlot, dpi/bg/width conventions.
- `references/clustering_choices.md` — `FindClusters` algorithm enum
  (Louvain / Leiden / SLM), resolution semantics, `dims` choice off the
  elbow plot, UMAP defaults changes (v4 → v5), `min.dist` / `spread`
  tuning, when to skip clustering for reference mapping.
- `references/markers_and_annotation.md` — every `test.use` option
  (Wilcoxon / MAST / DESeq2 / ROC / LR / negbinom / poisson) with when
  to pick each, marker interpretation heuristics, the manual vs
  reference-mapped annotation paths, why DESeq2 is wrong for per-cell.

## Install

Idempotent — re-running is a no-op. Seurat is installed via the ABA capability
catalogue.

```r
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(patchwork)
  library(tidyr)
})

# Wilcoxon accelerator — Seurat v5 auto-uses presto (presto::wilcoxauc) for
# FindMarkers/FindAllMarkers; WITHOUT it, Wilcoxon falls back to a slow looping
# base-R implementation (>10 min vs <1 s on a small dataset). Shipped with Seurat
# by default (r-environment.yml); self-heal if a deployment predates that.
if (!requireNamespace("presto", quietly = TRUE)) {
  try(ensure_capability("presto"), silent = TRUE)   # prebuilt r-presto (bioconda)
}

stopifnot(packageVersion("Seurat") >= "5.0.0")  # this recipe is v5-only
```

`library(Seurat)` does NOT attach `ggplot2`, `dplyr`, or `cowplot` — load them
by name. There is no `tidyverse` meta-package in this environment.

For the full dependency table (what each package is for + alternatives), see
`references/installation_and_io.md`.

## Decisions to surface up front

Tell the user these are the analysis-defining choices:

1. **Species / MT prefix** — `^MT-` (human), `^mt-` (mouse), `^Mt-` or
   organism-specific for others. The full per-organism table is in
   `references/installation_and_io.md`. If `n_mt == 0` after Step 1, stop
   and fix.
2. **QC thresholds** — `nFeature_RNA` floor and ceiling, plus `percent.mt`
   ceiling. Read them off the per-metric quantile tables in Step 2 (not from
   a prior dataset's defaults). PBMC defaults: `nFeature 200–5000`,
   `percent.mt < 15`. Tissue-specific guidance in
   `references/qc_and_thresholds.md`.
3. **Number of PCs (`dims`)** — feeds `FindNeighbors` and `RunUMAP`. The
   PBMC3k vignette uses `1:10`; v5 routine practice is `1:30`. Default
   `DIMS_CHOSEN = 30` and adjust off the elbow if there is a clear cliff.
   Picking heuristics in `references/clustering_choices.md`.
4. **Clustering resolution** — `FindClusters(resolution = …)` controls
   cluster count. `0.5` is moderate; coarser/finer mapping in
   `references/clustering_choices.md`.
5. **Mode** — INTERACTIVE (default) when this is the primary analysis;
   BATCH when invoked as one of many sub-runs in a larger plan (e.g. one
   sample at a time across a GEO accession). The orchestrator passes
   `args="batch"`; the agent declares the mode in `present_plan`.

Figures the user will see as the analysis proceeds:
`qc_violins_pre.png`, `qc_scatters_pre.png`, `qc_violins_post.png`,
`hvg_plot.png`, `pca_elbow.png`, `pca_heatmap.png`, `umap_clusters.png`,
`markers_dotplot.png`, `markers_featureplot.png`.

---

## Step 1 — Load the counts matrix and create the Seurat object

Three common input shapes; pick by what's on disk.

```r
# 10x CellRanger directory: barcodes.tsv[.gz], features.tsv[.gz] (v3) or
# genes.tsv[.gz] (v2), matrix.mtx[.gz]. Gzipped files are auto-detected.
counts <- Read10X(data.dir = "/path/to/filtered_feature_bc_matrix")
```

```r
# 10x .h5 (filtered_feature_bc_matrix.h5):
counts <- Read10X_h5("/path/to/filtered_feature_bc_matrix.h5")
```

```r
# Loose, GSM-prefixed GEO triplets (Read10X will NOT find these — non-standard
# names). Use ReadMtx with explicit paths; feature.column = 2 → gene symbols
# (column 1 is Ensembl).
counts <- ReadMtx(
  mtx      = "/path/to/GSMxxxxxxx_matrix.mtx.gz",
  cells    = "/path/to/GSMxxxxxxx_barcodes.tsv.gz",
  features = "/path/to/GSMxxxxxxx_features.tsv.gz",
  feature.column = 2
)
```

Build the Seurat object with **minimal** pre-filtering — keep it permissive so
Step 2's QC plots show the full unfiltered distribution:

```r
obj <- CreateSeuratObject(
  counts       = counts,
  project      = "sample",          # short label used by orig.ident
  min.cells    = 3,                  # drop genes seen in <3 cells (sparsity, not biology)
  min.features = 200                 # drop barcodes with <200 genes (empty droplets)
)

# Report the load and the implicit min.cells/min.features filter delta.
n_genes_raw <- nrow(counts); n_cells_raw <- ncol(counts)
cat(sprintf("Pre-filter: %d genes x %d cells\n", n_genes_raw, n_cells_raw))
cat(sprintf("After CreateSeuratObject (min.cells=3, min.features=200): %d genes (-%d) x %d cells (-%d)\n",
            nrow(obj), n_genes_raw - nrow(obj),
            ncol(obj), n_cells_raw - ncol(obj)))

# Sanity-check the loaded object — symbol vs Ensembl, MT prefix coverage.
stopifnot(inherits(obj, "Seurat"))
head(rownames(obj))
n_mt <- sum(grepl("^MT-", rownames(obj)))
cat(sprintf("MT genes matched by '^MT-': %d  (0 means wrong prefix; mouse uses '^mt-')\n", n_mt))
```

**Report:** input format, raw cells/genes, cells/genes lost to the load filter,
whether MT prefix matched. If `n_mt == 0`, stop and fix the prefix before
computing `percent.mt` in Step 2.

`min.features = 200` is an empty-droplet floor — **not** a QC threshold. Real
QC happens in Step 2 against the actual distributions.

For reader signatures, gz auto-detection, AnnData / h5Seurat bridges, and
per-organism MT prefix table, read `references/installation_and_io.md`.

---

## Step 2 — Compute QC metrics and read the distributions

Compute first, then read the quantile tables BEFORE plotting, so thresholds
are picked from data and not eyeballed off a plot.

```r
# Mitochondrial fraction — adapt the pattern to species.
obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
# Ribosomal protein content — '^Rp[sl]' for mouse, etc.
obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

# Per-metric quantiles — pick floors from qs_lo, ceilings from qs_hi.
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
```

Quick reading guidance (full table of tissue-specific cutoffs in
`references/qc_and_thresholds.md`):

- **`nFeature_RNA` floor** — pick from `qs_lo` (1st–5th percentile is the
  empty-droplet shoulder). 200 is a typical PBMC floor; lower for low-RNA
  cells (neutrophils, erythrocytes).
- **`nFeature_RNA` ceiling** — pick from `qs_hi` (99th percentile usually
  sits where the doublet shoulder begins). 5000–8000 for 10x PBMC; higher
  for neurons.
- **`percent.mt` ceiling** — PBMC/blood 5–10%; solid tissue 15–25%;
  nuclei (snRNA-seq) ~5%.
- **`percent.ribo`** — don't filter unless visibly bimodal; high ribo is
  often biology, not failure.

Write the chosen thresholds into one small list so Step 3 can draw them as
red dashed lines AND apply them via `subset()`:

```r
THRESH <- list(
  nFeature_low  = 200,    # REPLACE — read off the quantile table above
  nFeature_high = 5000,   # REPLACE
  mt_high       = 15      # REPLACE
)
```

**Report:** the chosen thresholds and the quantile each one corresponds to.

For metric definitions, tissue-specific cutoffs, doublet detection, and
cell-cycle scoring, read `references/qc_and_thresholds.md`.

---

## Step 3 — Apply QC: plot the decision, filter, re-check

### 3a. Pre-filter plots — show the cells about to be removed

Color each cell by whether it will be kept or filtered; draw the threshold
values as red dashed lines. This shows whether the threshold is intercepting
the population you wanted to remove.

```r
# Single boolean: a cell passes ALL active criteria.
obj$qc_kept <- with(obj@meta.data,
  ifelse(nFeature_RNA > THRESH$nFeature_low  &
         nFeature_RNA < THRESH$nFeature_high &
         percent.mt   < THRESH$mt_high,
         "kept", "filtered"))

# --- pre-filter violin: jitter colored by qc_kept, thresholds drawn as red lines ---
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

# --- pre-filter scatters: colored by qc_kept, threshold lines in red ---
df <- as.data.frame(obj@meta.data[, c("nCount_RNA","nFeature_RNA",
                                      "percent.mt","qc_kept")])

s1 <- ggplot(df, aes(nCount_RNA, percent.mt, colour = qc_kept)) +
  geom_point(alpha = 0.10, size = 0.9) +
  geom_hline(yintercept = THRESH$mt_high, colour = "red",
             linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = c(kept = "black", filtered = "red"),
                      guide = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
  labs(title = sprintf("nCount vs percent.mt (cutoff %g%%)", THRESH$mt_high),
       x = "nCount_RNA", y = "percent.mt (%)", colour = NULL) +
  theme_cowplot()

s2 <- ggplot(df, aes(nCount_RNA, nFeature_RNA, colour = qc_kept)) +
  geom_point(alpha = 0.10, size = 0.9) +
  geom_hline(yintercept = c(THRESH$nFeature_low, THRESH$nFeature_high),
             colour = "red", linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = c(kept = "black", filtered = "red"),
                      guide = guide_legend(override.aes = list(alpha = 1, size = 1.5))) +
  labs(title = sprintf("nCount vs nFeature (cutoffs %g, %g)",
                       THRESH$nFeature_low, THRESH$nFeature_high),
       x = "nCount_RNA", y = "nFeature_RNA", colour = NULL) +
  theme_cowplot()

p_sc <- (s1 | s2) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "QC scatters, pre-filter",
                  theme = theme(plot.title = element_text(size = 13, face = "bold"))) &
  theme(legend.position = "bottom")

ggsave("qc_scatters_pre.png", plot = p_sc,
       bg = "white", dpi = 120, width = 12, height = 5.3)
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

If you lose more than ~20% of cells the thresholds are too tight — return to
Step 2 and re-read the quantiles.

### 3c. Post-filter violin — confirm no truncated peak at the cutoff

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

**Report:** cells before/after, percent removed, dominant failure mode
(low-`nFeature`, high-`mt`, or both), and any visible truncation at the
threshold in the post-filter violin.

For doublet detection (scDblFinder / DoubletFinder) — when to add it and
why it runs AFTER QC but BEFORE normalization — read
`references/qc_and_thresholds.md`. For the figure conventions used here
(palette, alpha, layout), read `references/figure_style.md`.

---

## Step 4 — Normalize, find variable features, scale

The canonical Seurat log-normalize path. (For low-depth or
heterogeneous-depth data, use `seurat-sctransform` instead — that recipe
replaces this whole step with one `SCTransform()` call.)

```r
obj <- NormalizeData(obj,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000,
                     verbose = FALSE)

obj <- FindVariableFeatures(obj,
                            selection.method = "vst",
                            nfeatures = 2000,
                            verbose = FALSE)

# Sanity-check the HVG list — if MT-/RPL/RPS/MALAT1 dominate, QC was too loose.
top20 <- head(VariableFeatures(obj), 20)
print(top20)
suspect <- grep("^(MT-|RPS|RPL|MALAT1)", top20, value = TRUE)
cat("QC-warning HVGs:", if (length(suspect)) paste(suspect, collapse = ", ") else "none", "\n")
```

HVG mean-variance plot. `VariableFeaturePlot` has no `alpha` argument — we
walk `p$layers` and set `aes_params$alpha` on the `GeomPoint` (public ggplot2
API; alpha-poke pattern documented in `references/figure_style.md`).

```r
# ggrepel is a transitive dependency of Seurat's LabelPoints; install BEFORE
# we build the plot that uses it.
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")

top10 <- head(VariableFeatures(obj), 10)
p_hvg <- VariableFeaturePlot(obj)
for (k in seq_along(p_hvg$layers)) {
  if (inherits(p_hvg$layers[[k]]$geom, "GeomPoint")) {
    p_hvg$layers[[k]]$aes_params$alpha <- 0.35
  }
}
p_hvg <- LabelPoints(p_hvg, points = top10, repel = TRUE)
p_hvg <- p_hvg +
  ggtitle(sprintf("HVG selection (vst, top %d of %d genes); top 10 labeled",
                  length(VariableFeatures(obj)), nrow(obj))) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"))

ggsave("hvg_plot.png", plot = p_hvg,
       bg = "white", dpi = 120, width = 8, height = 5.5)
```

Then scale — the default (HVGs only) is what PCA needs:

```r
obj <- ScaleData(obj, verbose = FALSE)
# Pass features = rownames(obj) ONLY if you need scaled values for a non-HVG
# gene later (e.g. plotting a specific marker on a heatmap). Doubles ScaleData
# time without changing PCA results.
```

Regress out unwanted variation only when there's evidence it matters
(e.g. an explicit cell-cycle confound visible in Step 5's PC loadings):

```r
# obj <- ScaleData(obj, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
#                  verbose = FALSE)
```

**Report:** top 20 HVGs, any QC-warning HVGs flagged, total variable features.

### Variants — alternative normalization in the same slot

`NormalizeData(..., normalization.method = ...)` accepts:

- `"LogNormalize"` (default; the canonical block above).
- `"CLR"` — centered log ratio. Use for CITE-seq protein assays, NOT RNA.
- `"RC"` — relative counts (no log). Rarely useful for standard scRNA-seq.

For SCTransform (the full regularized-NB alternative that replaces
`NormalizeData` + `FindVariableFeatures` + `ScaleData`), use
`seurat-sctransform` — it changes the downstream assay slot and DE
preparation, so it's a separate recipe rather than a variant within this
step.

---

## Step 5 — PCA and choose `dims`

```r
obj <- RunPCA(obj,
              features = VariableFeatures(obj),
              npcs = 50,
              verbose = FALSE)

# Print top loadings on PCs 1-5 — sanity-check that PC1 is biology, not
# MT/ribo/cell-cycle.
print(obj[["pca"]], dims = 1:5, nfeatures = 8)
```

PCA elbow — per-PC variance (solid blue) and cumulative (grey dashed) on the
same axis, normalized to the **total HVG-matrix variance** (not the sum of
the 50 PCs' variance — the cumulative curve then plateaus at the true
fraction the PCs capture, ~15–30% typical for a single 10x sample).

```r
DIMS_CHOSEN <- 30   # default; override if the elbow shows a clear cliff

sd_mat    <- GetAssayData(obj, layer = "scale.data")
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
  labs(title = "PCA variance explained (of total HVG-matrix variance)",
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

DimHeatmap — visualize each PC's top-loading cells × top genes. Useful to
spot the PC where structure becomes noise and to detect MT/ribo/cell-cycle
domination of an early PC.

```r
# fast = FALSE returns real ggplots so the diverging scale applies (the
# fast = TRUE default uses image() and ignores ggplot scales).
p_dimheat <- DimHeatmap(obj, dims = 1:10, cells = 500,
                        balanced = TRUE, fast = FALSE, combine = TRUE,
                        ncol = 5) &
  scale_fill_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                       midpoint = 0) &
  theme(legend.position = "none")  # drop 10 identical per-panel legends

suppressMessages(
  ggsave("pca_heatmap.png", plot = p_dimheat,
         bg = "white", dpi = 180, width = 20, height = 9)
)
```

**Report:** `DIMS_CHOSEN`, the top-loading genes on PCs 1–5, and any PC that
looks technical (MT, ribo, cell cycle).

If PC1 is driven by `MT-` genes, ribosomal genes, or cell-cycle markers, the
options are: (a) tighten QC, (b) regress the offending score in `ScaleData`
(see Step 4), or (c) exclude that PC from `dims` (e.g. `dims = 2:30`).

For elbow-reading heuristics, the JackStraw alternative, and `dims` choice
in detail, read `references/clustering_choices.md`.

---

## Step 6 — Neighbor graph, clustering, UMAP

```r
obj <- FindNeighbors(obj, dims = 1:DIMS_CHOSEN, verbose = FALSE)

# FindClusters defaults to Louvain (algorithm = 1). For Leiden, pass
# algorithm = 4 (requires the leidenalg Python package via reticulate).
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

obj <- RunUMAP(obj, dims = 1:DIMS_CHOSEN, verbose = FALSE)

cat("Clusters:", length(levels(Idents(obj))), "\n")
print(table(Idents(obj)))
```

UMAP coloured by cluster, labels drawn on the embedding. DimPlot doesn't
expose `alpha` — alpha-poke its `GeomPoint` layer for readable density.

```r
p_umap <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE,
                  pt.size = 0.4) +
  ggtitle(sprintf("Louvain res=0.5 . dims=1:%d . n=%d cells . %d clusters",
                  DIMS_CHOSEN, ncol(obj), length(levels(Idents(obj))))) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold")) +   # 14pt default crops the title at width=7
  NoLegend() +              # MUST come after theme_cowplot — the theme overrides legend.position
  coord_fixed()

for (k in seq_along(p_umap$layers)) {
  if (inherits(p_umap$layers[[k]]$geom, "GeomPoint")) {
    p_umap$layers[[k]]$aes_params$alpha <- 0.6
  }
}

ggsave("umap_clusters.png", plot = p_umap,
       bg = "white", dpi = 120, width = 7, height = 6.5)
```

**Report:** algorithm + resolution used, cluster count, cluster-size table,
any cluster <1% of cells (often a doublet bridge or rare population to
flag).

For the algorithm enum (`1`/`2`/`3`/`4` = Louvain/multilevel/SLM/Leiden),
resolution semantics + cluster-count mapping, `dims` choice heuristics,
the v4 → v5 UMAP default changes, and `min.dist`/`spread` tuning, read
`references/clustering_choices.md`.

---

## Step 7 — Find cluster markers

`FindAllMarkers` runs a per-cluster vs.-rest Wilcoxon test by default.

```r
markers <- FindAllMarkers(obj,
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

> **Speed.** Seurat v5 auto-detects `presto` and uses it for Wilcoxon
> (~10× faster, identical results). `install.packages("presto")` if missing.

Two figures — dotplot of top 5 per cluster, and FeaturePlot of canonical
lineage markers. **Don't paint a wall of FeaturePlots over the full top-N
list** — that overwhelms more than it informs.

```r
# Figure 1: dotplot, top 5 markers per cluster.
genes_to_show <- unique(top5$gene)

p_dot <- DotPlot(obj, features = genes_to_show, cluster.idents = FALSE) +
  RotatedAxis() +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0, name = "avg expr") +
  ggtitle(sprintf("Top 5 markers per cluster (Wilcoxon, only.pos, n=%d genes)",
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
# Figure 2: FeaturePlot, 4-8 hand-picked canonical lineage markers (one per
# major population). REPLACE the panel with markers for YOUR tissue/organism
# from a curated reference (PanglaoDB, CellMarker, or a published atlas for
# the same organ). The list below is a PBMC placeholder.
canonical <- c("CD3D","CD8A","MS4A1","GNLY","LYZ","PPBP")
canonical <- canonical[canonical %in% rownames(obj)]

p_feat <- FeaturePlot(obj, features = canonical, order = TRUE,
                      pt.size = 0.3, ncol = 3) &
  scale_colour_gradient(low = "grey85", high = "#b2182b") &
  theme_cowplot() &
  theme(legend.position = "right", legend.key.size = unit(0.4, "cm"))

# Alpha-poke each panel's GeomPoint (FeaturePlot lacks alpha arg).
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

**Report:** total significant markers, top 5 per cluster (gene, avg_log2FC,
pct.1, pct.2, p_val_adj), any cluster with no genes meeting threshold (often
under-clustering — neighboring biology bleeds across the boundary).

### Annotation — biology, not numbers

Cluster labels come from biology. Two paths:

1. **Manual** — check top markers against a curated panel for the tissue
   (PanglaoDB, CellMarker, or a published atlas). Re-label with
   `obj <- RenameIdents(obj, "0" = "CD14+ mono", ...)`.
2. **Automated reference mapping** — `SingleR`, `Azimuth`, or `scArches` if
   a reference exists for the tissue/organism.

**Do not invent cell-type names from a single marker.** "FCGR3A+ monocytes"
needs FCGR3A *plus* the rest of the non-classical-monocyte signature.

For the full `test.use` enum (Wilcoxon / MAST / DESeq2 / ROC / LR /
negbinom / poisson) with when-to-use-each, marker interpretation
heuristics, and the manual vs reference-mapped annotation paths, read
`references/markers_and_annotation.md`.

---

## Step 8 — Save & deliver: object, viewer store, and interactive link

The closing step *delivers* the result — save the object, write the viewer store,
and hand the user a **clickable way to explore it**, not just files on disk.
Present this as the final "save & deliver" step of your plan; the interactive
viewer link comes from here, so don't drop it.

```r
# 1. Save the processed object.
saveRDS(obj, file = "seurat_processed.rds")
sz <- file.info("seurat_processed.rds")$size / 1e6
cat(sprintf("Wrote seurat_processed.rds (%.1f MB)\n", sz))
obj_check <- readRDS("seurat_processed.rds")
stopifnot(identical(dim(obj_check), dim(obj)))
rm(obj_check); invisible(gc())

# 2. Write the pagoda3 viewer store DIRECTLY from the live Seurat object with
#    lstar — pure R, highest fidelity (all reductions/metadata carried across).
#    Do NOT route through .h5ad: zellkonverter/sceasy spin up a basilisk/reticulate
#    Python env — slow and pointless just to view the result.
d <- lstar::read_seurat(obj)
lstar::lstar_write_viewer(d, "seurat_processed.lstar.zarr")   # viewer=TRUE: precomputes
# DE / variable-genes / cell-major counts in R so pagoda3 opens it *optimized* (no
# "Not viewer-optimized" banner). Needs the clustering set as Idents (Step 6). Plain
# lstar_write() also works but leaves the store un-optimized.
```

Then **call `open_viewer(file_path="seurat_processed.lstar.zarr")` and present the
returned link in your closing message** — a required part of delivering the result,
not optional. Notes:
- You *can* instead hand `open_viewer` the `seurat_processed.rds` directly (ABA
  converts on launch), but that's a lower-fidelity fallback for installs without
  the R stack — prefer the in-session `.lstar.zarr`.
- Export `.h5ad` only when the target is a *different* tool (scanpy, cellxgene),
  never as the route to the ABA viewer.
- If `open_viewer` returns `ok:false`, relay the error rather than handing out a
  dead link.

The `.rds` carries the full `RNA` assay (counts, data, scale.data layers), the
`pca`/`umap` reductions, the `RNA_nn`/`RNA_snn` graphs, `seurat_clusters`, and the
HVGs — a follow-up session resumes from markers/sub-clustering/annotation without
re-running 1–6.

> **Faster I/O on large objects.** `qs::qsave(obj, "seurat_processed.qs")` is
> 2–3× faster than `saveRDS` for big objects and produces smaller files.
> **Handing off to scanpy / Python?** `zellkonverter` (Bioconductor) via the
> `SingleCellExperiment` bridge gives a standard `.h5ad` — that's for *other tools*,
> not the ABA viewer. Details in `references/installation_and_io.md`.

---

## Batch variant — use INSTEAD of Steps 1–8 when invoked with args="batch"

Branch on `$ARGUMENTS == "batch"` at the top of the body. In batch mode the
canonical step-by-step + figure + report cadence is replaced by ONE
consolidated call:

- No `ggsave` / `png` calls (skip all per-step figures).
- No per-step quantile tables, no top-HVG printout, no top-marker print.
- Still save the canonical processed object (`seurat_processed.rds`) AND
  the marker table (`cluster_markers.csv`) — those are what the
  orchestrator's rollup steps consume.
- Print ONE final summary line of the form
  `"batch ok | <N> cells | <K> clusters | <M> markers | <MB> MB"`.

The batch caller is responsible for supplying QC thresholds (and any
`vars.to.regress`) up front — defaults below come from the interactive
recipe.

```r
# args="batch": single consolidated call, no figures, no narration.
# Caller passes input path + thresholds + dims/resolution via $ARGUMENTS.
counts <- Read10X(data.dir = "/path/to/filtered_feature_bc_matrix")
obj <- CreateSeuratObject(counts = counts, project = "sample",
                          min.cells = 3, min.features = 200)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

obj <- subset(obj, subset =
  nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, verbose = FALSE)

DIMS_CHOSEN <- 30
obj <- FindNeighbors(obj, dims = 1:DIMS_CHOSEN, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:DIMS_CHOSEN, verbose = FALSE)

markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)

saveRDS(obj, "seurat_processed.rds")
write.csv(markers, "cluster_markers.csv", row.names = FALSE)

sz <- file.info("seurat_processed.rds")$size / 1e6
cat(sprintf("batch ok | %d cells | %d clusters | %d markers | %.1f MB\n",
            ncol(obj), length(levels(Idents(obj))), nrow(markers), sz))
```

---

## Final response checklist

Summarize:

- input format (10x dir / .h5 / .mtx triplet), sample ID, and MT-prefix sanity
- raw cells/genes and the load filter delta
- chosen QC thresholds, percent of cells removed, dominant failure mode
- HVG list size, any QC-warning HVGs flagged
- `DIMS_CHOSEN` and any PC that looked technical
- clustering algorithm + resolution, cluster count, cluster-size summary
- total significant markers, top markers per cluster (subset for chat),
  any cluster with zero markers
- figures shown to the user (filenames)
- saved files (`seurat_processed.rds`, `cluster_markers.csv`, `seurat_processed.lstar.zarr`)
- **the interactive viewer link** from Step 8 (`open_viewer` on the
  `.lstar.zarr`) — always include it; if you couldn't produce one, say why
- caveats: doublet detection not run, batch effects (single-sample so
  irrelevant unless the sample is itself a multiplexed pool), weak markers,
  over/under-clustering

---

## See also

- `seurat-sctransform` — same workflow shape with `SCTransform()` replacing
  NormalizeData + FindVariableFeatures + ScaleData. Switch to it when the
  sample is low-depth or has heterogeneous sequencing depth across cell
  types, or when the user names SCTransform / regularized negative binomial.
- `seurat-integration` — multi-sample workflow (split layers, run
  `IntegrateLayers`, cluster on the integrated reduction). Run THIS recipe
  per-sample first to sanity-check each, then integrate.
- `scrna-qc-clustering` — scanpy / Python equivalent for a Python-native
  session.
