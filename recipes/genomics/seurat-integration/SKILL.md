---
name: seurat-integration
description: Multi-sample scRNA-seq integration with R/Seurat v5 — load N samples, merge, split the RNA assay into per-sample layers, preprocess (Normalize→HVG→Scale→PCA) on the layered object, then IntegrateLayers (CCAIntegration / RPCAIntegration / HarmonyIntegration / JointPCAIntegration / scVIIntegration / FastMNNIntegration) into a corrected reduction. Cluster + UMAP on the integrated reduction; JoinLayers before DE. The Seurat v5 layer-based replacement for v4's FindIntegrationAnchors / IntegrateData anchor flow.
when_to_use: Two or more scRNA-seq samples / donors / conditions / batches (10x lanes, GEO GSMs, stim-vs-ctrl) whose batch effect is visible in a plain PCA/UMAP, and you want one shared embedding before clustering, in R/Seurat. Use THIS when the session is R, or when the user names Seurat v5 / IntegrateLayers / CCA / RPCA / anchors / layer-based integration. For a Python/scanpy session use harmony-integration-scanpy; for atlas-scale deep-generative integration see scvi-integration; a single clean sample needs no integration — see seurat-scrna-v2.
avoid_when: You have one sample (use seurat-scrna-v2) or you must reproduce a pre-v5 anchor analysis (use the v4 FindIntegrationAnchors / IntegrateData flow noted in references/integration_methods.md).
invocation: interactive+batch
requires_tools: [run_r]
capabilities_needed: [Seurat]
keywords: [Seurat, Seurat v5, integration, IntegrateLayers, CCAIntegration, RPCAIntegration, HarmonyIntegration, JointPCAIntegration, FastMNNIntegration, scVIIntegration, batch correction, batch effect, multi-sample, multi sample, multiple samples, sample integration, layer-based, split layers, JoinLayers, integrated UMAP, integrated.cca, integrated.rpca, harmony, anchors, ReadMtx, Read10X, GEO, R, scRNA-seq, single cell]
produces: [umap_unintegrated_by_sample.png, umap_integrated_by_sample.png, umap_integrated_by_cluster.png, pca_elbow.png, seurat_integrated.rds]
domain: genomics
resource_profile: "medium (~3–10 min for 2–6 samples × ~5–15k cells each; CCA scales with cell count, RPCA/Harmony faster)"
source: "Seurat v5 layer-based integration vignette (Satija Lab) — satijalab.org/seurat/articles/integration_introduction; complemented by ?IntegrateLayers reference manual."
---

# Multi-sample scRNA-seq integration with R/Seurat (v5)

Integration aligns several scRNA-seq samples into **one shared low-dimensional
embedding** so that the same cell types from different samples co-embed, while
distinct cell types stay separate — letting you cluster, annotate, and compare
across the whole set instead of per-sample. This is the **Seurat v5**
layer-based flow: keep ONE merged object whose RNA counts are split into one
layer per sample, preprocess it like a single object, and let
`IntegrateLayers` learn a batch-corrected reduction. It replaces v4's
`FindIntegrationAnchors` / `IntegrateData` pipeline.

This recipe assumes you already know each sample loads cleanly — run
`seurat-scrna-v2` per sample first if you have not. Integrating before
per-sample QC hides sample-specific quality failures.

**Seurat version requirement: v5.0 or newer.** `IntegrateLayers`, layer-based
assays, and `JoinLayers` do not exist in v4. The recipe will fail loudly on v4
— the install step asserts the version.

## Bundled references — load on demand

This recipe is self-contained for the standard CCA-on-LogNormalize workflow.
For deeper detail, load the matching reference file ONLY when needed:

- `references/integration_methods.md` — full per-method comparison
  (CCAIntegration / RPCAIntegration / HarmonyIntegration / JointPCAIntegration /
  FastMNNIntegration / scVIIntegration), costs, when-to-pick rules, plus the
  classic v4 `FindIntegrationAnchors`/`IntegrateData` fallback.
- `references/layered_assay_model.md` — Seurat v5 split-layer semantics, why
  `JoinLayers` is required before `split()` on merged objects, layer naming
  conventions, JoinLayers timing relative to DE.
- `references/over_correction_diagnostics.md` — how to read the
  by-sample / by-cluster UMAPs, mixing metrics, marker preservation,
  patches when integration washed out real biology or left visible batch
  effect.
- `references/figure_style.md` — Seurat-collection figure conventions
  (palette, alpha-poke, `theme_cowplot()` placement vs `NoLegend()`, ggsave
  defaults, `coord_fixed()`). Same content across the seurat-* collection;
  duplicated for now.

## Install

Idempotent — re-running is a no-op once Seurat v5 is present.

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("Seurat",   quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2",  quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",    quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("cowplot",  quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("patchwork",quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("tidyr",    quietly = TRUE)) install.packages("tidyr")

# Hard pin: this recipe is v5-only. IntegrateLayers / split-RNA-assay / JoinLayers
# are Seurat-v5 APIs and have no v4 equivalent.
stopifnot(packageVersion("Seurat") >= "5.0.0")
```

Then attach the libraries once (the kernel session is persistent across
subsequent `run_r` calls; do not re-attach in every step):

```r
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(cowplot)
})
```

**Parallelism (`future`).** ABA defaults R `future` to **sequential** with a
generous `future.globals.maxSize`, so `IntegrateLayers` won't hit future's 500 MiB
globals cap on real data (the classic *"total size of globals exceeds the maximum
allowed size"* error — IntegrateLayers ships the layered object to workers). Most
runs should stay sequential. **To opt into parallelism** (large object, and you
were allocated multiple cores — CCA/RPCA anchor-finding benefits most):

```r
future::plan("multicore", workers = 4)   # N = cores allocated; needs RAM × workers
# … IntegrateLayers(...) …
future::plan("sequential")               # revert afterwards
```

The globals cap is already raised by ABA, so opting in won't trip it; if you ever
see the cap error anyway, `options(future.globals.maxSize = 8 * 1024^3)`.

Optional per-method dependencies — install ONLY if the user picks that
branch. For full per-method install + tradeoffs see
`references/integration_methods.md`.

```r
# HarmonyIntegration — CRAN harmony
# if (!requireNamespace("harmony", quietly = TRUE)) install.packages("harmony")
# FastMNNIntegration  — Bioc batchelor + SeuratWrappers
# scVIIntegration     — reticulate + scvi-tools Python env + SeuratWrappers
```

## Decisions to surface up front

Surface these in `present_plan` BEFORE running anything. Over-integration
silently erases real biology, so this is where you want the user (or an
advisor) to confirm.

1. **What to integrate over (the split variable)** — the per-sample / batch
   covariate you split layers on. Split on the **technical** nuisance
   (sample, donor, lane, batch, 10x run), **never** on the biological
   variable you intend to test. If `stim` vs `ctrl` IS the question,
   integrating over it washes out the very signal you're studying.
2. **Method** — `CCAIntegration` (default; most cell types shared across
   samples; can over-correct rare populations), `RPCAIntegration`
   (faster, more conservative; preferred for large or partly-shared
   cohorts), `HarmonyIntegration` (fast linear correction of the existing
   PCA), `JointPCAIntegration` (joint PCA across all samples), or
   `scVIIntegration` / `FastMNNIntegration` (companion-package wrappers).
   Full picker in `references/integration_methods.md`.
3. **Normalization path** — `LogNormalize` (default; `NormalizeData` +
   `FindVariableFeatures` + `ScaleData` + `RunPCA`) or `SCTransform`
   (single call replacing the first three). Pick LogNormalize unless the
   user names SCT or the dataset has very high cell-type diversity; SCT
   is heavier and requires `PrepSCTFindMarkers` before DE.
4. **Number of PCs (`dims`)** — feeds the integration, the neighbor
   graph, and the integrated UMAP. The vignette uses `1:30`. Read the
   elbow off `ElbowPlot(obj)` after PCA; 20–30 is safe when the elbow is
   unclear. Use the SAME `dims` everywhere downstream.
5. **Clustering resolution** — `FindClusters(resolution = …)` on the
   integrated SNN graph. 0.5–1.0 is the typical range; 1.0 is the
   vignette's value.
6. **Mode** — INTERACTIVE (default) when this is the primary analysis;
   BATCH when invoked as one of many sub-runs in a larger plan. The
   orchestrator passes `args="batch"`; the agent declares the mode in
   its `present_plan` so the user sees which one is in effect. See the
   `## Batch variant` section below.

Show the user these figures as the analysis proceeds:
- `umap_unintegrated_by_sample.png` — pre-integration, colored by sample
- `pca_elbow.png` — per-PC and cumulative variance
- `umap_integrated_by_sample.png` — post-integration, colored by sample
- `umap_integrated_by_cluster.png` — post-integration, colored by cluster

---

## Step 1 — Load each sample, merge, split into per-sample layers

Load every sample as its own Seurat object, then `merge()` into one and split
the RNA assay so v5 holds one `counts.<sample>` layer per sample. That split
is exactly what `IntegrateLayers` reads.

### Standard CellRanger output per sample (one filtered_feature_bc_matrix directory per sample)

```r
# One sample per CellRanger output directory. `Read10X` expects barcodes.tsv[.gz],
# features.tsv[.gz] (or genes.tsv[.gz] for V2), matrix.mtx[.gz] in the directory.
sample_dirs <- c(
  ctrl = "/path/to/ctrl/filtered_feature_bc_matrix",
  stim = "/path/to/stim/filtered_feature_bc_matrix"
)

objs <- lapply(names(sample_dirs), function(s) {
  counts <- Read10X(data.dir = sample_dirs[[s]])
  CreateSeuratObject(counts = counts, project = s,
                     min.cells = 3, min.features = 200)
})
names(objs) <- names(sample_dirs)
```

### Variant: GEO / SRA loose-triplet bundle — use INSTEAD of the standard block above when the filenames are GSM-prefixed (not a CellRanger directory)

GEO supplementary downloads typically deliver each sample as three
GSM-prefixed files in a shared directory. `Read10X` will NOT find these.
Use `ReadMtx` with explicit per-file paths.

```r
# feature.column = 2 reads symbols; column 1 is Ensembl IDs in a standard 10x triplet.
samples <- list(
  ctrl = list(
    mtx      = "/path/to/bundle/GSM5746268_ctrl_matrix.mtx.gz",
    cells    = "/path/to/bundle/GSM5746268_ctrl_barcodes.tsv.gz",
    features = "/path/to/bundle/GSM5746268_ctrl_features.tsv.gz"),
  stim = list(
    mtx      = "/path/to/bundle/GSM5746269_stim_matrix.mtx.gz",
    cells    = "/path/to/bundle/GSM5746269_stim_barcodes.tsv.gz",
    features = "/path/to/bundle/GSM5746269_stim_features.tsv.gz")
)

objs <- lapply(names(samples), function(s) {
  f <- samples[[s]]
  m <- ReadMtx(mtx = f$mtx, cells = f$cells, features = f$features,
               feature.column = 2)
  CreateSeuratObject(counts = m, project = s,
                     min.cells = 3, min.features = 200)
})
names(objs) <- names(samples)
```

### Merge into one object and split the RNA assay by sample

```r
# add.cell.ids prefixes barcodes so they stay unique across samples.
obj <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
obj$sample <- obj$orig.ident          # technical covariate to integrate over

# JoinLayers FIRST: merge() of v5 objects already produces per-input split
# layers; split() errors if any layers are already split. JoinLayers is a
# no-op when nothing is split — safe either way. See
# references/layered_assay_model.md for why.
obj <- JoinLayers(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
```

**Sanity-check the loaded object:**

```r
# After merge + split, the RNA assay should be v5 (Assay5) with one
# counts.<sample> layer per sample.
stopifnot(inherits(obj[["RNA"]], "Assay5"))
cat(sprintf("Merged %d samples | %d cells x %d genes total\n",
            length(objs), ncol(obj), nrow(obj)))
print(table(sample = obj$sample))                  # cells per sample
print(Layers(obj[["RNA"]]))                        # one counts.<sample> per sample
```

**Report:** number of samples, cells per sample, total cells × genes, and a
confirmation that the RNA assay split landed (one `counts.<sample>` layer per
sample). If a sample loaded with <500 cells or >50,000 cells, flag it.

**Pitfall — per-sample QC.** This recipe loads with permissive
`min.cells = 3` / `min.features = 200` only — that's a sparse-droplet floor,
NOT QC. Run real QC (`percent.mt`, doublets, cell-count cutoffs) per sample
BEFORE integrating, via `seurat-scrna-v2`. Integrating cells with bad QC is
the single most common way to get garbage integrated UMAPs.

For the v5 split-layer model and `JoinLayers` semantics in depth, read
`references/layered_assay_model.md`.

---

## Step 2 — Preprocess the merged layered object → PCA

With layers split, the four standard preprocessing steps run **per layer
automatically** — no per-sample loop needed. `ScaleData` MUST run before
`RunPCA`.

```r
obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                     scale.factor = 10000, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst",
                            nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
```

### Variant: SCTransform — use INSTEAD of NormalizeData + FindVariableFeatures + ScaleData when the user names SCT or the dataset is very heterogeneous

```r
# Replaces all three Log-normalize / HVG / Scale steps with one call.
# IntegrateLayers in Step 4 then needs normalization.method = "SCT".
obj <- SCTransform(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
```

### PCA elbow — pick `dims` once

Two curves on one axis: per-PC variance (solid blue) and cumulative variance
(grey dashed), as percent of total HVG-matrix variance. Default
`DIMS_CHOSEN = 30` matches the vignette and is robust; override only if the
elbow gives a strong reason.

```r
DIMS_CHOSEN <- 30   # default; override if the elbow says otherwise

sd_mat    <- GetAssayData(obj, layer = "scale.data")
total_var <- sum(apply(sd_mat, 1, var))
var_pc    <- Stdev(obj, reduction = "pca")^2
prop      <- var_pc / total_var
cum       <- cumsum(prop)
df_elb    <- data.frame(PC = seq_along(prop), prop = prop, cum = cum)
df_long   <- tidyr::pivot_longer(df_elb, c(prop, cum),
                                 names_to = "kind", values_to = "value")
df_long$kind <- factor(df_long$kind, levels = c("prop","cum"),
                       labels = c("per-PC","cumulative"))

p_elbow <- ggplot(df_long, aes(PC, value, colour = kind, linetype = kind)) +
  geom_vline(xintercept = DIMS_CHOSEN, colour = "red",
             linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = DIMS_CHOSEN - 0.7, y = max(df_long$value) * 0.95,
           label = sprintf("dims = 1:%d (chosen)", DIMS_CHOSEN),
           colour = "red", hjust = 1, size = 3.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 0.4, 0.05)) +
  scale_x_continuous(breaks = c(1, seq(5, 50, 5))) +
  scale_colour_manual(values = c("per-PC" = "#1f77b4",
                                 "cumulative" = "grey40")) +
  scale_linetype_manual(values = c("per-PC" = "solid",
                                   "cumulative" = "dashed")) +
  labs(title = "PCA variance explained (of total HVG-matrix variance)",
       x = "principal component", y = "proportion of total variance",
       colour = NULL, linetype = NULL) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = c(0.98, 0.55),
        legend.justification = c(1, 0.5),
        legend.background = element_rect(fill = alpha("white", 0.5), colour = NA),
        legend.title = element_blank())

ggsave("pca_elbow.png", p_elbow,
       width = 7, height = 4.8, units = "in", dpi = 120, bg = "white")
```

**Report:** PCs chosen, per-PC variance of PC1–PC5, cumulative variance at
`DIMS_CHOSEN`. Typical cumulative ranges 10–30% depending on cell-type
diversity (high heterogeneity spreads variance across more PCs). If PC1 is
dominated by `MT-` / ribosomal / cell-cycle genes (from
`print(obj[["pca"]], dims = 1:5, nfeatures = 8)`), tighten per-sample QC or
regress those scores in `ScaleData`.

For shared figure style (palette, alpha-poke, ggsave defaults) see
`references/figure_style.md`.

---

## Step 3 — Pre-integration UMAP: show the batch effect FIRST

A UMAP on the raw `pca` reduction, colored by sample. If samples form
separate islands per cell type, there's a batch effect worth integrating.
Keep this reduction under its own name (`umap.unintegrated`) so it does not
clobber the integrated one in Step 5.

```r
obj <- RunUMAP(obj, reduction = "pca", dims = 1:DIMS_CHOSEN,
               reduction.name = "umap.unintegrated", verbose = FALSE)

p_pre <- DimPlot(obj, reduction = "umap.unintegrated", group.by = "sample",
                 pt.size = 0.4) +
  ggtitle(sprintf("Before integration — colored by sample (n=%d cells, %d samples)",
                  ncol(obj), length(unique(obj$sample)))) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4)) +
  coord_fixed()

# DimPlot does not expose alpha; walk p$layers for GeomPoint and poke.
for (k in seq_along(p_pre$layers)) {
  if (inherits(p_pre$layers[[k]]$geom, "GeomPoint")) {
    p_pre$layers[[k]]$aes_params$alpha <- 0.6
  }
}

ggsave("umap_unintegrated_by_sample.png", p_pre,
       width = 7, height = 6, units = "in", dpi = 120, bg = "white")
```

**Report:** is the batch effect visible? Sample-segregated islands per cell
type → yes, integrate. Already well-mixed → integration is optional and may
be harmful.

---

## Step 4 — IntegrateLayers (the integration step)

One call learns a batch-corrected reduction from the split layers. `method`
is the **unquoted function reference** (not a string), `orig.reduction` is
the PCA computed in Step 2, and `new.reduction` is the name to assign the
corrected embedding.

### Default — CCAIntegration (LogNormalize path)

```r
# method is an UNQUOTED FUNCTION REFERENCE — NOT a string.
# Common args (formals(IntegrateLayers) verified against Seurat 5.5.0):
#   object         the layered Seurat object
#   method         function reference (CCAIntegration, RPCAIntegration, …)
#   orig.reduction "pca" (from Step 2; default in formals)
#   new.reduction  any name; downstream steps point at this
#   verbose        FALSE for quieter output
obj <- IntegrateLayers(
  object         = obj,
  method         = CCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.cca",
  verbose        = FALSE
)
```

### Variant: RPCAIntegration — use INSTEAD of the canonical block above when the dataset is large or populations are only partly shared

```r
# Faster and more conservative than CCA. No extra deps.
obj <- IntegrateLayers(
  object         = obj,
  method         = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.rpca",
  verbose        = FALSE
)
```

### Variant: HarmonyIntegration — use INSTEAD of the canonical block above when fast linear correction is sufficient

```r
# Requires the `harmony` package (install above).
# Harmony's own new.reduction default is "harmony"; we keep that name.
obj <- IntegrateLayers(
  object         = obj,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  verbose        = FALSE
)
```

For the full method picker (CCA / RPCA / Harmony / JointPCA / scVI /
FastMNN), per-method costs, when to pick which, and the SCTransform
normalization wiring, read `references/integration_methods.md`.

**Report:** which method was used, source reduction (`orig.reduction`),
destination reduction name (`new.reduction`), and a one-line confirmation
that the new reduction landed (`"integrated.cca" %in% Reductions(obj)`).

**Pitfall — `method` is a function, not a string.**
`method = "CCAIntegration"` fails with a confusing `could not find function`
error. Pass the unquoted name.

**Pitfall — `new.reduction` is the contract.** Whatever you set here is the
name every downstream step (`FindNeighbors`, `RunUMAP`) must point at. If
you switch methods, use a different `new.reduction` so the reductions sit
side-by-side and you can compare. Do NOT reuse `"pca"` as `new.reduction` —
that overwrites the raw PCA and you lose the pre-integration baseline.

---

## Step 5 — Cluster and UMAP on the INTEGRATED reduction

The only change vs an unintegrated run: every neighbor / UMAP / cluster
call now points at the integrated reduction from Step 4. Use the same
`dims` you picked in Step 2.

```r
INTEG_REDUC <- "integrated.cca"   # match the new.reduction from Step 4

obj <- FindNeighbors(obj, reduction = INTEG_REDUC, dims = 1:DIMS_CHOSEN,
                    verbose = FALSE)
obj <- FindClusters(obj, resolution = 1.0,
                    cluster.name = "integrated_clusters", verbose = FALSE)
obj <- RunUMAP(obj, reduction = INTEG_REDUC, dims = 1:DIMS_CHOSEN,
               reduction.name = "umap.integrated", verbose = FALSE)

cat(sprintf("Integrated clusters: %d (resolution = 1.0, dims = 1:%d)\n",
            length(unique(obj$integrated_clusters)), DIMS_CHOSEN))
print(table(cluster = obj$integrated_clusters))
```

**Pitfall — `FindClusters` defaults to Louvain (`algorithm = 1`).** For
Leiden pass `algorithm = 4` — but Leiden requires `leidenalg` via
`reticulate` and a configured Python env; on error fall back to Louvain by
dropping the `algorithm` arg. Don't silently switch between Louvain and
Leiden across runs — pick one and pin it.

**Pitfall — `dims` must be consistent.** Whatever `DIMS_CHOSEN` you used in
Step 2's PCA must be reused in `FindNeighbors`, in this `RunUMAP`, and
anywhere else you reference dimensions. Inconsistent `dims` is the #1
source of "the UMAP looks different from the clusters" confusion.

---

## Step 6 — Assess mixing (did integration actually help?)

Integration succeeds when cells from different samples **interleave**
within shared cell types yet **distinct cell types stay separate**. Read
the post-integration UMAP both by sample (should be mixed, contrast with
Step 3) and by cluster (should show clean populations).

```r
# By sample (the diagnostic for batch correction)
p_post_sample <- DimPlot(obj, reduction = "umap.integrated",
                         group.by = "sample", pt.size = 0.4) +
  ggtitle("After integration — colored by sample (should be mixed)") +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4)) +
  coord_fixed()
for (k in seq_along(p_post_sample$layers)) {
  if (inherits(p_post_sample$layers[[k]]$geom, "GeomPoint")) {
    p_post_sample$layers[[k]]$aes_params$alpha <- 0.6
  }
}
ggsave("umap_integrated_by_sample.png", p_post_sample,
       width = 7, height = 6, units = "in", dpi = 120, bg = "white")

# By cluster (the biological readout)
p_post_clust <- DimPlot(obj, reduction = "umap.integrated",
                        group.by = "integrated_clusters",
                        label = TRUE, repel = TRUE, pt.size = 0.4) +
  ggtitle(sprintf("After integration — %d clusters (res = 1.0, dims = 1:%d)",
                  length(unique(obj$integrated_clusters)), DIMS_CHOSEN)) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4)) +
  NoLegend() +              # AFTER theme_cowplot — cowplot 1.2.0 overrides legend.position
  coord_fixed()
for (k in seq_along(p_post_clust$layers)) {
  if (inherits(p_post_clust$layers[[k]]$geom, "GeomPoint")) {
    p_post_clust$layers[[k]]$aes_params$alpha <- 0.6
  }
}
ggsave("umap_integrated_by_cluster.png", p_post_clust,
       width = 7, height = 6, units = "in", dpi = 120, bg = "white")
```

**Report:** is the by-sample UMAP well-mixed (good) or do samples still
form sample-specific islands within shared cell types
(under-integrated)? Do biologically distinct cell types collapse into one
cluster (over-integrated)? Cells per sample per cluster:

```r
print(table(cluster = obj$integrated_clusters, sample = obj$sample))
```

If a cluster is dominated by one sample (>80% from one sample when the
cohort is balanced), that's either a sample-specific population (real
biology) or under-integration.

For mixing-vs-collapse diagnostics + patches when integration is
under/over-applied, read `references/over_correction_diagnostics.md`.

---

## Step 7 — JoinLayers before any DE / marker calling

Integration leaves counts split across per-sample layers. Before
`FindMarkers` / `FindAllMarkers` / pseudobulk DE, rejoin the layers into a
single matrix:

```r
# Rejoin counts / data layers across samples. Works on the assay or the object.
obj <- JoinLayers(obj)
```

DE runs on the **joined RNA expression**, never on the integrated
embedding (the embedding is for neighbors / UMAP only). For per-cluster
marker discovery see `seurat-scrna-v2`'s Step 7 (`FindAllMarkers` with
Wilcoxon defaults); for sample-level / condition-level comparisons use
**pseudobulk** (aggregate to per-sample-per-cluster counts, then `DESeq2`
/ `limma-voom`) — NOT direct per-cell DE which violates the test's
assumption that cells are independent samples.

---

## Step 8 — Save the integrated object

```r
saveRDS(obj, "seurat_integrated.rds")

# Verify the write
fi <- file.info("seurat_integrated.rds")
cat(sprintf("Wrote seurat_integrated.rds (%.1f MB)\n", fi$size / 1e6))
stopifnot(file.exists("seurat_integrated.rds"))
```

The saved `.rds` carries: the merged + joined RNA assay (counts + data +
scale.data), the `pca` reduction (pre-integration), the integrated
reduction (`integrated.cca` / `integrated.rpca` / `harmony` /
`integrated.scvi` / `integrated.dr` depending on the method/path used),
both UMAPs (`umap.unintegrated` + `umap.integrated`), the integrated SNN
graph, the `integrated_clusters` metadata column, and the 2000 HVGs. A
follow-up session can pick up at Step 7+ without re-running 1–6.

---

## Batch variant — use INSTEAD of Steps 1–8 when invoked with args="batch"

Branch on `$ARGUMENTS` at the top of the body. In batch mode:

- Skip all three diagnostic UMAPs and the elbow figure — the orchestrator
  produces the cross-batch comparative figures instead.
- Skip the per-step "Report" footers — the orchestrator does the rollup.
- Still produce the canonical `seurat_integrated.rds` (the orchestrator
  consumes it).
- Print ONE final summary line of the form
  `"batch ok | N samples | M cells | K clusters | <method>"` so the
  orchestrator can grep it from the log.

```r
# Batch path — consolidated, silent, one artifact, one summary line.
# Replace sample_dirs / samples / method / DIMS_CHOSEN with the values the
# orchestrator passes in (or with literals if calling solo with args="batch").

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
stopifnot(packageVersion("Seurat") >= "5.0.0")

# --- Load + merge + split ---
# (use ONE of the loader branches from Step 1; same paths the agent already has.)
objs <- lapply(names(sample_dirs), function(s) {
  CreateSeuratObject(Read10X(sample_dirs[[s]]),
                     project = s, min.cells = 3, min.features = 200)
})
obj <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(sample_dirs))
obj$sample <- obj$orig.ident
obj <- JoinLayers(obj)                              # merge() already splits v5; consolidate before re-splitting
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)

# --- Preprocess ---
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

# --- Integrate (default CCA; swap method as needed) ---
DIMS_CHOSEN <- 30
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)

# --- Cluster + UMAP on the integrated reduction (no plot save) ---
obj <- FindNeighbors(obj, reduction = "integrated.cca",
                     dims = 1:DIMS_CHOSEN, verbose = FALSE)
obj <- FindClusters(obj, resolution = 1.0,
                    cluster.name = "integrated_clusters", verbose = FALSE)
obj <- RunUMAP(obj, reduction = "integrated.cca",
               dims = 1:DIMS_CHOSEN, reduction.name = "umap.integrated",
               verbose = FALSE)

# --- JoinLayers + save the canonical artifact ---
obj <- JoinLayers(obj)
saveRDS(obj, "seurat_integrated.rds")

cat(sprintf("batch ok | %d samples | %d cells | %d clusters | CCAIntegration\n",
            length(unique(obj$sample)), ncol(obj),
            length(unique(obj$integrated_clusters))))
```

What batch mode keeps vs drops:
- KEEPS: the canonical processed object (`seurat_integrated.rds`) and the
  one summary line.
- DROPS: `ggsave` calls, per-step `print()` reports, the elbow figure,
  the cross-sample tabulation.

---

## Final response checklist

Summarize, in this order:

- samples loaded and integrated (count + names, cells per sample)
- normalization path (LogNormalize vs SCT) and integration method
  (`CCAIntegration` / `RPCAIntegration` / `HarmonyIntegration` /
  `JointPCAIntegration` / `scVIIntegration` / `FastMNNIntegration`)
- source reduction (`orig.reduction`) and integrated reduction
  (`new.reduction`) names
- `dims` used (must match across PCA / IntegrateLayers / FindNeighbors /
  RunUMAP)
- clustering algorithm (Louvain default; Leiden if `algorithm = 4`),
  resolution, and resulting cluster count
- cells-per-sample-per-cluster contingency (was a cluster sample-skewed?)
- visual assessment: pre-vs-post by-sample UMAP mixing; any populations
  collapsed by over-integration
- whether `JoinLayers` was called (required before DE)
- files written (figures + `seurat_integrated.rds`) and their location
- caveats: per-sample QC done? signal-vs-nuisance split correct?
  over-integration risk? Leiden availability?

---

## See also

- `seurat-scrna-v2` — single-sample QC → clustering → markers. Run this
  PER sample before integrating, and consult its Step 7 for the
  marker-calling code that runs on the joined RNA assay after `JoinLayers`.
- `seurat-reference-mapping` — once you have an integrated, annotated
  reference, this projects a new query dataset onto it via
  `FindTransferAnchors` + `TransferData` + `MapQuery`.
- `harmony-integration-scanpy` — same Harmony correction in a Python /
  scanpy session.
- `scvi-integration` — atlas-scale deep-generative integration; preferred
  over `scVIIntegration` here when the user explicitly wants scVI / scANVI.
