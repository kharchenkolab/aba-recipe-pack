---
name: seurat-cite-seq
description: CITE-seq (paired RNA + ADT surface-protein) analysis in R/Seurat v5 — load a 10x multimodal bundle, build a Seurat object with an RNA assay AND an ADT assay, CLR-normalize ADT per-feature (margin=2), run the standard RNA QC/PCA/UMAP/clustering pipeline, and characterize clusters with biaxial protein plots (FeatureScatter) + cross-modality FeaturePlots (`adt_CD3` vs `rna_CD3E`). For joint clustering driven by BOTH modalities, step up to seurat-wnn-multimodal.
when_to_use: A 10x Genomics CITE-seq run with Gene Expression + Antibody Capture libraries (a `filtered_feature_bc_matrix` directory or `.h5` whose feature table splits into both modalities), and you want per-cell protein levels alongside the RNA clustering. Use THIS when the user wants the per-cell ADT workup with biaxial flow-style scatters and ADT overlays on the RNA UMAP. For joint RNA+ADT clustering (cell-specific modality weighting), switch to `seurat-wnn-multimodal`. For mapping a CITE-seq query onto a pre-built multimodal reference, use `seurat-multimodal-reference-mapping`.
invocation: interactive+batch
requires_tools: [run_r]
capabilities_needed: [Seurat]
keywords: [Seurat, CITE-seq, ADT, antibody derived tags, surface protein, multimodal, Antibody Capture, Gene Expression, CLR, centered log ratio, NormalizeData, FeatureScatter, biaxial, adt_, rna_, DefaultAssay, CreateAssayObject, R, v5]
produces: [adt_qc_violins.png, umap_rna_clusters.png, adt_biaxial.png, adt_featureplot.png, adt_vs_rna.png, cite_processed.rds, cite_processed.lstar.zarr]
domain: genomics
source: "Seurat CITE-seq multimodal vignette (Satija Lab) — satijalab.org/seurat/articles/multimodal_vignette (Seurat 5.5.0); Stoeckius et al. 2017 (CITE-seq, Nat. Methods)"
---

# CITE-seq (paired RNA + ADT) analysis with R/Seurat (v5)

A CITE-seq library is one 10x run that captures **two feature classes per cell**:
gene-expression mRNAs (Gene Expression) and antibody-derived tags reading surface
proteins (Antibody Capture). Seurat stores these as **two assays on one
SeuratObject** — `RNA` and `ADT` — and lets you switch between them with
`DefaultAssay()` or query either via the `rna_<gene>` / `adt_<protein>` feature
keys. This recipe runs the standard RNA workflow (QC → norm → PCA → UMAP →
clusters) and adds the **ADT-specific** steps: CLR normalization per feature,
biaxial protein plots, and protein overlays on the RNA UMAP.

This recipe does NOT do joint RNA+ADT clustering — that's WNN
(`seurat-wnn-multimodal`). The clustering here is driven by RNA only; ADT is
used to *characterize* and *validate* the RNA clusters in protein space.

Pin **Seurat v5** (the vignette tests on `Seurat 5.5.0`). The `[[<-` assay
syntax and `CreateAssayObject`/`CreateAssay5Object` patterns below are v5; on
v4 the API works but layers differ. Confirm with `packageVersion("Seurat")`.

## Bundled references — load on demand

- `references/adt_normalization.md` — CLR `margin=1` vs `margin=2` semantics,
  panel-size effects, DSB alternative for isotype-equipped panels.
- `references/adt_qc_and_panel.md` — panel-suffix cleanup, isotype controls,
  per-protein QC gates, hashtag/HTO confounders.
- `references/biaxial_visualization.md` — FeatureScatter flow-cytometry-style
  biaxials (CLR vs raw), `adt_`/`rna_` feature-key conventions, FeaturePlot
  overlays on the RNA UMAP.
- `references/figure_style.md` — Seurat-collection figure-style cheatsheet
  (palettes, theme_cowplot, alpha-poke, ggsave conventions).

## Install

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("Seurat",    quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2",   quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",     quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("cowplot",   quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
# Pin v5 — the CLR margin=2 + adt_/rna_ key behavior is the v5 contract.
stopifnot(packageVersion("Seurat") >= "5.0.0")
```

Attach once per session (NOT per code block):

```r
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(cowplot); library(patchwork)
})
```

## Decisions to surface up front

Tell the user these are the analysis-defining decisions:

1. **Species / MT prefix** — `^MT-` (human) vs `^mt-` (mouse) for `percent.mt`.
   Same as any scRNA recipe — gets %MT to 0 if wrong, and dying cells survive QC.
2. **ADT panel size** — small panels (10–30 markers) cluster better on protein
   than large ones; very small panels (<10) should skip ADT-only clustering
   entirely and use ADT for marker overlay only. See
   `references/adt_qc_and_panel.md` for the per-protein QC + panel-size
   selection criteria.
3. **ADT normalization choice** — CLR with `margin=2` (per-feature, across
   cells) is the Seurat default and recommended for most panels. CLR with
   `margin=1` (per-cell, across features), DSB (when isotypes + empty droplets
   are available) — see `references/adt_normalization.md` for the
   decision tree.
4. **Cross-modality feature naming** — Seurat exposes `adt_<protein>` and
   `rna_<gene>` keys that bypass `DefaultAssay()`. Prefer these in plots; they
   make the assay explicit at the call site and prevent silent assay swaps.
5. **Mode** — INTERACTIVE (default) when this is the primary analysis;
   BATCH when invoked as one of many sub-runs in a larger plan (e.g. one
   CITE-seq run per donor). The orchestrator passes `args="batch"`; the
   agent declares the mode in `present_plan` so the user sees which one is
   in effect.

Figures to show as the analysis proceeds:
- `adt_qc_violins.png`
- `umap_rna_clusters.png`
- `adt_biaxial.png`
- `adt_featureplot.png`
- `adt_vs_rna.png`

---

## Step 1 — Load the multimodal 10x bundle and build the object

`Read10X()` on a CITE-seq `filtered_feature_bc_matrix` directory returns a
**named list** with one matrix per feature class. The class names match the
`features.tsv` `feature_type` column — typically `"Gene Expression"` and
`"Antibody Capture"`.

```r
# Multimodal 10x directory: features.tsv carries a feature_type column,
# Read10X returns a list named by that column.
counts10x <- Read10X(data.dir = "/path/to/filtered_feature_bc_matrix")

# Inspect the modality split — names() MUST contain both modalities.
print(names(counts10x))
stopifnot("Gene Expression"  %in% names(counts10x))
stopifnot("Antibody Capture" %in% names(counts10x))

rna_counts <- counts10x[["Gene Expression"]]
adt_counts <- counts10x[["Antibody Capture"]]

# Clean noisy panel suffixes (TotalSeq-B style "_TotalSeqB" / "_control_TotalSeqB")
# See references/adt_qc_and_panel.md for the full regex menu + symbol audit.
rownames(adt_counts) <- gsub("_[A-Za-z_]*TotalSeq[A-C]?$", "", rownames(adt_counts))
```

### Variant: 10x .h5 input — use INSTEAD of the directory block when the bundle is a single HDF5

```r
counts10x <- Read10X_h5("/path/to/filtered_feature_bc_matrix.h5")
# Same list shape as the directory loader; same downstream split.
rna_counts <- counts10x[["Gene Expression"]]
adt_counts <- counts10x[["Antibody Capture"]]
```

Now build a single Seurat object with both assays. RNA is the primary assay;
ADT is added as a second assay restricted to the cells the RNA object kept
(the per-modality cell sets differ because RNA's `min.features=200` filter
drops empty droplets that ADT may still see):

```r
obj <- CreateSeuratObject(
  counts       = rna_counts,
  project      = "cite_seq",
  min.cells    = 3,
  min.features = 200
)

# Restrict ADT to RNA-kept cells; CreateAssayObject builds a v5-compatible assay.
# (Use CreateAssay5Object if you explicitly want a v5 layered assay.)
obj[["ADT"]] <- CreateAssayObject(counts = adt_counts[, colnames(obj)])
```

**Report:** number of RNA cells / RNA genes / ADT proteins / ADT cells (must
match RNA cells after the subset). Confirm `Assays(obj)` returns
`c("RNA", "ADT")`.

```r
cat(sprintf("RNA: %d genes x %d cells | ADT: %d proteins x %d cells\n",
            nrow(obj[["RNA"]]), ncol(obj[["RNA"]]),
            nrow(obj[["ADT"]]), ncol(obj[["ADT"]])))
print(Assays(obj))
stopifnot(identical(Assays(obj), c("RNA", "ADT")))
```

For panel-suffix variants beyond TotalSeq-B (BioLegend `_C`, hashtag mixing),
isotype-control detection, and per-protein QC gates, read
`references/adt_qc_and_panel.md`.

---

## Step 2 — Run the RNA workflow (QC → norm → HVG → PCA → UMAP → clusters)

Same Seurat single-sample workflow as `seurat-scrna-v2`. The full QC plotting
+ marker discussion lives there; here we run the canonical chain with sensible
defaults so the recipe stays focused on the *protein* steps.

```r
DefaultAssay(obj) <- "RNA"

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")   # human; ^mt- for mouse
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)

cat(sprintf("RNA clusters: %d (n=%d cells)\n",
            length(levels(Idents(obj))), ncol(obj)))
```

UMAP colored by RNA cluster — same DimPlot styling as the canonical Seurat
recipe (label on the embedding, no side legend, alpha-poke to 0.6). For the
full figure-style cheatsheet (palettes, alpha-poke, theme_cowplot ordering
quirks), read `references/figure_style.md`.

```r
p_umap <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE,
                  pt.size = 0.4) +
  ggtitle(sprintf("RNA clusters · res=0.8 · n=%d cells · %d clusters",
                  ncol(obj), length(levels(Idents(obj))))) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank()) +
  NoLegend() +              # AFTER theme_cowplot — cowplot 1.2.0 overrides legend.position
  coord_fixed()
for (k in seq_along(p_umap$layers)) {
  if (inherits(p_umap$layers[[k]]$geom, "GeomPoint")) {
    p_umap$layers[[k]]$aes_params$alpha <- 0.6
  }
}
ggsave("umap_rna_clusters.png", p_umap,
       width = 7, height = 6.5, units = "in", dpi = 120, bg = "white")
```

**Report:** post-QC cells, RNA cluster count, cluster sizes (`table(Idents(obj))`).
For deeper QC / HVG / PC-count tuning, see `seurat-scrna-v2`.

---

## Step 3 — CLR-normalize the ADT assay

The ADT count distribution is unlike RNA: each protein has a low-background
peak (non-specific binding) and a specific-binding peak. Seurat's recommended
normalization is the **centered log-ratio (CLR) with `margin = 2`** — that
is, the CLR is computed *across cells, per feature*, which is what stabilizes
the per-protein bimodality.

```r
# margin = 2 → CLR across cells, per feature (Seurat's recommended default).
# margin = 1 → CLR across features, per cell (older convention; legacy).
# assay = "ADT" routes the call to the ADT assay specifically.
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, assay = "ADT")

# ScaleData on ADT — used for heatmaps and any RunPCA on ADT downstream.
obj <- ScaleData(obj, assay = "ADT")
```

**Report:** number of ADT features, fraction of cells with any ADT counts,
and a quick `summary()` of normalized values for two reference proteins so
the user sees the CLR landed sensibly (typical range ~0 to ~5).

```r
adt_norm <- GetAssayData(obj, assay = "ADT", layer = "data")
cat(sprintf("ADT: %d proteins x %d cells normalized (CLR margin=2)\n",
            nrow(adt_norm), ncol(adt_norm)))
# Spot-check two proteins — adapt names to your panel.
example_proteins <- intersect(c("CD3", "CD19", "CD4", "CD8"), rownames(adt_norm))
if (length(example_proteins) > 0) {
  print(t(apply(adt_norm[example_proteins, , drop = FALSE], 1, summary)))
}
```

For `margin=1` vs `margin=2` semantics, the DSB alternative (when isotypes +
empty droplets are available), and panel-size effects on normalization
choice, read `references/adt_normalization.md`.

---

## Step 4 — ADT QC plots (per-protein violins on the RNA clusters)

The first protein view is a stack of per-protein VlnPlots grouped by RNA
cluster — does the staining track the clusters? Use the `adt_<protein>`
feature key so the assay is explicit at the call site:

```r
# Replace `panel_show` with the 6-12 lineage markers in YOUR panel.
panel_show <- c("CD3", "CD4", "CD8", "CD19", "CD56", "CD14")
panel_show <- panel_show[panel_show %in% rownames(obj[["ADT"]])]

p_vln <- VlnPlot(obj,
                 features = paste0("adt_", panel_show),
                 pt.size  = 0,
                 ncol     = 3) &
  theme_cowplot() &
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none",
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor   = element_blank())

n_rows <- ceiling(length(panel_show) / 3)
ggsave("adt_qc_violins.png", p_vln,
       width = 12, height = max(4, 3 * n_rows), units = "in", dpi = 120, bg = "white")
```

**Assess and report:** which proteins cleanly mark single RNA clusters
(good — the lineage was captured), which proteins are diffusely positive
across many clusters (background staining or pan-marker), which proteins
show no signal (failed antibody, missing from cells, or panel-mapping
issue). For the diagnosis table mapping each pattern to a likely cause
+ remediation, see `references/adt_qc_and_panel.md`.

---

## Step 5 — Biaxial ADT scatters (flow-cytometry style)

Two-protein scatters on CLR-normalized ADT — the protein-domain analog of a
flow gate. `FeatureScatter` accepts the `adt_<protein>` key directly. For
the full FeatureScatter / FeaturePlot / `slot="counts"` reference (CLR vs
raw axes, `adt_`/`rna_` key conventions, panel adaptation), read
`references/biaxial_visualization.md`.

```r
# Three canonical PBMC biaxials — adapt to your panel:
b1 <- FeatureScatter(obj, feature1 = "adt_CD4",  feature2 = "adt_CD8")
b2 <- FeatureScatter(obj, feature1 = "adt_CD3",  feature2 = "adt_CD19")
b3 <- FeatureScatter(obj, feature1 = "adt_CD14", feature2 = "adt_CD16")

p_biax <- (b1 | b2 | b3) &
  theme_cowplot() &
  theme(plot.title = element_text(size = 10, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank()) &
  NoLegend()              # AFTER theme_cowplot, broadcast via & to all sub-panels

# Alpha-poke the points (FeatureScatter doesn't expose alpha)
for (i in seq_along(p_biax)) {
  pl <- p_biax[[i]]
  if (!is.null(pl$layers)) {
    for (k in seq_along(pl$layers)) {
      if (inherits(pl$layers[[k]]$geom, "GeomPoint")) {
        pl$layers[[k]]$aes_params$alpha <- 0.6
      }
    }
    p_biax[[i]] <- pl
  }
}

ggsave("adt_biaxial.png", p_biax,
       width = 15, height = 5, units = "in", dpi = 120, bg = "white")
```

**Assess and report:** are the canonical quadrants populated as expected
(CD4+CD8-, CD8+CD4-, CD3+CD19-, …)? Cells in the wrong quadrant relative to
their RNA cluster are the suspicious cells worth a second look.

---

## Step 6 — ADT FeaturePlots on the RNA UMAP + cross-modality comparison

Overlay each ADT marker on the RNA UMAP (using the `adt_` key) so the
protein gradients sit in the same embedding as the RNA clusters:

```r
panel_show <- c("CD3", "CD4", "CD8", "CD19", "CD56", "CD14")   # adapt to panel
panel_show <- panel_show[panel_show %in% rownames(obj[["ADT"]])]

p_adt_fp <- FeaturePlot(obj,
                        features = paste0("adt_", panel_show),
                        reduction = "umap",
                        order = TRUE,
                        pt.size = 0.3,
                        ncol = 3) &
  scale_colour_gradient(low = "grey85", high = "#b2182b") &
  theme_cowplot() &
  theme(plot.title = element_text(size = 11, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"))

for (i in seq_along(p_adt_fp)) {
  pl <- p_adt_fp[[i]]
  if (!is.null(pl$layers)) {
    for (k in seq_along(pl$layers)) {
      if (inherits(pl$layers[[k]]$geom, "GeomPoint")) {
        pl$layers[[k]]$aes_params$alpha <- 0.6
      }
    }
    p_adt_fp[[i]] <- pl
  }
}

n_rows <- ceiling(length(panel_show) / 3)
ggsave("adt_featureplot.png", p_adt_fp,
       width = 13, height = max(4, 4 * n_rows),
       units = "in", dpi = 120, bg = "white")
```

Then pair each protein with its RNA-message gene side-by-side — this is the
classic CITE-seq sanity check, looking for cells where protein > RNA (which
is biology — translation lag, surface stable) and cells where RNA > protein
(also biology — newly induced transcripts before staining):

```r
# Pair protein with its gene; gene names sometimes differ from protein labels
# (e.g. ADT 'CD3' <-> RNA 'CD3E', ADT 'CD19' <-> RNA 'MS4A1' for B cells).
# Adapt this list to YOUR panel — the protein-to-gene map is panel-specific.
pairs <- list(
  c("adt_CD3",  "rna_CD3E"),
  c("adt_CD4",  "rna_CD4"),
  c("adt_CD8",  "rna_CD8A"),
  c("adt_CD19", "rna_CD19")
)
# Drop pairs where either side is missing
have_feat <- function(f) {
  pfx <- substr(f, 1, 4); name <- substr(f, 5, nchar(f))
  if (pfx == "adt_") name %in% rownames(obj[["ADT"]]) else name %in% rownames(obj[["RNA"]])
}
pairs <- Filter(function(p) all(vapply(p, have_feat, logical(1))), pairs)
features_flat <- unlist(pairs)

p_avr <- FeaturePlot(obj, features = features_flat,
                     reduction = "umap", order = TRUE, pt.size = 0.3,
                     ncol = 2) &
  scale_colour_gradient(low = "grey85", high = "#b2182b") &
  theme_cowplot() &
  theme(plot.title = element_text(size = 11, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"))

n_rows <- ceiling(length(features_flat) / 2)
ggsave("adt_vs_rna.png", p_avr,
       width = 10, height = max(4, 3.5 * n_rows),
       units = "in", dpi = 120, bg = "white")
```

**Assess and report:** for each canonical lineage marker, does the ADT
gradient track the RNA gradient on the UMAP? Mismatches concentrated in one
cluster argue for re-examining either the antibody (dropout in that
population) or the RNA gene mapping (wrong symbol, e.g. CD19 surface protein
vs. MS4A1 transcript for B cells). For the protein-to-gene symbol audit and
how to read protein > RNA vs RNA > protein discordance, see
`references/biaxial_visualization.md`.

---

## Step 7 — Step up to joint clustering (when to leave this recipe)

This recipe clusters on RNA and overlays ADT. When the panel is large and
informative (>20 markers, lineage-relevant), **joint clustering** that weighs
both modalities cell-by-cell typically beats the RNA-only clustering at
resolving fine cell types. Symptoms that argue for the step up:

- ADT VlnPlots (Step 4) show lineage markers crossing several RNA clusters,
  but with sharper boundaries than the RNA itself.
- ADT biaxials (Step 5) reveal sub-populations within one RNA cluster.
- You're working with bone marrow, blood, or any tissue where the
  hematopoietic lineages have stronger surface markers than transcripts.

For joint clustering, save the object and switch:

```r
Skill(skill = "seurat-wnn-multimodal")
```

Pass the same object — WNN starts from the same per-modality reductions.

---

## Step 8 — Save the processed object

```r
saveRDS(obj, file = "cite_processed.rds")
cat(sprintf("Wrote cite_processed.rds (%.1f MB)\n",
            file.info("cite_processed.rds")$size / 1e6))
```

The saved `.rds` carries both assays (RNA + ADT), the RNA reductions
(`pca`, `umap`), the RNA neighbor graphs, the cluster identities, and the
CLR-normalized ADT data layer.

---

## Batch variant — use INSTEAD of Steps 1–7 when invoked with args="batch"

Branch on `$ARGUMENTS == "batch"` at the top of the body. In batch mode:

- Skip all per-step figures (no `ggsave` calls).
- Skip the per-step "Assess and report" footers — the orchestrator does the
  rollup.
- Still save the canonical `cite_processed.rds` so downstream orchestrator
  steps can pick it up.
- Print ONE final summary line.

```r
suppressPackageStartupMessages({ library(Seurat) })

counts10x <- Read10X(data.dir = "/path/to/filtered_feature_bc_matrix")
rna_counts <- counts10x[["Gene Expression"]]
adt_counts <- counts10x[["Antibody Capture"]]
rownames(adt_counts) <- gsub("_[A-Za-z_]*TotalSeq[A-C]?$", "", rownames(adt_counts))

obj <- CreateSeuratObject(counts = rna_counts, project = "cite_seq",
                          min.cells = 3, min.features = 200)
obj[["ADT"]] <- CreateAssayObject(counts = adt_counts[, colnames(obj)])

DefaultAssay(obj) <- "RNA"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)
obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE)

obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, assay = "ADT")
obj <- ScaleData(obj, assay = "ADT")

saveRDS(obj, file = "cite_processed.rds")
cat(sprintf("batch ok | %d cells | %d clusters | %d ADT proteins | %.1f MB\n",
            ncol(obj),
            length(levels(Idents(obj))),
            nrow(obj[["ADT"]]),
            file.info("cite_processed.rds")$size / 1e6))
```

---

## Final response checklist

Summarize:
- input format (10x directory vs .h5) and the modality split confirmed (RNA / ADT cells, genes, proteins)
- ADT panel size and which proteins were retained after the symbol cleanup
- RNA QC thresholds applied (nFeature low/high, percent.mt cutoff) and post-QC cell count
- RNA clustering settings (PC count, resolution) and cluster count
- ADT normalization (CLR `margin=2`) — and any DSB substitution if isotypes were available
- which ADT proteins cleanly track RNA clusters (lineage markers good signal)
- which proteins look diffuse / failed (background, dropped antibody)
- ADT-vs-RNA pair comparisons — concordance of each lineage marker on the UMAP
- figures shown to the user (adt_qc_violins, umap_rna_clusters, adt_biaxial, adt_featureplot, adt_vs_rna)
- saved files (cite_processed.rds)
- caveats: ADT panel may be too small for joint clustering (step up to `seurat-wnn-multimodal`), no doublet check run, hashtag/HTO de-multiplexing assumed already done upstream, panel-specific protein/gene name mismatches require human judgment

---

## See also

- `seurat-wnn-multimodal` — joint clustering driven by BOTH modalities with
  cell-specific weights; switch when the panel is large and informative and
  the RNA-only clusters look under-resolved.
- `seurat-multimodal-reference-mapping` — map this CITE-seq query onto a
  pre-built multimodal reference (e.g. PBMC 162k) for cell-type labels at
  multiple granularities.
- `seurat-scrna-v2` — the RNA-only single-sample workflow this recipe shares
  its QC / PCA / UMAP steps with.
- `bp-cite-seq` — best-practice surface-protein chapter (sample-aware MAD
  QC, DSB normalization, isotype controls); read when the panel is large
  and the dataset is from a multi-donor study.

## Offer an interactive view

Write a viewer-optimized store DIRECTLY from the live Seurat object with lstar (pure R,
highest fidelity — do NOT route through `.h5ad`) and **proactively offer to open it** (a
required part of delivering the result):
```r
DefaultAssay(obj) <- "RNA"        # RNA expression view (ADT was the last active assay)
d <- lstar::read_seurat(obj)
lstar::lstar_write_viewer(d, "cite_processed.lstar.zarr")   # precomputes DE / HVGs /
                                                            # cell-major counts (optimized)
```
Then call `open_viewer(file_path="cite_processed.lstar.zarr")` and present the returned link
so the user can explore the RNA clusters + metadata on the UMAP in pagoda3 — it opens
instantly (pre-optimized, no on-launch conversion). If `open_viewer` returns `ok:false`,
relay the error rather than a dead link. Format / sharing → **`scrna-viewing-and-interchange`**.
