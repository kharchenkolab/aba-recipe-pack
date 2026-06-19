---
name: seurat-wnn-multimodal
description: Weighted Nearest Neighbor (WNN) multimodal clustering in R/Seurat v5 — process each modality independently (RNA + ADT, or RNA + ATAC) to get its own dimensional reduction, then learn cell-specific modality weights with FindMultiModalNeighbors() and cluster on the joint wsnn graph. Produces one shared UMAP where the modalities' contributions are weighted per cell, plus a per-cell `RNA.weight` metadata column revealing which cells are RNA-driven vs protein/chromatin-driven.
when_to_use: A multimodal Seurat object whose two assays carry complementary biological signal — typically CITE-seq (RNA + ADT, panel ≥20 markers) or 10x Multiome (RNA + ATAC) — and the user wants ONE clustering and ONE UMAP that respects both modalities, with cell-specific weights instead of a fixed blend. Use THIS when the RNA-only clusters under-resolve populations the protein/chromatin signal would split, or when the user explicitly asks for WNN / weighted nearest neighbors / joint clustering. For ADT-as-overlay-on-RNA (no joint clustering), use `seurat-cite-seq`. For RNA + ATAC specifically (Signac LSI on the ATAC side), see `seurat-rna-atac-integration`.
invocation: interactive+batch
requires_tools: [run_r]
capabilities_needed: [Seurat]
keywords: [Seurat, WNN, weighted nearest neighbor, multimodal clustering, FindMultiModalNeighbors, wsnn, wknn, weighted.nn, modality weight, RNA.weight, joint clustering, CITE-seq, RNA ADT, RNA ATAC, multiome, bmcite, R, v5]
produces: [umap_wnn_clusters.png, umap_rna_only.png, umap_adt_only.png, modality_weights.png, wnn_processed.rds]
domain: genomics
source: "Seurat Weighted Nearest Neighbor (WNN) vignette (Satija Lab) — satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis (Seurat 5.5.0, bmcite dataset); Hao et al. 2021, Cell 184:3573 (the WNN paper)"
---

# Weighted Nearest Neighbor (WNN) joint clustering with R/Seurat (v5)

WNN takes a Seurat object with **two (or more) modalities already processed
independently** — each with its own dimensional reduction — and learns a
**per-cell weighting** between them. The weights are derived from how well
each modality's local neighborhood predicts the cell's identity: cells whose
protein neighborhood is sharper than their RNA neighborhood (e.g. T cells in
a CITE-seq dataset) get up-weighted on ADT; cells where RNA carries the
identity (e.g. early progenitors) get up-weighted on RNA. The output is one
shared neighbor graph (`weighted.nn`), one shared SNN graph for clustering
(`wsnn`), one joint UMAP, and a per-cell `RNA.weight` column (1.0 = pure
RNA, 0.0 = pure ADT/ATAC).

This recipe assumes you have already produced **per-modality reductions** —
either from `seurat-cite-seq` (RNA `pca` + ADT processing) or from
`seurat-rna-atac-integration` (RNA `pca` + ATAC `lsi`). If you don't, do
that first; WNN doesn't run from raw counts.

Pin **Seurat v5** (the vignette tests on `Seurat 5.5.0`). The
`FindMultiModalNeighbors()` API and the `weighted.nn` / `wknn` / `wsnn`
artifact names are stable v3.2+ but the recipe uses the v5 default modality
storage (assays-as-layers, `[[<-` for new assays).

## Bundled references — load on demand

- `references/wnn_internals.md` — FindMultiModalNeighbors arguments,
  artifact placement (`@neighbors` vs `@graphs`), the
  `modality.weight.name` 2-element-vector gotcha (single-string silent
  fallback to `<DefaultAssay>.weight`), k.nn and prune.SNN tuning.
- `references/modality_weight_diagnostics.md` — interpreting the
  per-cell `RNA.weight` column, the sorted VlnPlot pattern, what a
  bimodal vs flat weight distribution implies, expected ranges per
  tissue.
- `references/figure_style.md` — Seurat-collection figure-style
  cheatsheet (palettes, theme_cowplot, alpha-poke, ggsave
  conventions).

## Install

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("Seurat",      quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2",     quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",       quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("cowplot",     quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("patchwork",   quietly = TRUE)) install.packages("patchwork")
stopifnot(packageVersion("Seurat") >= "5.0.0")
```

Attach once per session:

```r
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(cowplot); library(patchwork)
})
```

## Decisions to surface up front

Tell the user these are the analysis-defining decisions:

1. **Which reductions to combine** — typically `pca` (RNA) and `apca`
   (ADT-PCA), or `pca` (RNA) and `lsi` (ATAC). The reductions must already
   exist on the object. Name them consistently — WNN matches by reduction
   name in `reduction.list`.
2. **Per-modality dimensions** — `dims.list = list(1:30, 1:18)` is the
   vignette default for bmcite (RNA 30 PCs, ADT 18 PCs). The protein
   reduction needs fewer dims because the panel is small (25 antibodies in
   bmcite → 24 informative PCs max). Read the per-modality elbow before
   committing. See `references/wnn_internals.md` for the dims-tuning
   table.
3. **Clustering algorithm + resolution** — `FindClusters(graph.name="wsnn",
   algorithm=3, resolution=2)` is the vignette pattern. `algorithm=3` is the
   SLM (Smart Local Moving) clusterer; `resolution=2` is high because the
   joint graph resolves more populations. Tune resolution; **never** use
   `algorithm=1` (Louvain) without explicit reason — SLM handles the
   weighted graph better.
4. **Modality-weight interpretation** — the `RNA.weight` column reports
   the per-cell weighting (0 to 1). Populations with `RNA.weight` near 1.0
   are RNA-driven; near 0.0 are ADT/ATAC-driven. Surface this as a VlnPlot
   grouped by cluster — it's the diagnostic that tells you whether WNN was
   worth it. For the per-tissue expected-range table and how to read a
   bimodal vs flat weight distribution, see
   `references/modality_weight_diagnostics.md`.
5. **Mode** — INTERACTIVE (default) when this is the primary analysis;
   BATCH when invoked as one of many sub-runs in a larger plan (e.g. one
   multimodal sample per donor). The orchestrator passes `args="batch"`;
   the agent declares the mode in `present_plan`.

Figures to show as the analysis proceeds:
- `umap_wnn_clusters.png`
- `umap_rna_only.png`
- `umap_adt_only.png`
- `modality_weights.png`

---

## Step 1 — Verify the object has per-modality reductions

WNN starts from an existing Seurat object with both modalities already
processed. Confirm before calling `FindMultiModalNeighbors`:

```r
# obj is the multimodal SeuratObject loaded from upstream (cite_processed.rds,
# multiome_processed.rds, etc.). Verify both reductions exist.
obj <- readRDS("/path/to/cite_processed.rds")

# WNN expects: at minimum a `pca` reduction (RNA) AND a per-modality reduction
# for the second modality. CITE-seq convention: `apca` for ADT. Multiome
# convention: `lsi` for ATAC.
print(Reductions(obj))
stopifnot("pca" %in% Reductions(obj))
# Adapt the second-modality name to your case:
stopifnot(any(c("apca", "lsi") %in% Reductions(obj)))
```

If the second reduction is missing, run the per-modality steps first. For
CITE-seq:

```r
# ADT pipeline: CLR-normalize, scale, PCA into a NEW reduction named "apca"
# (so it doesn't overwrite the default "pca" from RNA).
DefaultAssay(obj) <- "ADT"
VariableFeatures(obj) <- rownames(obj[["ADT"]])              # use all ADT features
obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, assay = "ADT")
obj <- ScaleData(obj, assay = "ADT")
# ADT panels are small (often 20-50 markers); RunPCA defaults to npcs=50 +
# approx=TRUE (irlba). Both fail loudly on a 25-protein panel: "did not converge"
# and "too large % of singular values". Cap npcs at nfeatures-1 AND switch to
# the full-SVD path (approx=FALSE) which doesn't warn at high coverage.
n_adt    <- nrow(obj[["ADT"]])
adt_npcs <- min(50, n_adt - 1)
obj <- RunPCA(obj, assay = "ADT", reduction.name = "apca",
              npcs = adt_npcs, approx = FALSE, verbose = FALSE)
DefaultAssay(obj) <- "RNA"
```

For RNA + ATAC (10x Multiome), the ATAC reduction is LSI, not PCA — see
`seurat-rna-atac-integration` for the Signac-based `RunTFIDF` → `RunSVD`
pipeline that yields the `lsi` reduction.

**Report:** the two reductions WNN will combine and their dimensionalities.

```r
cat(sprintf("RNA reduction: %s (%d dims) | second modality: %s (%d dims)\n",
            "pca",
            ncol(Embeddings(obj, "pca")),
            if ("apca" %in% Reductions(obj)) "apca" else "lsi",
            ncol(Embeddings(obj,
              if ("apca" %in% Reductions(obj)) "apca" else "lsi"))))
```

---

## Step 2 — Run FindMultiModalNeighbors

Compute the per-cell modality weights and the joint neighbor graph. For
the full argument reference (`k.nn`, `knn.range`, `prune.SNN`,
`l2.norm`) and the `modality.weight.name` silent-fallback gotcha, read
`references/wnn_internals.md`.

```r
# reduction.list — names of the per-modality dimensional reductions
# dims.list      — per-modality dimensions to include (read from per-modality elbows)
# modality.weight.name — name of the metadata column that will store the RNA weight
#                        (Seurat infers the partner column from this: ADT.weight, etc.)
# CITE-seq case (RNA + ADT):
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list       = list("pca", "apca"),
  dims.list            = list(1:30, 1:18),
  # 2-element vector, not a single string — passing a single string triggers
  # Seurat's silent paste0(DefaultAssay, ".weight") fallback, which produces
  # SCT.weight (not RNA.weight) when DefaultAssay = SCT, breaking downstream
  # obj$RNA.weight references.
  modality.weight.name = c("RNA.weight", "ADT.weight"),
  verbose              = FALSE
)

# Outputs created on obj (verify):
# - obj[["weighted.nn"]] : multimodal neighbor object (used by RunUMAP nn.name=)
# - obj[["wknn"]]        : k-NN graph
# - obj[["wsnn"]]        : SNN graph (used by FindClusters graph.name=)
# - obj$RNA.weight        : per-cell RNA modality weight (0..1)
stopifnot("weighted.nn" %in% names(obj@neighbors))
stopifnot(all(c("wknn", "wsnn") %in% names(obj@graphs)))
stopifnot("RNA.weight" %in% colnames(obj@meta.data))
```

### Variant: RNA + ATAC (10x Multiome) — use INSTEAD of the canonical block above when the second modality is ATAC

```r
# ATAC reduction is `lsi`. Drop LSI component 1 — it captures sequencing depth,
# not biology (Signac convention).
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list       = list("pca", "lsi"),
  dims.list            = list(1:30, 2:30),
  modality.weight.name = c("RNA.weight", "ATAC.weight"),  # 2-element vec — see CITE-seq block above
  verbose              = FALSE
)
```

**Report:** range of `RNA.weight` (min, median, max) — if everything is near
0.5 the modalities are equally weighted (typical for balanced panels); a
broad spread (0.1–0.9) means WNN is doing real work.

```r
cat(sprintf("RNA.weight: min=%.2f  median=%.2f  max=%.2f\n",
            min(obj$RNA.weight), median(obj$RNA.weight), max(obj$RNA.weight)))
```

For the per-tissue expected `RNA.weight` ranges and how to read a
narrow vs broad spread, see `references/modality_weight_diagnostics.md`.

---

## Step 3 — Joint UMAP on the WNN graph

`RunUMAP` accepts a pre-computed neighbor object via `nn.name=` instead of
`dims=`. Use a non-default `reduction.name` so the WNN UMAP coexists with any
per-modality UMAP already on the object.

```r
obj <- RunUMAP(
  obj,
  nn.name        = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key  = "wnnUMAP_",
  verbose        = FALSE
)
```

Also compute the per-modality UMAPs side-by-side for comparison — same
parameters, different `reduction=` and `dims=`:

```r
obj <- RunUMAP(obj, reduction = "pca",  dims = 1:30, assay = "RNA",
               reduction.name = "rna.umap", reduction.key = "rnaUMAP_",
               verbose = FALSE)
obj <- RunUMAP(obj,
               reduction = if ("apca" %in% Reductions(obj)) "apca" else "lsi",
               dims = if ("apca" %in% Reductions(obj)) 1:18 else 2:30,
               assay = if ("apca" %in% Reductions(obj)) "ADT" else "ATAC",
               reduction.name = "adt.umap", reduction.key = "adtUMAP_",
               verbose = FALSE)
```

> **`RunUMAP` requires `nn.name` OR `dims=`, not both.** The first call
> (WNN) uses `nn.name`; the per-modality calls use `dims=`. Mixing them
> errors out — Seurat's argument-routing detects the conflict.

---

## Step 4 — Cluster on the joint SNN graph

`FindClusters` accepts a non-default graph via `graph.name=`. The WNN
graph is `wsnn`.

```r
# algorithm = 3 → SLM (Smart Local Moving). Recommended for the WNN graph.
#   algorithm = 1 (Louvain) is the default for RNA-only graphs but is
#   noticeably worse on the weighted graph.
# resolution is independent of RNA-only clustering — start at 2 (vignette
# default) and tune; the joint graph typically supports higher resolution.
obj <- FindClusters(
  obj,
  graph.name = "wsnn",
  algorithm  = 3,
  resolution = 2,
  verbose    = FALSE
)

cat(sprintf("WNN clusters: %d (n=%d cells)\n",
            length(levels(Idents(obj))), ncol(obj)))
print(table(Idents(obj)))
```

---

## Step 5 — Plot the joint UMAP, per-modality UMAPs, and the weight distribution

Three panels: joint WNN UMAP (the headline), then RNA-only and ADT/ATAC-only
UMAPs colored by the SAME (WNN) cluster labels — so the reader can see which
modality drives which cluster. Figure-style discipline (alpha-poke, cowplot
ordering, palette choices) is in `references/figure_style.md`.

```r
make_dimplot <- function(obj, reduction, title) {
  p <- DimPlot(obj, reduction = reduction, label = TRUE, repel = TRUE,
               pt.size = 0.4) +
    ggtitle(title) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          panel.grid.minor = element_blank()) +
    NoLegend() +              # AFTER theme_cowplot — cowplot 1.2.0 overrides legend.position
    coord_fixed()
  for (k in seq_along(p$layers)) {
    if (inherits(p$layers[[k]]$geom, "GeomPoint")) {
      p$layers[[k]]$aes_params$alpha <- 0.6
    }
  }
  p
}

p_wnn <- make_dimplot(obj, "wnn.umap", sprintf("WNN joint · %d clusters",
                      length(levels(Idents(obj)))))
p_rna <- make_dimplot(obj, "rna.umap", "RNA only (same WNN clusters)")
p_adt <- make_dimplot(obj, "adt.umap",
                      if ("apca" %in% Reductions(obj))
                        "ADT only (same WNN clusters)"
                      else
                        "ATAC only (same WNN clusters)")

ggsave("umap_wnn_clusters.png", p_wnn,
       width = 7, height = 6.5, units = "in", dpi = 120, bg = "white")
ggsave("umap_rna_only.png", p_rna,
       width = 7, height = 6.5, units = "in", dpi = 120, bg = "white")
ggsave("umap_adt_only.png", p_adt,
       width = 7, height = 6.5, units = "in", dpi = 120, bg = "white")
```

Modality-weight VlnPlot — sorted by median weight so the populations driven
by ADT (or ATAC) land at one end and the RNA-driven populations at the
other:

```r
p_w <- VlnPlot(obj, features = "RNA.weight",
               group.by = "seurat_clusters",
               sort = TRUE,
               pt.size = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "red", linewidth = 0.4) +
  ggtitle("Per-cell RNA modality weight (sorted by median)") +
  # Use ggplot's default discrete fill (hue_pal): the joint graph at resolution=2
  # routinely produces 15-25 clusters, which exceeds RColorBrewer Set2's 8-color
  # cap (n too large warning + recycled colors). hue_pal generates as many
  # distinct hues as needed.
  scale_fill_discrete(guide = "none") +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor   = element_blank()) +
  NoLegend()              # AFTER theme_cowplot — cowplot 1.2.0 overrides legend.position

ggsave("modality_weights.png", p_w,
       width = 10, height = 5, units = "in", dpi = 120, bg = "white")
```

**Assess and report:** the WNN cluster count vs. the RNA-only cluster count
(higher count typically validates the joint clustering); which clusters
sit cleanly on the RNA UMAP and which only resolve on the WNN UMAP (those
are the wins); the RNA.weight distribution — clusters above the dashed
0.5 line are RNA-driven, below are ADT/ATAC-driven, populations spread
across 0.5 are the ones WNN actually balanced. For the diagnosis table
and the per-tissue expected-range reference, read
`references/modality_weight_diagnostics.md`.

---

## Step 6 — Save the joint object

```r
saveRDS(obj, file = "wnn_processed.rds")
cat(sprintf("Wrote wnn_processed.rds (%.1f MB)\n",
            file.info("wnn_processed.rds")$size / 1e6))
```

The saved `.rds` carries both assays, all three reductions
(`pca`/`apca`-or-`lsi`/`wnn.umap`/`rna.umap`/`adt.umap`), the WNN neighbors
and graphs (`weighted.nn`/`wknn`/`wsnn`), the `RNA.weight` metadata, and
the WNN cluster identities. Downstream marker calling can run on the RNA
assay (`FindAllMarkers(obj, assay="RNA")`) or on the ADT assay
(`FindAllMarkers(obj, assay="ADT")`) — the latter usually gives cleaner
markers for the lineage-defined clusters.

---

## Batch variant — use INSTEAD of Steps 1–6 when invoked with args="batch"

Branch on `$ARGUMENTS == "batch"` at the top of the body. In batch mode:

- Skip all per-step figures.
- Skip the per-step "Assess and report" footers.
- Still save the canonical `wnn_processed.rds`.
- Print ONE final summary line.

```r
suppressPackageStartupMessages({ library(Seurat) })

obj <- readRDS("/path/to/cite_processed.rds")    # or multiome_processed.rds

# Ensure the second-modality reduction exists; if not, compute it (CITE-seq case shown).
if (!"apca" %in% Reductions(obj) && "ADT" %in% Assays(obj)) {
  DefaultAssay(obj) <- "ADT"
  VariableFeatures(obj) <- rownames(obj[["ADT"]])
  obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, assay = "ADT")
  obj <- ScaleData(obj, assay = "ADT")
  adt_npcs <- min(50, nrow(obj[["ADT"]]) - 1)
  obj <- RunPCA(obj, assay = "ADT", reduction.name = "apca",
                npcs = adt_npcs, approx = FALSE, verbose = FALSE)
  DefaultAssay(obj) <- "RNA"
}

second_red  <- if ("apca" %in% Reductions(obj)) "apca" else "lsi"
second_dims <- if (second_red == "apca") 1:18 else 2:30

second_weight <- if (second_red == "apca") "ADT.weight" else "ATAC.weight"
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list       = list("pca", second_red),
  dims.list            = list(1:30, second_dims),
  modality.weight.name = c("RNA.weight", second_weight),
  verbose              = FALSE
)
obj <- RunUMAP(obj, nn.name = "weighted.nn",
               reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
               verbose = FALSE)
obj <- FindClusters(obj, graph.name = "wsnn",
                    algorithm = 3, resolution = 2, verbose = FALSE)

saveRDS(obj, file = "wnn_processed.rds")
cat(sprintf("batch ok | %d cells | %d wnn clusters | RNA.weight median=%.2f | %.1f MB\n",
            ncol(obj),
            length(levels(Idents(obj))),
            median(obj$RNA.weight),
            file.info("wnn_processed.rds")$size / 1e6))
```

---

## Final response checklist

Summarize:
- which two reductions were combined (pca + apca / pca + lsi) and their dim ranges
- WNN cluster count vs RNA-only cluster count (the win-over-RNA-only signal)
- RNA.weight distribution (min/median/max) and which clusters sit predominantly above vs below 0.5 — i.e. RNA-driven vs ADT/ATAC-driven populations
- clustering settings (graph.name="wsnn", algorithm=3, resolution)
- whether the joint UMAP separates populations the per-modality UMAPs do not (the headline finding)
- figures shown (umap_wnn_clusters, umap_rna_only, umap_adt_only, modality_weights)
- saved files (wnn_processed.rds)
- caveats: WNN is sensitive to per-modality dims.list — wrong dim counts skew the weighting; small ADT panels (<10 markers) may not benefit; the joint UMAP is a visualization, not a definition — read clusters off `seurat_clusters`, not off polygons on the UMAP

---

## See also

- `seurat-cite-seq` — single-modality protein workup; use when the user
  wants ADT as overlay on the RNA UMAP, not joint clustering. WNN is the
  step-up when that recipe's `adt_qc_violins.png` shows lineage markers
  resolving sub-populations the RNA clusters miss.
- `seurat-rna-atac-integration` — RNA + ATAC WNN specifically. Same
  `FindMultiModalNeighbors` API; the difference is on the ATAC side
  (`Signac::RunTFIDF` → `RunSVD` → `lsi` reduction) and the
  `dims.list = list(1:30, 2:30)` convention (LSI component 1 dropped).
- `seurat-multimodal-reference-mapping` — when a pre-built WNN reference
  exists for your tissue (e.g. PBMC 162k), map a new query onto it
  instead of computing WNN from scratch.
- `seurat-scrna-v2` — the RNA-only workflow that produces the `pca`
  reduction WNN consumes.
