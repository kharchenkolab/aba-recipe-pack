---
name: conos-integration
description: Joint/integrative analysis of multiple scRNA-seq samples with the conos R6 pipeline (joint graph, clustering, label propagation, DE)
when_to_use: Two or more scRNA-seq samples to integrate into one joint graph; want cross-sample clusters, a shared embedding, label transfer, or between-group differential expression
requires_tools: [run_r]
capabilities_needed: [conos, pagoda2, Matrix, igraph]
keywords: [conos, integration, batch correction, multi-sample, joint graph, mNN, label propagation, leiden, scRNA-seq, kharchenkolab, R6]
produces: [Conos object, joint clusters, joint 2D embedding, propagated labels, per-cluster markers, between-group DE tables, joint.lstar.zarr]
domain: genomics
source: github:kharchenkolab/conos + vignette https://github.com/kharchenkolab/conos/blob/main/doc/walkthrough.md
---
# conos multi-sample integration

conos is an R6 package for joint analysis of a panel of single-cell
samples. Each sample must first be preprocessed (a `Pagoda2` R6 object, a
Seurat object, or mixed). The `Conos` object holds the list of samples
and builds a joint kNN graph across them, on which you cluster, embed,
propagate labels, and run DE. Methods mutate the object in place via `$`.

**Version 2.0 (June 2026): API verbs renamed** to the `runX(method=)`
vocabulary shared with pagoda2 2.0. The old names still work but warn:
`buildGraph→runGraph`, `findCommunities→runClustering`,
`embedGraph→runEmbedding`, `getDifferentialGenes→runMarkers`. Methods
also accept method names as strings (`"leiden"`, `"walktrap"`, etc.)
instead of function refs. ComplexHeatmap is gone from the dependency
chain — heatmaps render via `sccore`'s native grid engine.

## Approach
1. `library(conos)`. Start from a **named list** of count matrices (genes x
   cells), one per sample. Ensure cell names are globally unique across
   samples (prefix by sample id) — `any(duplicated(unlist(lapply(panel, colnames))))`
   must be FALSE.
2. Preprocess each sample into a `Pagoda2` object with a PCA reduction. **pagoda2 has
   two API generations and this recipe supports BOTH** — 2.0 (`devel`) and 1.x
   (CRAN/`main`) — which construct differently; a snippet for one ERRORS on the other
   (`Pagoda2$from` / `$run()` don't exist on 1.x → "attempt to apply non-function").
   Detect the flavor once and route. conos 2.0 reads counts from either flavor
   internally (its accessor shim falls back `getExpressionBlock()` → `$counts`), so
   ONLY per-sample construction branches — requires **conos >= 2.0**.
   ```r
   pagoda2_is_devel <- function() is.function(tryCatch(pagoda2::Pagoda2$from, error = function(e) NULL))
   preprocess_p2 <- function(cm, n.cores = 1) {
     if (pagoda2_is_devel()) {        # pagoda2 >= 2.0 (devel): unified constructor + step pipeline
       Pagoda2$from(cm, n.cores = n.cores, verbose = FALSE)$run(steps = c("variance", "pca"), verbose = FALSE)
     } else {                          # pagoda2 1.x (CRAN/main): classic constructor + explicit steps
       p <- Pagoda2$new(cm, n.cores = n.cores, log.scale = TRUE, verbose = FALSE)
       p$adjustVariance(plot = FALSE, verbose = FALSE)
       p$calculatePcaReduction(nPcs = 30, n.odgenes = 2000, verbose = FALSE)
       p
     }
   }
   panel.preprocessed <- lapply(panel, preprocess_p2)   # variance + PCA = the minimum conos needs
   ```
   File readers differ by flavor: on **devel**, `Pagoda2$from/from10x/from10xH5/fromAnnData/`
   `fromH5Seurat/fromLoom/fromLstar` read any format (no Python); on **CRAN 1.x** only 10x is
   native (`read10xMatrix(path)` → `preprocess_p2()`), other formats need devel pagoda2 (or
   convert upstream then `preprocess_p2(cm)`). For Seurat, pass Seurat objects directly (or a
   Pagoda2 + Seurat mix); conos aligns on the PCA.
3. `con <- Conos$new(panel.preprocessed, n.cores = 1)`.
4. (Optional, recommended for multimodal panels) `con$planIntegration()` — polls
   each sample's modalities (pagoda2 facets / Seurat assays), reports per-modality
   shared-feature count and per-pair overlap with a usable / marginal / not-integrable
   verdict, and records the resolved default modality on the object. Surfaces e.g.
   independently-called scATAC peaks that won't reconcile without consistent
   pre-processing.
5. Build the joint graph (new verb `runGraph`; `buildGraph` is a deprecated alias):
   ```r
   con$runGraph(k = 30, k.self = 5, space = 'PCA', ncomps = 30,
                 n.odgenes = 2000, matching.method = 'mNN', metric = 'angular',
                 score.component.variance = TRUE)
   ```
   `space=` may be 'PCA' (fast default), 'CPCA' (more distortion-robust),
   'CCA' (low-similarity / cross-species), or 'genes' (same platform).
   To recompute a space, clear its cache: `con$pairs$PCA <- NULL`. For very
   large panels (memory) pass `pairs.storage = "drop"` (recompute on next
   build) or `"disk"` (offload to disk, restored on next build) — default is
   `"keep"` (in-memory).
6. Joint clustering (new verb `runClustering`; `findCommunities` is deprecated):
   ```r
   con$runClustering(method = 'leiden', resolution = 1)
   ```
   `method` accepts strings: `'leiden'`, `'walktrap'` (hierarchical, also takes
   `steps = 8-10`), `'louvain'`/`'multilevel'`, `'infomap'`, `'fastgreedy'`,
   `'labelprop'`, `'leadingeigen'`. (You can still pass a function, e.g.
   `leiden.community` or `igraph::walktrap.community`.) Results stored as a list
   under `con$clusters$<name>` (e.g. `con$clusters$leiden$groups`). To sweep
   resolutions: `con$scanResolution(method = 'leiden', resolutions = seq(0.3, 2, 0.1))`
   returns cluster count + modularity per value.
7. Embed the joint graph (new verb `runEmbedding`; `embedGraph` is deprecated):
   ```r
   con$runEmbedding(method = 'largeVis')           # default
   # or:
   con$runEmbedding(method = 'UMAP', min.dist = 0.01, spread = 15)
   ```
   You MUST call this before plotting. Name multiple embeddings with
   `embedding.name=`.
8. Visualize: `con$plotPanel(clustering = 'leiden', use.common.embedding = TRUE, font.size = 4)`
   (per-sample small multiples on the JOINT embedding). **`use.common.embedding = TRUE` is
   required** — the preprocessing above builds only variance + PCA, so samples have no
   per-sample embedding and the default `FALSE` errors ("No 'tSNE' embedding presented in the
   samples"). (`use.local.clusters = TRUE` shows each sample's own clusters.) Also
   `con$plotGraph(color.by = 'sample', alpha = 0.1)` /
   `con$plotGraph(gene = 'GZMK')` (the joint embedding). Both wrap
   `sccore::embeddingPlot`.
9. Label propagation: `info <- con$propagateLabels(labels = cellannot, verbose = TRUE)`
   transfers a named factor of annotations from labeled to unlabeled cells.
   Returns `$labels`, `$uncertainty`, `$label.distribution`.
10. Cluster markers (new verb `runMarkers`; `getDifferentialGenes` is deprecated):
    ```r
    de <- con$runMarkers(groups = new.annot, append.auc = TRUE)
    ```
    Per-group table with `M`, `Z`, `PValue`, `PAdj`, `AUC`, `Specificity`,
    `Precision`. Mixed Seurat + pagoda2 panels are now scored uniformly through
    `sccore::matrixDE` (numerically equivalent to the per-sample pagoda2 path).
    Heatmap (no ComplexHeatmap dep — uses sccore's grid engine):
    ```r
    plotDEheatmap(con, as.factor(groups), de, n.genes.per.cluster = 5,
                  column.metadata = list(samples = con$getDatasetPerCell()))
    ```
    `return.details = TRUE` now returns a heatmap *spec* (input to
    `sccore::drawHeatmap`), not a ComplexHeatmap object.
11. Specific-marker dot plot (new): `con$plotMarkerDotPlot(n.genes.per.group = 5, min.auc = 0.6)`.
    (It computes the markers internally from the joint clustering — there is NO `de` argument;
    the first positional arg is `clustering`, so do not pass a precomputed `de` table here.)
    Ranks by a "balanced" rule (precision × expression-fraction harmonic mean + `min.auc`
    discrimination floor) and assigns each gene to one best cluster, so housekeeping /
    mitochondrial genes don't dominate.
12. Between-group DE: define `samplegroups <- list(bm = c(...), cb = c(...))`, then
    `getPerCellTypeDE(con, groups = as.factor(new.annot), sample.groups = samplegroups,
    ref.level = 'bm')` (pseudobulk per cluster, DESeq2 under the hood). For custom
    models pull meta-cell counts with `con$getClusterCountMatrices()`. Whole joint
    matrix: `con$getJointCountMatrix()`.

## Key decisions / parameters
- `space`: PCA is the fast default; escalate to CPCA/CCA only when samples
  are highly dissimilar. With 'angular' metric keep ~30 components (don't cut
  to the variance elbow); with 'L2' fewer components can be better.
- `runClustering` `resolution` (leiden) trades cluster granularity; walktrap
  `steps` higher = finer + slower. Use `scanResolution()` to choose.
- Force tighter alignment with `runGraph(alignment.strength = 0.3, ...)` and
  `balance.edge.weights = ` to rebalance by a factor (e.g. tissue).
- `greedyModularityCut(con$clusters$walktrap$result, N)` cuts the walktrap
  dendrogram to N clusters for hierarchical exploration.
- Memory on large panels: `runGraph(pairs.storage = "drop")` or `"disk"`.

## Caveats
- Inputs are preprocessed per-sample objects, not raw matrices — run
  pagoda2 (R6) or Seurat first.
- Cell names must be unique across the whole panel or the graph is corrupt.
- Plotting needs an explicit prior `runEmbedding` call.
- Joint clusters are comparable across samples; per-sample local clusters
  (`use.local.clusters = TRUE`) are not.
- Removed in 2.0: `p2app4conos()` (web app — pagoda2 1.x's app is gone) and
  `saveConosForScanPy()` (the bespoke "corrected pseudo-expression" export).
  For interchange use the **lstar** / Zarr path or `con$getJointCountMatrix()`.
- `buildGraph / findCommunities / embedGraph / getDifferentialGenes` still
  work in 2.0 but emit deprecation warnings — prefer the `runX` verbs.

## In ABA
- Install via `ensure_capability("conos")` — both conos and pagoda2 2.0 are
  pure-CRAN-ish builds now: PPM binaries cover Matrix/igraph; conos no
  longer Suggests ComplexHeatmap (the heatmap renderer moved to `sccore`'s
  native grid backend), so there's no Bioconductor dep on the critical path
  anymore. `igraph` still needs the GLPK system lib (conda `glpk`). Run all
  steps in `run_r`.
- Pagoda2's R6 API supports a multimodal "facet" design — each modality
  (RNA, ADT/CITE-seq, ATAC) is its own facet with a view model
  (plain / CLR / TF-IDF). For single-modality scRNA-seq the recipe above is
  what you want; cross-modality WNN-style work is a separate recipe.
- For zarr/HDF5 interchange (round-tripping with AnnData / Seurat, or feeding
  ABA's viewer), the `lstar` package is the canonical path. Build an lstar
  Dataset from the **live** object, then write a store:
  ```r
  d <- lstar::write_conos(con)                 # joint Conos object   -> lstar Dataset
  d <- lstar::read_pagoda2(p2)                 # a single Pagoda2 obj -> lstar Dataset
  lstar::lstar_write_viewer(d, "joint.lstar.zarr")  # -> optimized .lstar.zarr store
                                                    #    (precomputes DE/HVG; no banner)
  ```
  (lstar's `.rds`/CLI converter reads only Seurat/SCE — pagoda2/conos ingest
  goes through these R functions on the live object, not a saved `.rds`.)
### Offer an interactive view

**Required final step — not optional.** After writing `joint.lstar.zarr`,
**proactively offer to open it in ABA's interactive viewer** — call
`open_viewer(file_path="joint.lstar.zarr")` and present the returned link in your
closing message. Point the viewer at the store (or an `.h5ad`), never the raw
conos `.rds`. If `open_viewer` returns `ok:false`, relay the error.
