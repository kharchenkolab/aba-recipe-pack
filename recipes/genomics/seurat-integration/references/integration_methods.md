# Integration methods — picker + per-method semantics

When to load this: choosing between `CCAIntegration`,
`RPCAIntegration`, `HarmonyIntegration`, `JointPCAIntegration`,
`scVIIntegration`, `FastMNNIntegration`; deciding whether to switch
methods after a first pass; needing the v4 anchor-based fallback for
reproducing a pre-v5 analysis.

The SKILL.md body shows the canonical CCAIntegration call and three
inline variants (RPCA, Harmony, plus the SCTransform-path tweak). This
reference adds the full method set, per-method cost vs behavior,
recommended defaults by dataset shape, and the v4 `IntegrateData`
fallback.

## What `IntegrateLayers` actually dispatches on

Per `?IntegrateLayers` (Seurat 5.5.0), the signature is:

```
IntegrateLayers(object, method, orig.reduction = "pca", assay = NULL,
                features = NULL, layers = NULL,
                scale.layer = "scale.data", ...)
```

`method` is a **function reference**, not a string. The function-arg
pass-through (`...`) routes every other argument into the chosen
method's own function. So `k.anchor`, `k.filter`, `normalization.method`,
`dims`, etc. all live on the inner method's signature — read
`formals(CCAIntegration)` (or whichever method you chose) to see what
you can tune.

Built-in methods packaged inside Seurat itself
(`?IntegrateLayers` → "Integration Method Functions"):

- `CCAIntegration`
- `RPCAIntegration`
- `HarmonyIntegration`
- `JointPCAIntegration`

`scVIIntegration` and `FastMNNIntegration` are documented in the v5
introductory vignette but ship in companion packages, NOT in Seurat
itself. On a fresh Seurat install `exists("scVIIntegration")` returns
`FALSE` — installing `SeuratWrappers` (devtools from
`satijalab/seurat-wrappers`) is what adds them. Always gate with
`stopifnot(exists("scVIIntegration"))` before calling.

## Per-method picker

| Method | Speed | Memory | Best when | Caveats |
|---|---|---|---|---|
| `CCAIntegration` | slow | high | Most cell types shared across samples; vignette default for ifnb / pbmcsca | Can over-correct rare or sample-specific populations |
| `RPCAIntegration` | fast | low | Large cohorts (≥10 samples) or partly-shared populations | More conservative — may leave some residual batch effect |
| `HarmonyIntegration` | fastest | low | Quick first pass; iterative correction on existing PCA | Linear method; can under-correct strong nonlinear batch effects. Default `new.reduction = "harmony"`. |
| `JointPCAIntegration` | fast | medium | Joint PCA across samples (alternative to per-sample PCA + correction) | Newer; fewer published comparisons than CCA/RPCA |
| `FastMNNIntegration` | fast | medium | MNN-style anchor matching across samples | Needs `SeuratWrappers` + Bioc `batchelor` |
| `scVIIntegration` | slow (GPU helpful) | high | Atlas scale; deep-generative integration | Needs `SeuratWrappers` + `reticulate` + a Python env with `scvi-tools`. For atlas work prefer the standalone `scvi-integration` recipe. |

### How to pick at first run

1. **Default = CCAIntegration.** Matches the vignette; widely benchmarked.
2. **Switch to RPCA if** cohort is ≥10 samples / ≥100k cells, or the
   CCA pass collapses a population the user knows is distinct.
3. **Switch to Harmony if** you need a fast iterative pass (e.g.
   exploratory; harmonize across many timepoints quickly).
4. **Switch to scVI/FastMNN only if** the user explicitly asks for them,
   or the project already standardized on one.

Compare side-by-side by re-running `IntegrateLayers` with a different
`method` and a different `new.reduction` name — the corrections sit
together on the object and you can `DimPlot` each.

## Per-method args worth knowing

### `CCAIntegration` / `RPCAIntegration` (`?CCAIntegration`)

Both share the same inner signature. Common tunables (Seurat 5.5.0):

| Arg | Default | Effect |
|---|---|---|
| `new.reduction` | `"integrated.dr"` | Override per call (`"integrated.cca"`, `"integrated.rpca"`) so reductions stack side-by-side. |
| `normalization.method` | `c("LogNormalize","SCT")` | Set to `"SCT"` (case-sensitive) when Step 2 used `SCTransform`. |
| `dims` | `1:30` | Inner dims for anchor finding — usually leave at default; the recipe's `dims = 1:DIMS_CHOSEN` is the OUTER dims for neighbors/UMAP. |
| `k.anchor` | (passes to `FindIntegrationAnchors`, default 5) | Raise to 15–20 to strengthen integration when CCA leaves visible batch effect. |
| `k.weight` | 100 | Number of neighbors for anchor weighting. |
| `k.filter` | NA | Anchor filtering; NA disables. |
| `sample.tree` | NULL | Matrix encoding the order pairwise integrations are performed. NULL = auto. |

Example — increasing `k.anchor`:

```r
obj <- IntegrateLayers(
  object         = obj,
  method         = CCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.cca",
  k.anchor       = 20,                 # passes through ... to FindIntegrationAnchors
  verbose        = FALSE
)
```

### `HarmonyIntegration` (`?HarmonyIntegration`)

Different inner signature (Harmony is a re-implementation, not an
anchor flow). Common tunables:

| Arg | Default | Effect |
|---|---|---|
| `new.reduction` | `"harmony"` | Default name; downstream `FindNeighbors/RunUMAP` point at this. |
| `theta` | NULL (auto) | Diversity penalty; higher = stronger correction. |
| `lambda` | NULL (auto) | Ridge regularization. |
| `sigma` | 0.1 | Kernel bandwidth. |
| `max.iter.harmony` | 10 | Outer iterations. |
| `key` | `"harmony_"` | Reduction key prefix. |

```r
obj <- IntegrateLayers(
  object         = obj,
  method         = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction  = "harmony",
  theta          = 2,                  # stronger correction across samples
  verbose        = FALSE
)
```

### `JointPCAIntegration`

Same outer call shape; computes a joint PCA across samples instead of
correcting an existing per-sample PCA.

```r
obj <- IntegrateLayers(
  object         = obj,
  method         = JointPCAIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.jpca",
  verbose        = FALSE
)
```

### `scVIIntegration` (companion: `SeuratWrappers`)

```r
# Verify the function exists; SeuratWrappers must be installed AND attached.
# scvi-tools must be reachable via reticulate.
stopifnot(exists("scVIIntegration"))   # FALSE on fresh Seurat install

obj <- IntegrateLayers(
  object         = obj,
  method         = scVIIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.scvi",
  verbose        = FALSE
)
```

For atlas-scale work (>500k cells, label transfer via scANVI, GPU
pinning), prefer the standalone `scvi-integration` recipe instead of
this wrapper — it handles GPU env, batch keys, and label-transfer
plumbing the wrapper does not.

### `FastMNNIntegration` (companion: `SeuratWrappers` + Bioc `batchelor`)

```r
stopifnot(exists("FastMNNIntegration"))

obj <- IntegrateLayers(
  object         = obj,
  method         = FastMNNIntegration,
  orig.reduction = "pca",
  new.reduction  = "integrated.mnn",
  verbose        = FALSE
)
```

## SCTransform normalization wiring

If Step 2 took the SCTransform variant, `IntegrateLayers` must know.
Two routes — pick whichever matches the call style:

```r
# Route 1: pass normalization.method through ...
obj <- IntegrateLayers(
  object               = obj,
  method               = CCAIntegration,
  orig.reduction       = "pca",
  normalization.method = "SCT",         # CASE-SENSITIVE — "SCT", not "sct"
  new.reduction        = "integrated.cca",
  verbose              = FALSE
)
```

```r
# Route 2: pass the assay name explicitly (the SCT assay is named "SCT")
obj <- IntegrateLayers(
  object         = obj,
  method         = CCAIntegration,
  orig.reduction = "pca",
  assay          = "SCT",               # see ?IntegrateLayers — assay arg
  new.reduction  = "integrated.cca",
  verbose        = FALSE
)
```

After SCT-path integration, DE downstream needs
`PrepSCTFindMarkers(obj)` before `FindAllMarkers` / `FindMarkers`.

## Classic v4 anchor pipeline — fallback for reproducing a v4 analysis

The pre-v5 route uses a LIST of separate objects, finds anchors across
them, and materializes a new `integrated` assay. **Prefer the v5
`IntegrateLayers` flow in SKILL.md; use this only to reproduce a v4
analysis exactly.**

```r
# Starts from a LIST of per-sample objects (NOT a merged-and-split object).
obj.list <- lapply(objs, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000,
                       verbose = FALSE)
})
features <- SelectIntegrationFeatures(obj.list)
anchors  <- FindIntegrationAnchors(
  object.list     = obj.list,
  anchor.features = features,
  reduction       = "cca"       # CASE-SENSITIVE: "cca" / "rpca" / "jpca" / "rlsi"
                                #   — NEVER "pca"; not a valid anchor reduction.
)
combined <- IntegrateData(anchorset = anchors)
DefaultAssay(combined) <- "integrated"
combined <- ScaleData(combined, verbose = FALSE)   # MUST ScaleData BEFORE RunPCA
combined <- RunPCA(combined, verbose = FALSE)      #   on the integrated assay.
combined <- FindNeighbors(combined, dims = 1:30, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)
combined <- RunUMAP(combined, dims = 1:30, verbose = FALSE)
# DE in the v4 flow runs on the `RNA` assay, not the `integrated` assay:
DefaultAssay(combined) <- "RNA"
```

Key v4-vs-v5 contrasts:

| v4 | v5 |
|---|---|
| List of N objects in, integrated assay out | One merged object, RNA assay split into N layers, corrected reduction out |
| `FindIntegrationAnchors(reduction = "cca" \| "rpca" \| ...)` (string, case-sensitive) | `IntegrateLayers(method = CCAIntegration)` (function reference) |
| Creates a new `integrated` assay | Adds a corrected reduction next to `pca` |
| DE on the `RNA` assay (set `DefaultAssay`) | DE on the joined RNA layers (call `JoinLayers` first) |
| `PrepSCTIntegration(obj.list)` then `FindIntegrationAnchors(... normalization.method = "SCT")` | `SCTransform(obj)` then `IntegrateLayers(... normalization.method = "SCT")` |

## Comparing two methods on the same object

```r
# CCA pass
obj <- IntegrateLayers(obj, method = CCAIntegration,
                       orig.reduction = "pca",
                       new.reduction  = "integrated.cca", verbose = FALSE)
# RPCA pass on the same object — reductions stack
obj <- IntegrateLayers(obj, method = RPCAIntegration,
                       orig.reduction = "pca",
                       new.reduction  = "integrated.rpca", verbose = FALSE)

# UMAP each into its own reduction
obj <- RunUMAP(obj, reduction = "integrated.cca",  dims = 1:30,
               reduction.name = "umap.cca",  verbose = FALSE)
obj <- RunUMAP(obj, reduction = "integrated.rpca", dims = 1:30,
               reduction.name = "umap.rpca", verbose = FALSE)

# Compare visually — same axis scale, side by side
p_cca  <- DimPlot(obj, reduction = "umap.cca",  group.by = "sample")
p_rpca <- DimPlot(obj, reduction = "umap.rpca", group.by = "sample")
# combine via patchwork: p_cca | p_rpca
```
