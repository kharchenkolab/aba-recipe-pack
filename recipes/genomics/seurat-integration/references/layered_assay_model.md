# Seurat v5 split-layer assay model

When to load this: you need to understand how the v5 `Assay5` class
stores split counts/data/scale.data layers; you hit `The following
layers are already split: …` from `split()`; you're wondering when to
call `JoinLayers`; you need to inspect layer state on a merged object.

The SKILL.md body works fine without reading this — the canonical
`JoinLayers(obj) → split(obj[["RNA"]], f = sample)` is enough for the
standard flow. Load this when something looks wrong with the layer
shape, or when adapting the flow to a non-standard input.

## The two assay classes

A Seurat v5 object's RNA assay can be one of:

- **`Assay`** (v4 legacy class) — one counts matrix, one data matrix,
  one scale.data matrix. No per-sample split possible. Older saved
  `.rds` objects often arrive in this shape.
- **`Assay5`** (v5 class) — supports multiple **named layers** per
  matrix slot. Layers are named `counts.<x>`, `data.<x>`,
  `scale.data.<x>` after a split; they're named `counts`, `data`,
  `scale.data` (no suffix) when consolidated.

Check which one you have:

```r
class(obj[["RNA"]])              # "Assay5" if v5, "Assay" if v4 legacy
Layers(obj[["RNA"]])             # list of layer names in current state
```

If `class(obj[["RNA"]])` is `"Assay"` (v4), convert before integrating:

```r
obj[["RNA"]] <- as(obj[["RNA"]], "Assay5")
```

## Layer naming after `merge()` of N objects

When you `merge()` v5 objects, the result's RNA assay automatically
ends up split by `orig.ident` — Seurat preserves each input's data as
its own layer:

```r
obj <- merge(objs[[1]], y = objs[-1], add.cell.ids = names(objs))
Layers(obj[["RNA"]])
# -> "counts.ctrl"  "counts.stim"
#    (and "data.ctrl" "data.stim" if any input was already normalized)
```

This is the **already-split state**. `split(obj[["RNA"]], f = ...)`
**refuses to operate on already-split layers** — it errors with:

```
Error in .local(x, f, ...) :
  The following layers are already split: 'counts.ctrl', 'counts.stim'
```

The fix is `JoinLayers(obj)` BEFORE `split()`. It consolidates back
to single `counts`/`data`/`scale.data` layers; then `split(... , f = ...)`
can re-split on YOUR chosen factor (e.g. `obj$sample` if you've renamed
`orig.ident`).

## `JoinLayers` is the safe-by-default operation

`JoinLayers` is a no-op when nothing is split, and consolidates when
something is. So calling it before `split()` always works, regardless
of the input shape:

```r
# Safe — works whether obj is fresh, post-merge, post-integration, or already joined.
obj <- JoinLayers(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
```

This pattern is what the SKILL.md Step 1 uses, and it's why.

## When `JoinLayers` is REQUIRED

| Operation | Split layers OK? | Why |
|---|---|---|
| `NormalizeData`, `FindVariableFeatures`, `ScaleData`, `RunPCA` | YES (per-layer auto) | Each operates per layer; that's the whole point of the v5 layout. |
| `IntegrateLayers` | YES (REQUIRED split) | Needs the split layers as the per-sample views to align. |
| `FindNeighbors`, `FindClusters`, `RunUMAP` (on a reduction) | YES — works either way | These operate on the reduction, not on the assay layers. |
| `FindAllMarkers`, `FindMarkers` (Wilcoxon, default) | **NO — JoinLayers first** | DE iterates over the data layer; a split layout means it can't find a single `data` slot, so it errors with "Layer 'data' is not present" or similar. |
| `AverageExpression`, `AggregateExpression` | **NO — JoinLayers first** | Same — operates over a single `counts` / `data` layer. |
| `saveRDS` | YES (either) | Layers are preserved as-is on disk; consumer can re-join on load. |

The pattern in the SKILL.md: split before Step 4 (IntegrateLayers
needs it), JoinLayers after Step 6 (Step 7 / downstream DE needs it).

## Operating on a specific layer

Most Seurat functions take a `layer` argument. After splitting, you
can address a single sample's layer explicitly:

```r
# Get the raw counts for one sample
ctrl_counts <- GetAssayData(obj, assay = "RNA", layer = "counts.ctrl")

# Set a layer (rare; usually you set via NormalizeData etc.)
LayerData(obj, layer = "data.ctrl") <- some_matrix
```

To list all layers and their dimensions:

```r
sapply(Layers(obj[["RNA"]]), function(L) {
  m <- GetAssayData(obj, layer = L)
  paste0(nrow(m), " x ", ncol(m))
})
```

## SCTransform interaction with layers

`SCTransform` on a split-layer object operates per layer and produces
a NEW assay called `SCT`. The RNA assay's layers are unaffected. After
SCT, the SCT assay has its own per-sample variance-stabilized counts
+ Pearson residuals; `IntegrateLayers(method = ..., assay = "SCT")` (or
`normalization.method = "SCT"`) integrates THOSE.

```r
obj <- SCTransform(obj, verbose = FALSE)
Assays(obj)                  # "RNA" "SCT"
DefaultAssay(obj)            # set to "SCT" by SCTransform
Layers(obj[["SCT"]])         # one per sample
```

After SCT-path integration, before DE call `PrepSCTFindMarkers(obj)`
to fix the per-cell variance estimates across the merged samples;
`FindAllMarkers` then runs on the joined SCT assay.

## When you legitimately want to KEEP layers split

- Per-sample marker discovery — explicit per-sample Wilcoxon on one
  layer at a time. Rare. Pseudobulk DE on per-sample aggregates is
  the standard practice instead.
- Pseudobulk DE — you ALREADY operate on per-sample aggregates, so
  the layer split is incidental; aggregate first
  (`AggregateExpression(obj, group.by = c("sample", "cluster"), …)`,
  which itself joins layers internally).

In practice: split for the upstream integration, join before any DE,
don't re-split unless you have a specific per-sample analysis in
mind.

## Inspecting a freshly-loaded `.rds`

If you receive a saved Seurat object from a colleague:

```r
obj <- readRDS("/path/to/obj.rds")
class(obj[["RNA"]])              # Assay vs Assay5
Layers(obj[["RNA"]])             # split state
Reductions(obj)                  # what reductions exist
table(obj$orig.ident)            # how many samples — does it have a sample column?
```

A common arrival shape: v5 Assay5 with `counts.<sample>` +
`data.<sample>` layers (already split), with `pca` already computed. To
integrate from there: `JoinLayers(obj)` then re-`split` by your chosen
covariate, redo PCA on the joined+re-split, then `IntegrateLayers`.
You don't have to recompute normalization — only PCA, since splitting
changes the per-layer scale.data inputs.
