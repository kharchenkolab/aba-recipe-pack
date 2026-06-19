# Multimodal mapping internals

What `FindTransferAnchors` + `MapQuery` actually do under the hood, the
supervised-PCA (`spca`) construction, SCT vs LogNormalize anchor finding,
and dims-range guidance. Load this when the agent needs to defend a
parameter choice or unpack a step that looks like a black box.

Sources consulted: Hao et al. *Cell* 184, 3573–3587 (2021) — "Integrated
analysis of multimodal single-cell data" (the WNN + multimodal mapping
paper); Seurat source `R/integration.R` (functions
`FindTransferAnchors`, `MapQuery`, `TransferData`, `IntegrateEmbeddings`,
`ProjectUMAP`), tested against Seurat 5.5.0; the Satija Lab vignette at
`satijalab.org/seurat/articles/multimodal_reference_mapping`.

## What MapQuery composes

`MapQuery()` is a thin wrapper that calls three Seurat functions in
sequence. Knowing this lets you (a) understand which outputs come from
which step, and (b) call the sub-functions individually when you need
finer control (e.g. you want anchors but not the UMAP projection).

```r
# Conceptually (simplified):
query <- TransferData(anchorset = anchors, refdata = refdata,
                       reference = reference, query = query, ...)
query <- IntegrateEmbeddings(anchorset = anchors, reference = reference,
                              query = query,
                              new.reduction.name = "ref.spca", ...)
query <- ProjectUMAP(query = query, query.reduction = "ref.spca",
                     reference = reference,
                     reference.reduction = "spca",
                     reduction.model = "wnn.umap",
                     reduction.name = "ref.umap", ...)
```

What each does:

- **`TransferData`** — for each entry in `refdata`, transfers labels (if
  the value names a metadata column) or imputes assays (if the value
  names a reference assay). Categorical transfers produce a
  `predicted.<key>` column plus a per-cell `predicted.<key>.score`
  (maximum prediction probability across labels). Continuous transfers
  produce a new assay on the query named by the LIST KEY (so
  `predicted_ADT = "ADT"` writes `query[["predicted_ADT"]]`).
- **`IntegrateEmbeddings`** — projects the query cells into the
  reference's `spca` space (creates `query[["ref.spca"]]`). Uses the
  anchor pairs to learn a transform from query expression to reference
  spca coordinates.
- **`ProjectUMAP`** — runs `uwot::umap_transform` to place the projected
  query cells onto the reference's `wnn.umap` (creates
  `query[["ref.umap"]]`). The reference's `wnn.umap` must have been
  saved with `return.model = TRUE` for this to work.

If `query[["ref.umap"]]` is empty after `MapQuery`, the reference's
`wnn.umap` was likely saved without the UMAP model — re-run
`RunUMAP(..., return.model = TRUE)` on the reference and re-save.

## What supervised PCA (`spca`) is and why it matters

`RunSPCA()` (Signac/Seurat) computes a PCA on an assay *constrained* by a
neighbor graph. The "supervision" comes from the graph: typically the
`wsnn` graph from `FindMultiModalNeighbors`, which encodes joint cell-cell
similarities across modalities. Components are chosen so they preserve
relationships in that graph as well as possible.

Why this beats vanilla PCA for reference mapping:

- Vanilla PCA on RNA preserves variance but mixes modality-specific
  signals (e.g. RNA components dominated by RNA-only contamination —
  cycling, MT-content — push apart cells that the joint analysis would
  group together).
- `spca` was trained against the reference's WNN graph, so its top
  components already separate the cell types the reference's labels
  encode. Anchors found in `spca` space are therefore "cell-type
  aware" — much more useful than anchors found in modality-imbalanced
  RNA PCA.

You don't compute `spca` in this recipe — it must exist on the
reference before mapping. Building one requires the WNN graph; see the
`seurat-wnn-multimodal` recipe + `RunSPCA(assay = "RNA", graph = "wsnn")`.

```r
# Confirm spca is present and has the expected width:
print(dim(Embeddings(reference, "spca")))         # cells x ncomp
print(head(Loadings(reference, "spca"))[, 1:3])   # gene x ncomp
```

## SCT vs LogNormalize anchor finding

`FindTransferAnchors(normalization.method = ...)` is the single switch
that selects the matching mode. Both reference and query must be in the
same normalization regime — this is the highest-frequency source of
silent wrongness in this workflow.

**SCT mode** (`normalization.method = "SCT"`):
- Both reference and query must carry an `SCT` assay.
- `recompute.residuals = TRUE` (default) recomputes the query's SCT
  residuals using the reference's SCT model when SCTransform versions
  differ. Seurat skips the recompute if the versions match — so it's
  harmless to leave on.
- Anchors are found in the reference's `spca` if you pass
  `reference.reduction = "spca"`, otherwise in the reference's standard
  `pca`.

**LogNormalize mode** (`normalization.method = "LogNormalize"`):
- Both reference and query must have been `NormalizeData`'d (LogNorm on
  the `RNA` assay).
- `recompute.residuals` is harmless but **irrelevant** — there are no
  SCT residuals to recompute.
- Anchors are still found in whatever reference reduction you pass via
  `reference.reduction`. If the LogNormalize-d reference has an `spca`
  (computed against its WNN graph), pass it; otherwise pass `"pca"`.

How to tell which mode the reference is in:

```r
print(Assays(reference))
# Has "SCT" -> SCTransform-normalized -> use normalization.method = "SCT"
# Only "RNA" (no "SCT") -> LogNormalize-d -> use normalization.method = "LogNormalize"
```

A reference can carry BOTH `RNA` and `SCT` assays; in that case SCT is
the active normalization for mapping (the recipe targets SCT).

## `dims` range — what 1:50 means

`dims = 1:50` selects the first 50 components of the reference's chosen
reduction (here `spca`). The PBMC 162k reference was computed with 50
supervised components; using fewer truncates resolution, using more is
impossible (the components don't exist).

For OTHER references:

```r
# Inspect the available width of the supervised PCA:
print(dim(Embeddings(reference, "spca")))   # rows = cells, cols = max dims
# Use 1:cols (or a truncation if you have a reason — e.g. the elbow plot
# of the reference's spca eigenvalues plateaus earlier).
```

Truncating below the reference's full width is rarely worth it; the
later components carry the long-tail discriminative signal that
separates fine cell types (e.g. CD8 naive vs CD8 effector). For
custom-built references where you computed `spca` with `npcs = 100`,
using all 100 is fine.

## Per-step alternatives — when you'd skip MapQuery

You'd call the sub-functions individually when:

- **You want to filter on anchor quality before label transfer.**
  ```r
  # Inspect anchor scores
  hist(anchors@anchors$score)
  # Filter to high-confidence anchors and rebuild the anchorset (advanced)
  ```
- **You only want labels, not the UMAP projection.** Drop the
  `ProjectUMAP` step (skip `MapQuery` and call `TransferData` directly).
- **You want labels at MULTIPLE granularities AND custom downstream
  embedding logic.** Run `TransferData` once with a long `refdata` list,
  then handle embedding yourself.

For most workflows `MapQuery` is the right call — it bakes in the
defaults from the Hao et al. vignette and saves a lot of plumbing.

## Anchor-count sanity heuristic

After `FindTransferAnchors`, check `nrow(anchors@anchors)`:

| Ratio (anchors / query cells) | Interpretation |
|---|---|
| > 100% (i.e. multiple anchors per query cell on average) | Strong match — every query cell has multiple reference neighbors |
| 50–100% | Healthy — typical for tissue-matched references |
| 10–50% | Marginal — partial match, some query populations may not be represented |
| < 10% (especially < 1%) | Poor — wrong reference, wrong species, or wrong gene-symbol convention |

The recipe's heuristic ("if anchors / query cells < 1%, stop and
reconsider") is calibrated to the lower bound. Above that, the
`predicted.celltype.l*.score` distribution per Step 5 is the next
sanity gate.
