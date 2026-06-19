---
name: seurat-reference-mapping
description: Map a scRNA-seq query onto a labeled Seurat reference (R/Seurat v5) — find transfer anchors between reference and query, transfer cell-type labels and metadata via TransferData, and project query cells into the reference's PCA + UMAP space via MapQuery. The reference's UMAP model must be saved with return.model = TRUE. Covers the generic case (user-supplied reference Seurat object) and points at Azimuth for pre-built tissue references.
when_to_use: You have a labeled reference Seurat object (annotated cell types, a saved PCA, and a UMAP model) plus an unlabeled query scRNA-seq dataset, and you want to (a) transfer the reference's cell-type labels onto the query and/or (b) project query cells into the reference's UMAP space for visualization on the SAME embedding. Use when the user names FindTransferAnchors / TransferData / MapQuery / reference mapping / label transfer / project onto atlas in an R/Seurat session. For pre-built human tissue references (PBMC, lung, kidney, heart, …) use Azimuth — see references/azimuth_alternative.md.
avoid_when: You want to integrate query and reference symmetrically (use seurat-integration — IntegrateLayers — to learn a joint embedding instead). The reference Seurat object lacks a saved UMAP model (`return.model = TRUE` was never set) — re-run `RunUMAP` on the reference with `return.model = TRUE` first, or use `TransferData` only (no UMAP projection).
invocation: interactive+batch
requires_tools: [run_r]
capabilities_needed: [Seurat]
keywords: [Seurat, Seurat v5, reference mapping, label transfer, FindTransferAnchors, TransferData, MapQuery, IntegrateEmbeddings, ProjectUMAP, query, reference, predicted.id, prediction.score, ref.umap, ref.pca, return.model, Azimuth, cell type annotation, atlas, R, scRNA-seq, single cell]
produces: [umap_reference_celltype.png, umap_query_predicted.png, umap_reference_query_side_by_side.png, query_predicted_labels.csv, query_mapped.rds]
domain: genomics
resource_profile: "small (~1–5 min for a 5–20k cell query against a 10–50k cell reference); MapQuery memory ~scales with cells in query × ref PCs)"
source: "Seurat v5 reference mapping vignette (Satija Lab) — satijalab.org/seurat/articles/integration_mapping; complemented by ?FindTransferAnchors, ?TransferData, ?MapQuery reference pages."
---

# Reference-based mapping of a scRNA-seq query onto a labeled Seurat reference (v5)

Reference mapping projects a new query dataset onto an existing **labeled
reference** so the query inherits the reference's coordinate system and
cell-type labels. Two outputs in one workflow:

- **Label transfer** — every query cell gets a predicted cell type (with a
  per-cell score) sourced from the reference's annotations.
- **UMAP projection** — query cells land on the reference's UMAP via the
  reference's PCA basis + the saved UMAP model. Reference and query become
  visually comparable on the SAME embedding.

This is **asymmetric**: only the query is altered. The reference's
expression, PCA, UMAP, and labels are untouched. If you want a symmetric
joint embedding that re-learns from both datasets, use `seurat-integration`
(`IntegrateLayers`) instead.

**Seurat version requirement: v5.0 or newer.** `MapQuery` exists in v4 too
but the v5 flow and `return.model` semantics are the canonical path.

## Bundled references — load on demand

This recipe is self-contained for the generic case (user-supplied
reference + query). For deeper detail, load the matching reference file
ONLY when needed:

- `references/mapping_internals.md` — what `FindTransferAnchors` /
  `TransferData` / `MapQuery` actually do; arg semantics, the
  `pcaproject` vs `cca` reduction options, score interpretation
  (`prediction.score.max`, `mapping.score`), and how to debug
  low-anchor or low-confidence outcomes.
- `references/azimuth_alternative.md` — pre-built human tissue
  references via Azimuth (`pbmcref`, `lungref`, …); when to use the
  Azimuth wrapper instead of the manual flow; how its outputs map to
  this recipe's outputs.
- `references/figure_style.md` — Seurat-collection figure conventions
  (palette, alpha-poke, `theme_cowplot()` placement vs `NoLegend()`,
  ggsave defaults, `coord_fixed()`, patchwork side-by-side).
  Duplicated across the seurat-* collection.

## Install

Idempotent:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("Seurat",   quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2",  quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",    quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("cowplot",  quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("patchwork",quietly = TRUE)) install.packages("patchwork")

stopifnot(packageVersion("Seurat") >= "5.0.0")
```

Attach once:

```r
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(cowplot); library(patchwork)
})
```

For the Azimuth pre-built reference path (optional), see
`references/azimuth_alternative.md` — different install.

## Decisions to surface up front

1. **Reference choice** — a user-supplied labeled Seurat object (the
   generic case) or a pre-built Azimuth reference (PBMC, lung, kidney,
   etc., see `references/azimuth_alternative.md`). The reference defines
   the label vocabulary the query will inherit; pick one whose tissue /
   organism / annotation depth matches the query.
2. **Label column on the reference** — which metadata column on the
   reference holds the cell-type labels to transfer (commonly `celltype`,
   `cell_type`, `seurat_annotations`). Confirm by
   `colnames(reference[[]])`.
3. **Reference reduction (`reference.reduction`)** — the PCA / integrated
   reduction on the reference that the query will be projected through.
   Usually `"pca"` (the reference's PCA); if the reference was integrated,
   may be `"integrated.cca"` or similar. **The reference's UMAP must have
   been built on this same reduction with `return.model = TRUE`** — Step
   1 checks this.
4. **Number of PCs (`dims`)** — vignette uses `1:30`. Same `dims` flow
   through `FindTransferAnchors`, `TransferData`, and `MapQuery`.
5. **Normalization** — log-normalized query/reference is the default. If
   the reference was built with SCT, set `normalization.method = "SCT"`
   in `FindTransferAnchors` (and normalize the query with `SCTransform`
   first). For details + verified flow see `references/mapping_internals.md`.
6. **Path** — Path A (TransferData only — labels, no UMAP projection) vs
   Path B (MapQuery — labels + UMAP projection). Path B is the standard;
   pick A only when the reference lacks a saved UMAP model.
7. **Mode** — INTERACTIVE (default) when this is the primary analysis;
   BATCH when invoked as one of many sub-runs in a larger plan (e.g.
   map a cohort of N query samples against one reference). The
   orchestrator passes `args="batch"`. See the `## Batch variant`
   section below.

Show the user these figures as the analysis proceeds:
- `umap_reference_celltype.png` — reference UMAP, colored by cell type
- `umap_query_predicted.png` — query projected onto the reference UMAP,
  colored by predicted cell type
- `umap_reference_query_side_by_side.png` — side-by-side comparison

---

## Step 1 — Load the reference and the query; verify the reference has a saved UMAP model

The reference is a Seurat object with (a) raw counts + normalized data,
(b) a PCA reduction, (c) a UMAP reduction built with `return.model =
TRUE`, and (d) a cell-type label column in `meta.data`. Without
`return.model = TRUE` the `MapQuery` UMAP projection in Step 4 cannot
run — it needs the stored `uwot` model to project query cells.

```r
# Reference and query are Seurat objects on disk. Substitute the user's paths.
reference <- readRDS("/path/to/labeled_reference.rds")
query     <- readRDS("/path/to/unlabeled_query.rds")

# 1. The reference needs PCA and a saved-model UMAP.
stopifnot("pca"  %in% Reductions(reference))
stopifnot("umap" %in% Reductions(reference))

# 2. The UMAP model must have been saved via return.model = TRUE.
#    Without it MapQuery cannot project. If this fails, re-run RunUMAP on
#    the reference with return.model = TRUE BEFORE proceeding.
ref_umap_model <- reference[["umap"]]@misc$model
if (is.null(ref_umap_model)) {
  stop("Reference UMAP has no saved model. Re-run on the reference:\n",
       "  reference <- RunUMAP(reference, dims = 1:30, reduction = 'pca',\n",
       "                       return.model = TRUE)\n",
       "and save it back.")
}

# 3. Identify the label column to transfer. Confirm it exists.
LABEL_COL <- "celltype"   # CASE-SENSITIVE; check colnames(reference[[]]) first
stopifnot(LABEL_COL %in% colnames(reference[[]]))

# 4. Reference must be normalized (data layer present).
#    If the reference has counts only, normalize it once now.
if (!"data" %in% Layers(reference[["RNA"]])) {
  reference <- NormalizeData(reference, verbose = FALSE)
}

# 5. Query must also be normalized; reference mapping does not do this for you.
if (!"data" %in% Layers(query[["RNA"]])) {
  query <- NormalizeData(query, verbose = FALSE)
}

cat(sprintf("Reference: %d cells x %d genes | %d cell types in '%s'\n",
            ncol(reference), nrow(reference),
            length(unique(reference[[LABEL_COL, drop = TRUE]])), LABEL_COL))
cat(sprintf("Query:     %d cells x %d genes\n", ncol(query), nrow(query)))
print(table(reference[[LABEL_COL, drop = TRUE]]))
```

### Variant: Azimuth pre-built reference — use INSTEAD of the canonical block above when the user names Azimuth or wants a pre-built tissue atlas

The full Azimuth path lives in `references/azimuth_alternative.md`.
The one-line replacement for Steps 1–4:

```r
# library(Azimuth)
# query <- RunAzimuth(query, reference = "pbmcref")   # or lungref / kidneyref / ...
```

After `RunAzimuth`, jump to Step 5 — Azimuth has already produced the
projected `ref.umap` reduction and the `predicted.celltype.l1/l2/l3`
columns on the query.

**Pitfall — reference UMAP without `return.model`.** This is the single
most common failure on the first run. The reference's `RunUMAP` must
have been called with `return.model = TRUE`. If the reference was built
by someone else and you cannot retrain the UMAP, you can still run
label transfer (`FindTransferAnchors` + `TransferData`) — just skip the
UMAP projection (Path A in Step 3).

**Pitfall — gene-name mismatch between query and reference.**
`FindTransferAnchors` finds shared features automatically, but if
symbol-vs-Ensembl conventions differ between the two objects (one uses
`CD3D`, the other uses `ENSG…`), the intersection will be small or
empty and anchors will not be found. Run
`length(intersect(rownames(reference), rownames(query)))` before Step
2; a sane intersection is >5,000 shared features.

**Report:** reference cells/genes, label vocabulary (cell types + counts
per type), query cells/genes, shared feature count, and confirmation
that the reference UMAP has a saved model.

---

## Step 2 — Find transfer anchors between reference and query

`FindTransferAnchors` learns the correspondence between query and
reference cells through the reference's PCA. The `reference.reduction`
arg names the reduction; `dims` controls how many PCs feed the search.

```r
DIMS_USE <- 1:30   # match the reference's chosen dim count

# Verified against formals(FindTransferAnchors) — Seurat 5.5.0.
# Key args:
#   reference            the labeled reference Seurat object
#   query                the unlabeled query Seurat object
#   reference.reduction  name of the reduction on reference to project through
#                        ("pca" in the vignette; could be "integrated.cca" for
#                        an integrated reference)
#   normalization.method "LogNormalize" (default) or "SCT" — match the reference
#   dims                 PCs to use; vignette uses 1:30
anchors <- FindTransferAnchors(
  reference            = reference,
  query                = query,
  dims                 = DIMS_USE,
  reference.reduction  = "pca",
  normalization.method = "LogNormalize",
  verbose              = FALSE
)
```

**Pitfall — `reference.reduction` is the contract.** Whatever reduction
name you pass must exist on the reference (`"pca"` in the vignette;
could be `"integrated.cca"` for an integrated reference). Verify with
`Reductions(reference)`.

**Pitfall — `normalization.method` must match the reference.**
`"LogNormalize"` (default) for log-normalized references; `"SCT"` for
SCTransform-built references. The case matters: `"SCT"` not `"sct"`. If
the reference was SCT-built, `SCTransform(query)` must run before
`FindTransferAnchors`.

**Report:** number of anchors found (`nrow(anchors@anchors)`), number
of shared genes used, `reference.reduction` name. A reasonable anchor
count is >500 for a typical 5–20k-cell query; very low anchor counts
(<100) signal a gene-name mismatch, a reference-vs-query technology
mismatch, or non-overlapping cell types.

For arg-by-arg semantics (`k.anchor`, `k.filter`, `k.score`,
`max.features`, `npcs`, `approx.pca`) and what each does to the anchor
set, read `references/mapping_internals.md`.

---

## Step 3 — Path A: label transfer only (no UMAP projection)

Use `TransferData` if you want predicted labels without projecting
query cells into the reference UMAP space. This is the lightweight
path. **Skip to Step 4 if the reference UMAP has `return.model = TRUE`
— Path B (MapQuery) gives you Path A's output PLUS the UMAP
projection.**

**Path A and Path B are exclusive.** Don't `AddMetaData` from
`TransferData` into `query` and then also call `MapQuery(query, …)` on
the same object — the weights matrix `TransferData` leaves on the query
is captured into `MapQuery`'s parallel closure and trips `future`'s
`future.globals.maxSize` ceiling (default 500 MiB), with an obscure
error: `The total size of the N globals exported for future expression
… exceeds the maximum allowed size 500.00 MiB`. Either pick ONE path,
run the two on independent query copies, or pre-raise the limit:
`options(future.globals.maxSize = 2 * 1024^3)`.

```r
# refdata is the per-reference-cell label vector (or a named list of multiple
# vectors to transfer several columns in one call).
predictions <- TransferData(
  anchorset = anchors,
  refdata   = reference[[LABEL_COL, drop = TRUE]],
  dims      = DIMS_USE
)

# predictions is a data.frame indexed by query cell barcode with columns:
#   predicted.id           — the top-scoring label per cell
#   prediction.score.<lbl> — per-cell score for every label in the reference
#                            (one column per cell type)
#   prediction.score.max   — the max score across labels (per-cell confidence)

# Add to the query object's meta.data.
query <- AddMetaData(query, metadata = predictions)

cat(sprintf("Transferred labels: %d query cells x %d unique predicted types\n",
            ncol(query), length(unique(query$predicted.id))))
print(table(predicted = query$predicted.id))
cat(sprintf("Median prediction.score.max: %.3f\n",
            median(query$prediction.score.max, na.rm = TRUE)))
```

**Report:** predicted-label distribution, median `prediction.score.max`
(confidence), and the count of low-confidence cells (e.g.
`prediction.score.max < 0.5`). A bimodal score histogram — a confident
peak near 1.0 + a low-confidence tail — is healthy; a flat distribution
centered ~0.3 means the reference doesn't cover the query's cell types
well.

---

## Step 4 — Path B: full mapping with UMAP projection (MapQuery)

`MapQuery` wraps `TransferData` + `IntegrateEmbeddings` + `ProjectUMAP`
into one call. It produces the predicted labels (as in Step 3) AND
projects query cells into the reference's PCA + UMAP space. Use this
when you want the query overlaid on the reference UMAP.

```r
# Verified against formals(MapQuery) — Seurat 5.5.0.
# refdata as a NAMED LIST — names become column prefixes on query meta.data.
# Here we transfer the celltype column; on the query it lands as
# "predicted.celltype" (and "predicted.celltype.score" for the max-score).
query <- MapQuery(
  anchorset           = anchors,
  reference           = reference,
  query               = query,
  refdata             = list(celltype = LABEL_COL),
  reference.reduction = "pca",
  reduction.model     = "umap",       # name of the UMAP on the reference
                                       # with return.model = TRUE
  verbose             = FALSE
)
```

After `MapQuery`:

- `query$predicted.celltype` — top-scoring label per cell (one column
  per entry in `refdata`; the column name comes from the LIST NAME, not
  the reference column name).
- `query$predicted.celltype.score` — per-cell confidence for that label.
- `query[["ref.pca"]]` — query cells projected into the reference's
  PCA space.
- `query[["ref.umap"]]` — query cells projected onto the reference's
  UMAP via the saved UMAP model.

```r
# Sanity-check the projection landed.
stopifnot("ref.umap" %in% Reductions(query))
stopifnot("ref.pca"  %in% Reductions(query))
stopifnot("predicted.celltype" %in% colnames(query[[]]))

cat(sprintf("Mapped %d query cells | %d unique predicted celltypes\n",
            ncol(query), length(unique(query$predicted.celltype))))
print(table(predicted = query$predicted.celltype))
cat(sprintf("Median predicted.celltype.score: %.3f\n",
            median(query$predicted.celltype.score, na.rm = TRUE)))
```

**Pitfall — `reduction.model` names a UMAP that lives on the
REFERENCE.** The string `"umap"` here refers to `reference[["umap"]]`,
not to anything on the query. If the reference's UMAP was saved under a
different name (e.g. `"umap.integrated"`), pass that name.

**Pitfall — `refdata` as a named list controls the OUTPUT column name
on the query.** `refdata = list(celltype = LABEL_COL)` produces
`query$predicted.celltype`. If you pass `refdata = list(cell_type =
LABEL_COL)`, the column is `query$predicted.cell_type`. The LIST NAME
wins, not the reference column name. Confirm by `colnames(query[[]])`
right after the call.

**Report:** predicted-label distribution, median
`predicted.celltype.score`, fraction of cells with score < 0.5
(low-confidence), and that the `ref.umap` reduction was created.

For score interpretation and what to do about low-confidence cells, see
`references/mapping_internals.md`.

---

## Step 5 — Visualize: reference, query, and the two side-by-side

Three figures: the reference UMAP colored by cell type (so the user
sees the label vocabulary), the query on its projected reference-UMAP
colored by predicted label, and the two side-by-side for direct
comparison.

```r
# 1. Reference UMAP, colored by the labels you transferred from.
p_ref <- DimPlot(reference, reduction = "umap", group.by = LABEL_COL,
                 label = TRUE, repel = TRUE, pt.size = 0.4) +
  ggtitle(sprintf("Reference UMAP — %s (%d cells, %d types)",
                  LABEL_COL, ncol(reference),
                  length(unique(reference[[LABEL_COL, drop = TRUE]])))) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4)) +
  NoLegend() +              # AFTER theme_cowplot — cowplot 1.2.0 overrides legend.position
  coord_fixed()
for (k in seq_along(p_ref$layers)) {
  if (inherits(p_ref$layers[[k]]$geom, "GeomPoint")) {
    p_ref$layers[[k]]$aes_params$alpha <- 0.6
  }
}
ggsave("umap_reference_celltype.png", p_ref,
       width = 7, height = 6, units = "in", dpi = 120, bg = "white")

# 2. Query, projected onto the reference UMAP, colored by predicted label.
p_qry <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype",
                 label = TRUE, repel = TRUE, pt.size = 0.4) +
  ggtitle(sprintf("Query projected onto reference UMAP — %d cells, %d predicted types",
                  ncol(query),
                  length(unique(query$predicted.celltype)))) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4)) +
  NoLegend() +              # AFTER theme_cowplot — cowplot 1.2.0 overrides legend.position
  coord_fixed()
for (k in seq_along(p_qry$layers)) {
  if (inherits(p_qry$layers[[k]]$geom, "GeomPoint")) {
    p_qry$layers[[k]]$aes_params$alpha <- 0.6
  }
}
ggsave("umap_query_predicted.png", p_qry,
       width = 7, height = 6, units = "in", dpi = 120, bg = "white")

# 3. Side by side (patchwork). Same axes scale via coord_fixed() upstream.
p_side <- (p_ref | p_qry) +
  plot_annotation(title = "Reference (left) vs query mapped onto reference (right)",
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))
ggsave("umap_reference_query_side_by_side.png", p_side,
       width = 14, height = 6, units = "in", dpi = 120, bg = "white")
```

**Report:** predicted labels with very few query cells (e.g. <5) are
likely artifacts. Cell-type fractions on reference vs query — large
shifts can be real (different tissue composition) or label-transfer
error.

For shared figure style (palette, alpha-poke, ggsave defaults) see
`references/figure_style.md`.

---

## Step 6 — Persist results: predictions CSV + mapped query RDS

```r
# Write a per-cell predictions CSV the user can open in Excel / pandas.
pred_df <- data.frame(
  barcode                   = colnames(query),
  predicted.celltype        = query$predicted.celltype,
  predicted.celltype.score  = query$predicted.celltype.score,
  row.names                 = NULL
)
write.csv(pred_df, "query_predicted_labels.csv", row.names = FALSE)

# Save the full mapped query (predictions + ref.pca + ref.umap baked in).
saveRDS(query, "query_mapped.rds")

# Verify
cat(sprintf("Wrote query_predicted_labels.csv (%d rows)\n", nrow(pred_df)))
cat(sprintf("Wrote query_mapped.rds (%.1f MB)\n",
            file.info("query_mapped.rds")$size / 1e6))
```

---

## Batch variant — use INSTEAD of Steps 1–6 when invoked with args="batch"

When an orchestrator maps N query datasets against one reference (e.g.
a cohort study, atlas projection over many GSMs), the per-query
figures and reports become noise. Branch on `$ARGUMENTS` at the top:

- Skip all three UMAPs.
- Skip the per-step `Report` footers and per-call confidence summaries.
- Still write the canonical artifacts the orchestrator consumes:
  `query_mapped.rds` and `query_predicted_labels.csv`.
- Print ONE final summary line: `"batch ok | <N> query cells | <K>
  predicted types | median score <s>"`.

```r
suppressPackageStartupMessages({ library(Seurat) })
stopifnot(packageVersion("Seurat") >= "5.0.0")

# Inputs (orchestrator-provided, or literal for solo args="batch"):
#   reference (already loaded with PCA + UMAP w/ return.model = TRUE)
#   query     (already loaded + normalized)
#   LABEL_COL (label column on reference)
LABEL_COL <- "celltype"

# Reference checks (silent fail if the reference is wrong shape)
stopifnot("pca"  %in% Reductions(reference),
          "umap" %in% Reductions(reference),
          !is.null(reference[["umap"]]@misc$model),
          LABEL_COL %in% colnames(reference[[]]))

DIMS_USE <- 1:30

anchors <- FindTransferAnchors(
  reference            = reference, query = query, dims = DIMS_USE,
  reference.reduction  = "pca",
  normalization.method = "LogNormalize", verbose = FALSE
)
query <- MapQuery(
  anchorset           = anchors, reference = reference, query = query,
  refdata             = list(celltype = LABEL_COL),
  reference.reduction = "pca", reduction.model = "umap",
  verbose             = FALSE
)

# Persist canonical artifacts
write.csv(
  data.frame(barcode = colnames(query),
             predicted.celltype = query$predicted.celltype,
             predicted.celltype.score = query$predicted.celltype.score),
  "query_predicted_labels.csv", row.names = FALSE
)
saveRDS(query, "query_mapped.rds")

cat(sprintf("batch ok | %d query cells | %d predicted types | median score %.3f\n",
            ncol(query),
            length(unique(query$predicted.celltype)),
            median(query$predicted.celltype.score, na.rm = TRUE)))
```

What batch mode keeps vs drops:
- KEEPS: `query_mapped.rds`, `query_predicted_labels.csv`, the one summary line.
- DROPS: every `ggsave`, per-step `Report` footer, the confidence histogram.

---

## Final response checklist

Summarize, in this order:

- reference object: cells, genes, label column used, number of label
  categories, whether UMAP model was saved (`return.model = TRUE`)
- query object: cells, genes, normalization status
- shared genes between reference and query (anchor-finding feasibility)
- anchor count, `reference.reduction`, `dims` used,
  `normalization.method`
- path taken: TransferData only (labels) vs MapQuery (labels + UMAP)
- predicted-label distribution (top types + their query-cell counts)
- median `predicted.celltype.score` (or `prediction.score.max` for the
  TransferData path) — overall confidence
- fraction of low-confidence cells (e.g. `score < 0.5`)
- figures written (3 if MapQuery; 0 if batch) and their paths
- saved artifacts (`query_mapped.rds`, `query_predicted_labels.csv`)
- caveats: gene-name mismatches, label-vocabulary fit, technology /
  tissue mismatch between reference and query, low-confidence cells
  the user should review

---

## See also

- `seurat-integration` — when you want a symmetric joint embedding
  instead of projecting query onto a fixed reference (re-learns from
  both). Use this if the query is large enough to influence the
  embedding or if reference and query come from different protocols and
  you want them to inform each other.
- `seurat-scrna-v2` — single-sample QC + clustering. Run this on the
  query BEFORE reference mapping to confirm the query is clean; mapping
  bad-QC cells gives confidently-wrong predictions for dying / empty
  droplets.
- `annotate-celltype-scrna` — manual marker-based annotation as an
  alternative when no labeled reference exists.
- Azimuth — see `references/azimuth_alternative.md` for pre-built
  tissue references that wrap this whole flow into one call.
