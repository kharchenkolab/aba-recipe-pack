---
name: seurat-multimodal-reference-mapping
description: Map a scRNA-seq (or CITE-seq) query onto a pre-built multimodal reference in R/Seurat v5 — find SCT-flavoured anchors against the reference's supervised PCA (spca) with FindTransferAnchors(), then MapQuery() to transfer cell-type labels at MULTIPLE granularities (l1/l2/l3) AND impute predicted protein levels, all projected onto the reference's pre-computed WNN UMAP. Produces a query Seurat object with `predicted.celltype.l*` columns, per-cell `.score` confidences, a `predicted_ADT` assay, and a `ref.umap` reduction.
when_to_use: A scRNA-seq or CITE-seq query Seurat object AND a pre-built multimodal reference (the Azimuth PBMC 162k reference, the bone-marrow CITE-seq reference, or a custom WNN reference with `spca` + `wnn.umap`) when the user wants cell-type labels at multiple resolutions (broad → fine) WITHOUT running unsupervised clustering. Use THIS when the reference exists and is appropriate to the query's tissue. For RNA-only reference mapping (no spca, no predicted ADT), use `seurat-reference-mapping`. For building a NEW reference from your own data, run `seurat-wnn-multimodal` first (then export it as a reference via `RunSPCA` + the existing `wnn.umap`).
invocation: interactive+batch
requires_tools: [run_r]
capabilities_needed: [Seurat]
keywords: [Seurat, reference mapping, Azimuth, multimodal reference, FindTransferAnchors, MapQuery, supervised PCA, spca, reference.reduction, predicted celltype, celltype.l1, celltype.l2, predicted_ADT, ProjectUMAP, IntegrateEmbeddings, TransferData, ref.umap, R, v5]
produces: [umap_ref_celltype_l1.png, umap_ref_celltype_l2.png, predicted_score_hist.png, predicted_adt_featureplot.png, query_mapped.rds]
domain: genomics
source: "Seurat multimodal reference mapping vignette (Satija Lab) — satijalab.org/seurat/articles/multimodal_reference_mapping. Reference building: Hao et al. Cell 2021 (PBMC 162k CITE-seq atlas)."
---

# Multimodal reference mapping (Seurat v5)

A "multimodal reference" is a pre-built Seurat object that carries (a) a
**supervised PCA** (`spca`) — a PCA on the RNA assay trained to preserve
cell-type structure jointly with the protein modality — and (b) a
**WNN UMAP** (`wnn.umap`) — the joint embedding the reference was clustered
on. Mapping a query onto it transfers labels and embeddings WITHOUT running
unsupervised clustering on the query at all. Labels come at multiple
granularities (`celltype.l1`, `celltype.l2`, `celltype.l3` for the PBMC
reference), and the reference can also impute the ADT panel onto an RNA-only
query (the `predicted_ADT` assay).

The query can be RNA-only or CITE-seq; both go through the same anchor
finder. The reference must already exist on disk — this recipe does NOT
build references (use `seurat-wnn-multimodal` first if you need to).

Pin **Seurat v5** (recipe targets `Seurat 5.5.0`; the vignette was originally
drafted in v4 but the v5 API is backward-compatible for `FindTransferAnchors`
/ `MapQuery`).

## Bundled references — load on demand

This recipe is self-contained for the standard workflow. For deeper detail
on any aspect, load the matching reference file with `read_file` ONLY when
the task needs it — don't pre-load everything:

- `references/multimodal_mapping_internals.md` — what FindTransferAnchors +
  MapQuery actually do (TransferData + IntegrateEmbeddings + ProjectUMAP),
  the supervised-PCA (spca) construction, SCT vs LogNormalize anchor
  finding, dims-range guidance.
- `references/azimuth_references.md` — catalogue of pre-built references
  (PBMC 162k, bone marrow, fetal, kidney, …), the `Azimuth::RunAzimuth`
  wrapper, how to obtain reference RDS files, naming conventions per
  reference (the dash-suffixed ADT names on PBMC 162k, etc.).
- `references/prediction_qc.md` — score interpretation, when to flag
  low-confidence cells, bimodal vs unimodal score distributions, what
  ambiguity means biologically, the failure mode where a tissue-mismatched
  reference yields confident wrong labels.
- `references/figure_style.md` — Seurat plot styling shared with the rest
  of the collection (theme, palettes, alpha-poke, save dimensions).

## Install

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("Seurat",      quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("ggplot2",     quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",       quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("cowplot",     quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("patchwork",   quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("glmGamPoi",   quietly = TRUE)) {                # speeds up SCTransform
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("glmGamPoi", update = FALSE, ask = FALSE)
}
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

1. **Reference choice** — the reference MUST match the query's tissue (PBMC
   reference → PBMC query, BM reference → bone-marrow query). A
   tissue-mismatched reference will confidently assign wrong cell types
   with high `prediction.score` — the scores reflect anchor density, not
   correctness. See `references/azimuth_references.md` for the catalogue
   and the per-reference tissue scope.
2. **Normalization method match** — `FindTransferAnchors(normalization.method =
   "SCT")` requires the QUERY to be `SCTransform`-normalized AND the
   reference to be `SCT`-normalized. Most modern references use SCT; verify
   via `Assays(reference)`. If the reference is `LogNormalize`d (older
   references, custom-built ones via `NormalizeData`), pass
   `normalization.method = "LogNormalize"` and `NormalizeData` the query
   instead.
3. **`dims` range** — the supervised PCA on PBMC 162k uses 50 components
   (`dims = 1:50`); other references may differ. Use the reference's full
   `spca` width unless you have a reason to truncate. See
   `references/multimodal_mapping_internals.md`.
4. **Which labels to transfer** — pass a NAMED LIST to `refdata=` listing
   every metadata column you want imputed. The keys become the new
   `predicted.<key>` columns; the special key `predicted_ADT` (or any name
   pointing at an assay name) imputes that assay onto the query.
5. **`prediction.score` cutoff** — every transferred label carries a
   `predicted.<label>.score` column (max prediction probability). The
   common convention is to flag cells with `score < 0.5` as
   low-confidence. There is no universal cutoff — see
   `references/prediction_qc.md` for distribution-aware guidance.
6. **Mode** — INTERACTIVE (default) when this is the primary analysis;
   BATCH when invoked over many query samples in a larger plan. The
   orchestrator passes `args="batch"`; the agent declares the mode in
   `present_plan`.

Figures to show as the analysis proceeds:
- `umap_ref_celltype_l1.png`
- `umap_ref_celltype_l2.png`
- `predicted_score_hist.png`
- `predicted_adt_featureplot.png`

---

## Step 1 — Load the reference and inspect it

The reference is a Seurat object with two non-default reductions: `spca`
(supervised PCA on the RNA assay) and `wnn.umap` (the joint UMAP the
reference was clustered on). Both are REQUIRED by `MapQuery`.

```r
# The PBMC 162k reference (or whichever multimodal reference you have).
reference <- readRDS("/path/to/pbmc_multimodal_2023.rds")

# Inspect what the reference carries — this drives the refdata= list in Step 3.
print(Assays(reference))                  # expect e.g. c("SCT", "ADT")
print(Reductions(reference))              # MUST include "spca" and "wnn.umap"
print(head(colnames(reference@meta.data)))
# Look for the celltype.l1 / l2 / l3 columns (or whatever the reference names them).
celltype_cols <- grep("^celltype\\.", colnames(reference@meta.data), value = TRUE)
print(celltype_cols)

stopifnot("spca"     %in% Reductions(reference))
stopifnot("wnn.umap" %in% Reductions(reference))
```

**Report:** reference name + size (cells, genes, proteins), the celltype
columns available (l1/l2/l3 etc.), and which assays carry which modality.

For where to obtain the reference + alternatives, read
`references/azimuth_references.md`.

### Variant: use Azimuth's RunAzimuth wrapper — use INSTEAD of this whole recipe when you don't need anchor-finding control

```r
# install.packages("BiocManager"); BiocManager::install("satijalab/azimuth")
# library(Azimuth)
# query <- RunAzimuth(query, reference = "pbmcref")
# RunAzimuth fetches the reference, runs FindTransferAnchors + MapQuery,
# and returns the query with predicted.celltype.l1/l2/l3 already populated.
# Use this when the goal is just labels; use THIS recipe when you need to
# inspect anchors, tune dims, or filter on prediction.score before MapQuery.
```

---

## Step 2 — Prepare the query (SCTransform-normalize to match the reference)

The reference uses `SCT` normalization; the query must too. If the query is
already in a Seurat object with raw `RNA` counts (the typical case), run
`SCTransform` directly.

```r
# query is the Seurat object the user wants mapped. Read it in from upstream:
query <- readRDS("/path/to/query.rds")        # e.g. a CITE-seq object from seurat-cite-seq
DefaultAssay(query) <- "RNA"

# SCTransform must run on raw counts. method="glmGamPoi" is 5-10x faster than
# the poisson default; install glmGamPoi (Bioc) — already in the Install block.
query <- SCTransform(query,
                     method  = "glmGamPoi",
                     verbose = FALSE)
```

### Variant: LogNormalize query — use INSTEAD of the SCTransform block above when the reference is LogNormalize-d

```r
# Verify first: if rownames(reference[["SCT"]]) is empty / no SCT assay, the
# reference is LogNormalize-d. Many custom-built references and older
# vignettes go this route.
DefaultAssay(query) <- "RNA"
query <- NormalizeData(query, verbose = FALSE)
# In Step 3, also pass normalization.method = "LogNormalize" to FindTransferAnchors.
```

**Report:** query cell count, query gene overlap with the reference's
analysis features. Low overlap = wrong species or wrong reference choice;
should be >80% for a sensible mapping.

```r
# Use whichever assay the reference uses ("SCT" or "RNA"); the principle is the same.
ref_assay  <- if ("SCT" %in% Assays(reference)) "SCT" else "RNA"
qry_assay  <- if ("SCT" %in% Assays(query))     "SCT" else "RNA"
ref_features <- rownames(reference[[ref_assay]])
qry_features <- rownames(query[[qry_assay]])
ovl <- length(intersect(ref_features, qry_features))
cat(sprintf("Query: %d cells. Feature overlap with reference %s: %d/%d (%.1f%%)\n",
            ncol(query), ref_assay, ovl, length(ref_features),
            100 * ovl / length(ref_features)))
if (ovl / length(ref_features) < 0.8) {
  warning("Low feature overlap with reference — check that species + gene-symbol convention match.")
}
```

---

## Step 3 — Find transfer anchors against the reference's supervised PCA

`FindTransferAnchors` with `reference.reduction = "spca"` and
`normalization.method = "SCT"` is the multimodal-mapping mode. The
supervised PCA was trained on the reference to preserve cell-type structure
discriminatively — anchors found in `spca` space are higher quality than
anchors found in a vanilla PCA.

```r
# normalization.method = "SCT"  -> both reference and query must be SCT-normalized.
#                                  Pass "LogNormalize" instead if the reference is LogNormalize-d.
# reference.reduction = "spca"  -> use the reference's supervised PCA (not standard PCA).
# dims = 1:50                   -> use all 50 supervised components (vignette default for PBMC 162k).
# recompute.residuals = TRUE    -> SCT-only: recompute SCT residuals on the query using the
#                                  reference's model when SCTransform versions differ. Harmless but
#                                  IRRELEVANT under normalization.method = "LogNormalize".
anchors <- FindTransferAnchors(
  reference            = reference,
  query                = query,
  normalization.method = "SCT",          # swap to "LogNormalize" for log-normalized references
  reference.reduction  = "spca",
  dims                 = 1:50,
  recompute.residuals  = TRUE,
  verbose              = FALSE
)

cat(sprintf("Found %d anchors (query %d cells, reference %d cells)\n",
            nrow(anchors@anchors), ncol(query), ncol(reference)))
```

If the anchor count is very low (<1% of query cells), the reference is a
poor match for the query tissue — stop and reconsider.

For SCT vs LogNormalize semantics, the spca construction, and how anchors
are scored, read `references/multimodal_mapping_internals.md`.

---

## Step 4 — MapQuery: transfer labels + impute proteins + project onto ref UMAP

`MapQuery()` is the umbrella call. It runs three things internally:
- `TransferData()` for each entry in `refdata=` (categorical labels or
  continuous assays);
- `IntegrateEmbeddings()` to project the query's expression into the
  reference's `spca` space;
- `ProjectUMAP()` to place the query cells on the reference's `wnn.umap`.

```r
# refdata is a NAMED LIST. Keys become the new "predicted.<key>" columns on
# the query. Values are either:
#   - a string naming a reference metadata column (e.g. "celltype.l1") -> label transfer
#   - a string naming a reference ASSAY (e.g. "ADT")                    -> assay imputation
#     (the assay name on the query side will be the LIST KEY, e.g. "predicted_ADT").
query <- MapQuery(
  anchorset           = anchors,
  query               = query,
  reference           = reference,
  refdata             = list(
    celltype.l1   = "celltype.l1",
    celltype.l2   = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model     = "wnn.umap"
)
```

Outputs added to `query`:
- `query$predicted.celltype.l1`  + `query$predicted.celltype.l1.score`
- `query$predicted.celltype.l2`  + `query$predicted.celltype.l2.score`
- `query[["predicted_ADT"]]` — a new assay with imputed protein levels
- `query[["ref.umap"]]` — the projection of the query onto the reference's WNN UMAP

```r
# Verify the artifacts landed:
stopifnot("predicted.celltype.l1"      %in% colnames(query@meta.data))
stopifnot("predicted.celltype.l1.score" %in% colnames(query@meta.data))
stopifnot("predicted_ADT" %in% Assays(query))
stopifnot("ref.umap"      %in% Reductions(query))
```

**Report:** counts per `predicted.celltype.l1` (and `.l2`), median
`prediction.score` per label, fraction of cells below 0.5 score (the
low-confidence tail).

```r
cat("predicted.celltype.l1 counts:\n")
print(table(query$predicted.celltype.l1))

cat("\nlow-confidence fraction (l1 score < 0.5):\n")
cat(sprintf("%.1f%%\n",
            100 * mean(query$predicted.celltype.l1.score < 0.5)))
```

For what each MapQuery sub-step does and when to call them individually,
read `references/multimodal_mapping_internals.md`.

---

## Step 5 — Visualize labels at multiple granularities

Plot the query cells on the reference UMAP, colored by the transferred
labels. l1 (broad) and l2 (fine) side-by-side gives the user the resolution
spectrum.

```r
make_label_plot <- function(query, label, title) {
  p <- DimPlot(query, reduction = "ref.umap",
               group.by = label,
               label = TRUE, label.size = 3, repel = TRUE,
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

p_l1 <- make_label_plot(query, "predicted.celltype.l1",
                        sprintf("Predicted celltype.l1 (%d cells)", ncol(query)))
p_l2 <- make_label_plot(query, "predicted.celltype.l2",
                        sprintf("Predicted celltype.l2 (%d cells)", ncol(query)))

ggsave("umap_ref_celltype_l1.png", p_l1,
       width = 8, height = 6.5, units = "in", dpi = 120, bg = "white")
ggsave("umap_ref_celltype_l2.png", p_l2,
       width = 9, height = 6.5, units = "in", dpi = 120, bg = "white")
```

Prediction-score distribution — a histogram per l1 label catches the
low-confidence tail:

```r
df_scores <- data.frame(
  celltype = query$predicted.celltype.l1,
  score    = query$predicted.celltype.l1.score
)

p_scores <- ggplot(df_scores, aes(x = score, fill = celltype)) +
  geom_histogram(binwidth = 0.02, alpha = 0.85, colour = "black",
                 linewidth = 0.2) +
  geom_vline(xintercept = 0.5, colour = "red",
             linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ celltype, scales = "free_y") +
  # scale_fill_discrete handles any number of celltypes; Set2 caps at 8 and
  # silently recycles on tissues with more populations.
  scale_fill_discrete(guide = "none") +
  labs(title = "prediction.score by celltype.l1 (red line = 0.5 confidence threshold)",
       x = "prediction.score", y = "cells") +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor   = element_blank())

ggsave("predicted_score_hist.png", p_scores,
       width = 11, height = 7, units = "in", dpi = 120, bg = "white")
```

**Assess and report:** which l1 labels have a clean score distribution
(>0.8 median, narrow) — those mappings are robust. Which labels have a
bimodal score distribution — those are populations where the query is
ambiguous between two reference cell types. Cells below the 0.5 line: do
they cluster in one region of `ref.umap` (a population the reference
doesn't have) or scatter randomly (noise / doublets)?

For score-interpretation depth, read `references/prediction_qc.md`.

---

## Step 6 — Visualize predicted protein levels (RNA-only query)

If the query was RNA-only, `predicted_ADT` is the new assay carrying imputed
protein levels. Plot them on the same `ref.umap` to show the lineage
discrimination the protein panel adds:

```r
DefaultAssay(query) <- "predicted_ADT"

# Replace with the proteins your reference's ADT panel actually carries —
# check `rownames(reference[["ADT"]])` for the canonical names. The PBMC
# 162k reference uses dash-suffixed names like "CD3-1", "CD45RA", "IgD"
# (Hao et al. 2021 conventions); other references (bmcite, BM CITE) use
# plain names like "CD3", "CD56". See references/azimuth_references.md.
proteins_show <- intersect(c("CD3-1", "CD4-1", "CD8a", "CD19", "CD14", "CD56-1"),
                           rownames(query[["predicted_ADT"]]))
if (length(proteins_show) == 0) {
  # Fallback: take the first 6 features so the recipe doesn't crash on a
  # different panel; the user will rerun with the right names.
  proteins_show <- head(rownames(query[["predicted_ADT"]]), 6)
}

p_pa <- FeaturePlot(query,
                    features = proteins_show,
                    reduction = "ref.umap",
                    order = TRUE, pt.size = 0.3, ncol = 3,
                    cols = c("lightgrey", "#1b7837")) &      # green = imputed protein
  theme_cowplot() &
  theme(plot.title = element_text(size = 11, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.4, "cm"))

for (i in seq_along(p_pa)) {
  pl <- p_pa[[i]]
  if (!is.null(pl$layers)) {
    for (k in seq_along(pl$layers)) {
      if (inherits(pl$layers[[k]]$geom, "GeomPoint")) {
        pl$layers[[k]]$aes_params$alpha <- 0.6
      }
    }
    p_pa[[i]] <- pl
  }
}

n_rows <- ceiling(length(proteins_show) / 3)
ggsave("predicted_adt_featureplot.png", p_pa,
       width = 13, height = max(4, 4 * n_rows),
       units = "in", dpi = 120, bg = "white")

DefaultAssay(query) <- "RNA"   # leave the query on RNA for downstream work
```

> **Caveat: predicted ADT is imputed, not measured.** Use it for *labeling
> sanity-checking* and *visualization*, NOT for downstream differential
> protein expression — the values are inferred from the reference's
> RNA↔ADT correlations, so they collapse onto the reference's structure.
> If the query has a population the reference doesn't, the imputed
> proteins will be wrong for that population (and the prediction.score
> low — flag those cells via `references/prediction_qc.md`).

**Assess and report:** for each major lineage marker, does the imputed
protein gradient track the predicted labels on `ref.umap`? This is the
self-consistency check; large mismatches indicate the query has
populations outside the reference's training distribution.

---

## Step 7 — Save the mapped object

```r
saveRDS(query, file = "query_mapped.rds")
cat(sprintf("Wrote query_mapped.rds (%.1f MB)\n",
            file.info("query_mapped.rds")$size / 1e6))
```

The saved `.rds` carries the original query assays, the `predicted_ADT`
assay, the `ref.umap` reduction, and all `predicted.celltype.*` /
`predicted.celltype.*.score` metadata columns. Downstream work can pick
this up directly.

---

## Batch variant — use INSTEAD of Steps 1–7 when invoked with args="batch"

Branch on `$ARGUMENTS == "batch"` at the top of the body. In batch mode:

- Skip all per-step figures.
- Skip the per-step "Assess and report" footers.
- Still save the canonical `query_mapped.rds`.
- Print ONE final summary line.

```r
suppressPackageStartupMessages({ library(Seurat) })

reference <- readRDS("/path/to/pbmc_multimodal_2023.rds")
query     <- readRDS("/path/to/query.rds")

DefaultAssay(query) <- "RNA"
query <- SCTransform(query, method = "glmGamPoi", verbose = FALSE)

anchors <- FindTransferAnchors(
  reference            = reference,
  query                = query,
  normalization.method = "SCT",
  reference.reduction  = "spca",
  dims                 = 1:50,
  recompute.residuals  = TRUE,
  verbose              = FALSE
)

query <- MapQuery(
  anchorset           = anchors,
  query               = query,
  reference           = reference,
  refdata             = list(
    celltype.l1   = "celltype.l1",
    celltype.l2   = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model     = "wnn.umap"
)

saveRDS(query, file = "query_mapped.rds")
cat(sprintf("batch ok | %d cells | %d l1 labels | median l1 score=%.2f | low-conf=%.1f%% | %.1f MB\n",
            ncol(query),
            length(unique(query$predicted.celltype.l1)),
            median(query$predicted.celltype.l1.score),
            100 * mean(query$predicted.celltype.l1.score < 0.5),
            file.info("query_mapped.rds")$size / 1e6))
```

---

## Final response checklist

Summarize:
- reference used (name, cell count, label columns available)
- query cell count and feature overlap with the reference SCT features
- anchors found (count + as fraction of query cells)
- counts per predicted.celltype.l1 (and l2 if used) — the headline cell-type breakdown
- median prediction.score per l1 label and overall low-confidence fraction (<0.5)
- which labels are robust (high, narrow score distribution) vs ambiguous (bimodal)
- whether predicted_ADT was generated and whether it tracks the labels
- figures shown (umap_ref_celltype_l1, umap_ref_celltype_l2, predicted_score_hist, predicted_adt_featureplot)
- saved files (query_mapped.rds)
- caveats: prediction scores reflect anchor density not biological correctness (tissue mismatch can yield confident wrong labels), predicted_ADT is imputed not measured (don't use for protein DE), populations not in the reference will map to their nearest neighbor with low score — flag these explicitly

---

## See also

- `seurat-reference-mapping` — RNA-only reference mapping (uses standard
  PCA + reference UMAP, no `spca`, no `predicted_ADT`). Switch when the
  reference is RNA-only or the query has no need for protein imputation.
- `seurat-cite-seq` — per-cell protein workup; useful as preprocessing
  BEFORE this recipe when the query is a fresh CITE-seq dataset (so you
  have a clean RNA assay to SCTransform), and useful AFTER for biaxial
  protein inspection of the mapped populations.
- `seurat-wnn-multimodal` — build a custom multimodal reference from
  scratch (run WNN on a curated dataset, then `RunSPCA` to produce the
  `spca` reduction this recipe consumes).
- `seurat-scrna-v2` — the RNA-only single-sample workflow; run for
  quality control of the query BEFORE mapping (low-quality cells map
  confidently to nonsense labels).
