---
name: seurat-de-testing
description: Differential expression testing in R/Seurat (v5) — the practical method-choice recipe. Three paths — Path A FindAllMarkers (cluster markers, Wilcoxon), Path B FindMarkers (two-group test, optional MAST/LR with covariates), Path C AggregateExpression + DESeq2 (PSEUDOBULK for multi-sample condition effects). Includes the explicit ban on running bulk DE tools (DESeq2, edgeR, limma, pydeseq2) directly on per-cell scRNA matrices — those are bulk + pseudobulk tools only.
when_to_use: A processed Seurat object with clusters / cell-type labels and (a) the user wants cluster markers (Wilcoxon FindAllMarkers); OR (b) the user wants to compare two cell groups (FindMarkers); OR (c) the user wants to test a CONDITION effect (stim vs ctrl, treated vs untreated, disease vs healthy) across MULTIPLE samples / donors — that's the pseudobulk path, not per-cell DE. Use THIS when the session is R/Seurat. For Python/scanpy use the scanpy DE recipes; for full bulk RNA-seq DESeq2 see deseq2-r.
avoid_when: The session is Python/scanpy (use the scanpy DE recipes), the input is bulk RNA-seq with one column per sample (use `deseq2-r` / `bulk_rnaseq_de`), or you have only one biological sample for a condition effect (no replication ⇒ no robust DE).
invocation: interactive+batch
requires_tools: [run_r]
capabilities_needed: [Seurat]
keywords: [Seurat, differential expression, DE, FindMarkers, FindAllMarkers, Wilcoxon, MAST, DESeq2, edgeR, pseudobulk, AggregateExpression, condition effect, stim vs ctrl, latent.vars, presto, marker genes, cluster markers, v5, R, p_val_adj, avg_log2FC, pct.1, pct.2]
produces: [cluster_markers.csv, de_volcano.png, de_topgenes_dotplot.png, pseudobulk_de.csv]
domain: genomics
source: "Seurat v5 DE vignette - satijalab.org/seurat/articles/de_vignette (Seurat 5.5.0). Pseudobulk methodology grounded in Squair et al. 2021 (Nat Commun, DOI 10.1038/s41467-021-25960-2) and Soneson & Robinson 2018 (Nat Methods). Per-cell test references: Finak et al. 2015 (MAST, Genome Biol). DESeq2 / edgeR vignettes for pseudobulk design."
---

# Differential expression with R/Seurat (v5)

DE in scRNA-seq is a **method-choice problem** before it is a code problem.
This recipe picks one of three paths up front and runs it. Misclassifying the
question is the single biggest source of wrong DE results in scRNA — more
costly than any parameter choice.

This recipe expects a processed Seurat object with clusters / `seurat_annotations`
already set — i.e. the output of `seurat-scrna-v2` or `seurat-integration` (or
the joint-cluster output of `seurat-rna-atac-integration`). For Path C
(condition effect across samples), the object must also carry a `sample` /
`donor_id` column AND a `condition` column in `meta.data` — one row per cell.

## The three paths

| You want… | Path | Function | When |
|---|---|---|---|
| Markers for every cluster (one vs rest) | **A** | `FindAllMarkers(obj, only.pos = TRUE)` | "What marks each cluster?" — the standard tutorial call. |
| Compare TWO cell groups within one object | **B** | `FindMarkers(obj, ident.1, ident.2)` | "Cluster A vs B", "CD8 vs CD4". Wilcoxon (default) or MAST/LR for per-cell covariates. |
| Condition effect across MULTIPLE samples | **C** | `AggregateExpression(return.seurat=T)` → `FindMarkers(test.use="DESeq2")` | "stim vs ctrl", "responders vs non-responders". **PSEUDOBULK is required** — never run bulk DE tools directly on per-cell counts. |

> **Hard rule — bulk DE tools on per-cell counts is WRONG.** DESeq2, edgeR,
> limma, and pydeseq2 model **bulk** count distributions and assume each
> column is an independent biological replicate. A single cell is NOT a
> replicate of its sample. Running these tools directly on the per-cell
> matrix inflates the false-positive rate **10–100×** (Squair et al. 2021).
> Always pseudobulk first for condition effects. For the long form of this
> argument, citations, and the full anti-pattern catalogue, read
> `references/anti_patterns.md`.

## Bundled references — load on demand

This recipe carries the load-bearing code for all three paths. For the depth
behind any method, load the matching reference file with `read_file` ONLY when
the task needs it — don't pre-load everything:

- `references/per_cell_tests.md` — every `test.use` value in `FindMarkers`
  (wilcox / wilcox_limma / MAST / LR / ROC / t / bimod / negbinom / poisson),
  what each is, when it wins, when it fails, source citations.
- `references/pseudobulk_workflow.md` — `AggregateExpression` internals,
  donor/sample/cell-type structure, the limitation of the
  `FindMarkers(test.use="DESeq2")` wrapper (only ~group design), how to drop
  to native DESeq2 / edgeR for multi-factor / paired / interaction designs.
- `references/anti_patterns.md` — the per-cell-DESeq2 ban with Squair 2021 +
  Soneson & Robinson 2018 references, the pydeseq2-on-AnnData mistake, the
  MAST-with-donor-as-latent.vars mistake, why FPR inflation is the failure
  mode, and the right answer for each.
- `references/output_interpretation.md` — what every column means
  (`p_val`, `p_val_adj`, `avg_log2FC`, `pct.1`, `pct.2`, `myAUC`), how
  Bonferroni adjustment plays out on the scRNA gene-set, biological-coherence
  sanity checks (lineage markers in their own cluster, IFN signature UP in
  STIM).
- `references/figure_style.md` — palette, theme, dpi, alpha-poke, and the
  divergent-scale convention shared across the Seurat recipes.

## Install

`Seurat` is a CRAN binary (already provisioned). `presto` is a GitHub package
that Seurat v5 picks up automatically when present — install it once for a
10–100× speed-up on the Wilcoxon default with identical p-values.

```r
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(cowplot)
})

# Speed-up for the Wilcoxon default — Seurat v5 auto-uses presto when installed.
if (!requireNamespace("presto", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("immunogenomics/presto")
}

# Bioconductor packages — install ONLY for the path you'll use.
# Path B with MAST covariates:
# if (!requireNamespace("MAST",   quietly = TRUE)) BiocManager::install("MAST")
# Path C pseudobulk DESeq2:
# if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
# Path C native edgeR (multi-factor designs):
# if (!requireNamespace("edgeR",  quietly = TRUE)) BiocManager::install("edgeR")
```

## Decisions to surface up front

Surface these in `present_plan` BEFORE picking a test.

1. **The QUESTION** — markers (one cluster vs rest), two-group test (cluster A
   vs cluster B), or condition effect (stim vs ctrl across multiple samples)?
   Each maps to a different path below. See the "Three paths" table at the top.
2. **Sample replication** — for a CONDITION effect, how many biological
   samples do you have per condition? If `n_per_condition < 3`, you cannot do
   robust pseudobulk DE — fall back to a description, not a test. Pseudobulk
   design detail: `references/pseudobulk_workflow.md`.
3. **Covariates** — donor, sex, batch, cell-cycle score. For per-cell tests
   covariates go into `latent.vars` (MAST or LR test); for pseudobulk they
   go into the DESeq2 design formula (requires native DESeq2, not the Seurat
   wrapper — see `references/pseudobulk_workflow.md` §"Multi-factor designs").
4. **Adjusted p-value cutoff** — `p_val_adj < 0.05` is the Seurat default;
   reporting any "gene of interest" with `p_val_adj > 0.05` is reporting
   noise. Column semantics: `references/output_interpretation.md`.
5. **Effect-size cutoff** — `avg_log2FC` magnitude. `0.25` is permissive —
   use `≥ 0.5` for a stricter "biologically meaningful" set.
6. **Mode** — INTERACTIVE (default) when this is the primary analysis; BATCH
   when invoked as one of many DE comparisons in a larger plan. The
   orchestrator passes `args="batch"`; the agent declares mode in its
   `present_plan`.

Show the user these figures / tables as the analysis proceeds:
- `cluster_markers.csv` — full one-vs-rest markers, all clusters (Path A)
- `de_volcano.png` — volcano for the highlighted comparison (Path B / C)
- `de_topgenes_dotplot.png` — top hits per cluster / group (Path A / B)
- `pseudobulk_de.csv` — condition-effect table from the pseudobulk path (Path C)

For figure conventions (palette, theme, dpi), see `references/figure_style.md`.

---

## Path A — Cluster markers (one vs rest) with FindAllMarkers

The standard "what marks each cluster?" call. Wilcoxon, one cluster vs all
others, per cluster.

```r
# obj is your processed Seurat object with Idents() set to the cluster labels
# (set explicitly if it's not already the active identity):
#   Idents(obj) <- "seurat_clusters"      # or "wnn_clusters" / "seurat_annotations"

markers <- FindAllMarkers(
  obj,
  only.pos        = TRUE,    # one-sided: genes UP in each cluster vs rest
  min.pct         = 0.25,    # gene expressed in >=25% of either group
  logfc.threshold = 0.25,    # |avg_log2FC| >= 0.25 prefilter
  test.use        = "wilcox",
  verbose         = FALSE
)
# NOTE: Seurat 5.5.0's default is `logfc.threshold = 0.1, min.pct = 0.01`
# (loosened from earlier versions). We pin the historical PBMC3k-tutorial
# defaults (0.25 / 0.25) above for stability across Seurat patch releases.

# Top 5 markers per cluster by avg_log2FC
top5 <- markers %>%
        filter(p_val_adj < 0.05) %>%
        group_by(cluster) %>%
        slice_max(order_by = avg_log2FC, n = 5) %>%
        ungroup()

write.csv(markers, "cluster_markers.csv", row.names = FALSE)
cat(sprintf("Total significant markers (p_adj<0.05): %d\n",
            sum(markers$p_val_adj < 0.05)))
print(top5 %>% select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj),
      n = Inf)
```

> **`min.pct` and `logfc.threshold` are SPEED prefilters, not significance
> tests.** They reduce the gene set Seurat tests; the actual significance is
> `p_val_adj`. Set them low (e.g. `min.pct = 0.1`, `logfc.threshold = 0.1`)
> if you're looking for SUBTLE differences and don't mind a longer runtime.

### Dotplot of top markers

```r
genes_to_show <- unique(top5$gene)

p_dot <- DotPlot(obj, features = genes_to_show, cluster.idents = FALSE) +
  RotatedAxis() +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0, name = "avg expr") +
  ggtitle(sprintf("Top 5 markers per cluster (Wilcoxon, only.pos, n=%d genes)",
                  length(genes_to_show))) +
  theme_cowplot() +
  theme(plot.title  = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank())

suppressMessages(
  ggsave("de_topgenes_dotplot.png", p_dot,
         width = max(12, 0.18 * length(genes_to_show)),
         height = 6.5, dpi = 120, bg = "white")
)
```

**Report:** number of significant markers per cluster, the top 3 markers for
each cluster (by `avg_log2FC`), and any cluster with <5 markers (weakly
distinguished — may merge with a neighbor). For column semantics and
biological-coherence checks, see `references/output_interpretation.md`.

---

## Path B — Two-group FindMarkers (with optional covariates via MAST)

When the user names a specific comparison ("cluster 3 vs cluster 5",
"CD8 T-cells vs CD4 T-cells", "this cluster vs everything except cluster 7").

```r
# Pick the comparison — pass cluster IDs (or annotation labels) as ident.1 / ident.2.
# ident.2 = NULL means "compare ident.1 to ALL OTHER cells".
de <- FindMarkers(
  obj,
  ident.1         = "3",          # or "CD8 T cell"
  ident.2         = "5",          # or NULL for vs rest
  test.use        = "wilcox",
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = FALSE
)
de$gene <- rownames(de)
head(de[order(de$p_val_adj), ], 20)
```

### When to switch to MAST (per-cell, with covariates)

If the comparison is confounded by a per-cell covariate — donor effect within
a single-sample object, cell-cycle phase, percent.mt — adjust with MAST. MAST
is a zero-inflated two-part regression (Finak et al. 2015) that takes
`latent.vars` to condition out per-cell effects.

```r
# Install once: BiocManager::install("MAST")
de_mast <- FindMarkers(
  obj,
  ident.1     = "3",
  ident.2     = "5",
  test.use    = "MAST",
  latent.vars = c("percent.mt", "S.Score", "G2M.Score"),    # the per-cell confounds
  min.pct     = 0.25,
  logfc.threshold = 0.25,
  verbose     = FALSE
)
de_mast$gene <- rownames(de_mast)
```

`latent.vars` MUST be column names that exist in `obj@meta.data` AND be
numeric (Seurat's MAST wrapper coerces them into the `zlm` design as-is —
categorical covariates need one-hot encoding first). For the full `test.use`
catalogue (wilcox, wilcox_limma, MAST, LR, ROC, t, bimod, negbinom, poisson)
with one-line trade-offs per test and source citations, see
`references/per_cell_tests.md`.

> **MAST is still per-cell.** It treats cells as units of analysis. For a
> CONDITION effect across samples, go to Path C (pseudobulk) — MAST does NOT
> account for sample-level biological variance, just cell-level covariates.
> See `references/anti_patterns.md` §"MAST with donor as latent.vars" for the
> common confusion.

### Volcano plot

```r
de$nlp <- -log10(pmax(de$p_val_adj, .Machine$double.xmin))
de$sig <- with(de, ifelse(p_val_adj < 0.05 & abs(avg_log2FC) >= 0.5,
                          "sig", "ns"))

p_volc <- ggplot(de, aes(x = avg_log2FC, y = nlp, colour = avg_log2FC)) +
  geom_point(alpha = 0.5, size = 0.9) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0, name = "log2FC") +
  geom_hline(yintercept = -log10(0.05), colour = "red",
             linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = c(-0.5, 0.5), colour = "red",
             linetype = "dashed", linewidth = 0.5) +
  labs(title = "DE volcano (FindMarkers, Wilcoxon)",
       x = "avg_log2FC", y = "-log10(p_adj)") +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank())

ggsave("de_volcano.png", p_volc,
       width = 7, height = 5.5, dpi = 120, bg = "white")
```

**Report:** total tested genes, significant up / down, the top 10 by
`avg_log2FC` (filtered by `p_val_adj < 0.05`), and any covariates passed via
`latent.vars` if MAST was used.

---

## Path C — Condition effect across samples (PSEUDOBULK, the right way)

This is the bioinformatically correct path for "stim vs ctrl", "responders vs
non-responders", "disease vs healthy". The procedure:

1. **Aggregate** per-cell counts up to sample × cell-type pseudobulks with
   `AggregateExpression`.
2. **Test** with DESeq2 (or edgeR) on the sample-level pseudobulk matrix —
   each pseudobulk column IS a biological replicate now.

This is the Seurat v5 vignette pattern (satijalab.org/seurat/articles/de_vignette).
Methodological grounding: Squair et al. 2021 (Nat Commun) shows pseudobulk
methods match bulk-DE behavior on scRNA data; per-cell DE tools without
sample-level aggregation give inflated false-positive rates.

### Prerequisites

`obj@meta.data` must carry, per cell:
- a sample / donor identifier (e.g. `donor_id`),
- a condition label (e.g. `stim` with levels `"STIM"` / `"CTRL"`),
- a cell-type / cluster label (e.g. `seurat_annotations`).

Each donor × cell-type combination becomes ONE pseudobulk sample. You need
at least **3 donors per condition** for robust DE — fewer is descriptive at
best. If your fixture lacks a real `donor_id` (e.g. `ifnb` has only one
sample per condition), the Seurat v5 DE vignette mocks it for demonstration:
`set.seed(42); obj$donor_id <- sample(c('d1','d2','d3'), ncol(obj), replace=TRUE)`.
That exercises the code path but is NOT a biological claim — note it
explicitly in your report.

### Aggregate to pseudobulk

```r
# AggregateExpression is the canonical Seurat v5 pseudobulk function.
# return.seurat = TRUE wraps the aggregated counts back into a Seurat object
# whose columns are the (sample x celltype) pseudobulks (sum-of-counts by
# default; one column per group.by combination).
pseudo <- AggregateExpression(
  obj,
  assays      = "RNA",
  return.seurat = TRUE,
  group.by    = c("stim", "donor_id", "seurat_annotations")
)
# Joined cell-state ident, e.g. "CD14 Mono_STIM" / "CD14 Mono_CTRL":
pseudo$celltype.stim <- paste(pseudo$seurat_annotations, pseudo$stim, sep = "_")
Idents(pseudo) <- "celltype.stim"

cat(sprintf("Pseudobulk samples: %d (donors x celltypes x conditions)\n",
            ncol(pseudo)))
print(table(pseudo$stim, pseudo$seurat_annotations))
```

> **Verify replication.** Print `table(donor_id, stim)` — you need at least
> 3 donors per condition. If any cell-type has <3 donors in either condition
> after aggregation, EXCLUDE that cell-type from the test (the table call
> above shows you which to drop).

### Test with DESeq2 on pseudobulk

```r
# Install once: BiocManager::install("DESeq2")
# FindMarkers on the pseudobulk object routes through DESeq2 because each column
# is a biological replicate (donor x condition). NOTE the Seurat wrapper builds
# a `~group` design only — for multi-factor / paired / interaction designs
# (e.g. ~donor + condition) drop to native DESeq2 (see
# references/pseudobulk_workflow.md §"Multi-factor designs").
de_pb <- FindMarkers(
  object   = pseudo,
  ident.1  = "CD14 Mono_STIM",            # cell-type _ condition (one side)
  ident.2  = "CD14 Mono_CTRL",            # cell-type _ condition (other side)
  test.use = "DESeq2"
)
de_pb$gene <- rownames(de_pb)

write.csv(de_pb, "pseudobulk_de.csv", row.names = FALSE)
cat(sprintf("Pseudobulk DESeq2: %d significant (p_adj<0.05), %d up, %d down\n",
            sum(de_pb$p_val_adj < 0.05, na.rm = TRUE),
            sum(de_pb$p_val_adj < 0.05 & de_pb$avg_log2FC >  0, na.rm = TRUE),
            sum(de_pb$p_val_adj < 0.05 & de_pb$avg_log2FC <  0, na.rm = TRUE)))
```

### Loop over cell-types (one DE per cell-type)

For a thorough condition analysis, run pseudobulk DESeq2 per cell-type and
collate. Skip cell-types with <3 donors per condition.

```r
ct_list <- unique(pseudo$seurat_annotations)
de_all  <- list()
for (ct in ct_list) {
  ident_a <- paste0(ct, "_STIM");  ident_b <- paste0(ct, "_CTRL")
  cells_a <- WhichCells(pseudo, idents = ident_a)
  cells_b <- WhichCells(pseudo, idents = ident_b)
  if (length(cells_a) < 3 || length(cells_b) < 3) next  # need ≥3 donors / side
  d <- tryCatch(
    FindMarkers(pseudo, ident.1 = ident_a, ident.2 = ident_b, test.use = "DESeq2"),
    error = function(e) NULL
  )
  if (is.null(d)) next
  d$gene <- rownames(d); d$celltype <- ct
  de_all[[ct]] <- d
}
de_all <- do.call(rbind, de_all)
write.csv(de_all, "pseudobulk_de.csv", row.names = FALSE)
cat(sprintf("Pseudobulk DESeq2 (all cell-types): %d total sig (p_adj<0.05)\n",
            sum(de_all$p_val_adj < 0.05, na.rm = TRUE)))
```

### When the `~group` design isn't enough

The Seurat `FindMarkers(test.use = "DESeq2")` wrapper builds the design
`~group` (two-level factor; ident.1 vs ident.2). For a richer design
(`~donor + condition`, paired designs, interaction terms), drop to native
DESeq2 / edgeR — pull the matrix off the pseudobulk object via
`GetAssayData(pseudo, assay="RNA", layer="counts")` and build the
`DESeqDataSet` / `DGEList` yourself. The pattern is in
`references/pseudobulk_workflow.md` §"Multi-factor designs".

### Report

For the pseudobulk path, report:
- the design: # donors per condition, # cell-types tested, total pseudobulk columns
- which cell-types had insufficient replication and were excluded
- per cell-type: # significant genes (`p_val_adj < 0.05`), top up / down
- whether the DE direction makes biological sense — stim vs ctrl on `ifnb`
  should put **IFN-stimulated genes** (ISG15, IFIT1/2/3, CXCL10, RSAD2, MX1,
  OAS1) UP in STIM. Biological-coherence sanity-checks are documented in
  `references/output_interpretation.md`.

---

## Anti-pattern quick reference

For the full Squair-2021-grounded catalogue with citations, see
`references/anti_patterns.md`. The three most common live-session mistakes:

| Tempting wrong thing | Why it's wrong | Right answer |
|---|---|---|
| `DESeqDataSetFromMatrix(per_cell_counts, design = ~condition)` | Cells aren't replicates; FPR inflated 10–100× (Squair 2021) | Path C: `AggregateExpression` → DESeq2 |
| `FindMarkers(obj, test.use = "DESeq2")` directly on per-cell object for condition effect | Same — DESeq2 needs replicates, not cells | Path C: pseudobulk THEN DESeq2 |
| `MAST(latent.vars = "donor")` as a condition test across donors | MAST treats donor as per-cell covariate, not replicate; still anti-conservative | Path C |

---

## Batch variant — use INSTEAD of the path-specific code when invoked with args="batch"

Branch on `$ARGUMENTS` at the top of the body. In batch mode:

- Skip per-comparison volcano / dotplot figures.
- Drop the "Report" narrative blocks; print ONE final summary line per call.
- Still save the canonical CSV(s) — `cluster_markers.csv` for Path A,
  `pseudobulk_de.csv` for Path C.

Path A (cluster markers) batch — one call per Seurat object:

```r
obj <- readRDS("/path/to/seurat_processed.rds")
markers <- FindAllMarkers(obj, only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25,
                          test.use = "wilcox", verbose = FALSE)
write.csv(markers, "cluster_markers.csv", row.names = FALSE)
cat(sprintf("batch ok | %d clusters | %d significant markers (p_adj<0.05)\n",
            length(unique(markers$cluster)),
            sum(markers$p_val_adj < 0.05)))
```

Path C (pseudobulk DE) batch — one call per (cell-type, comparison):

```r
obj    <- readRDS("/path/to/seurat_processed.rds")
pseudo <- AggregateExpression(obj, assays = "RNA", return.seurat = TRUE,
                              group.by = c("stim", "donor_id",
                                           "seurat_annotations"))
pseudo$celltype.stim <- paste(pseudo$seurat_annotations, pseudo$stim, sep = "_")
Idents(pseudo) <- "celltype.stim"

de_pb <- FindMarkers(pseudo,
                     ident.1 = "CD14 Mono_STIM",
                     ident.2 = "CD14 Mono_CTRL",
                     test.use = "DESeq2")
de_pb$gene <- rownames(de_pb)
write.csv(de_pb, "pseudobulk_de.csv", row.names = FALSE)
cat(sprintf("batch ok | pseudobulk | %d sig (p_adj<0.05)\n",
            sum(de_pb$p_val_adj < 0.05, na.rm = TRUE)))
```

---

## Final response checklist

Summarize:
- which PATH you ran (A: cluster markers / B: two-group / C: pseudobulk
  condition effect) — and WHY that path matched the question
- input object: cells, clusters, sample / condition columns present
- test used (`test.use = ...`), any `latent.vars`, prefilter cutoffs
- effect-size threshold and adjusted p-value cutoff applied
- per-comparison: # genes tested, # significant, top 5 by `avg_log2FC`
- biological sanity: do the top hits make sense for the labelled groups
  (e.g. interferon-response genes UP in stim, lineage markers UP in their
  own cluster) — see `references/output_interpretation.md`
- saved files: `cluster_markers.csv` / `pseudobulk_de.csv` / volcano /
  dotplot
- caveats:
  - per-cell DE tools (Wilcoxon, MAST) treat cells as units — p-values are
    anti-conservative for cross-sample comparisons; use pseudobulk for
    condition effects (Squair 2021)
  - `min.pct` / `logfc.threshold` are speed prefilters, not significance
  - covariate adjustment via `latent.vars` only on `MAST` / `LR` / `negbinom`
    / `poisson` — Wilcoxon does not accept covariates
  - DESeq2 is for **pseudobulk only**, never per-cell
  - the Seurat `FindMarkers(test.use="DESeq2")` wrapper is `~group` design
    only; richer designs need native DESeq2 / edgeR

## See also

- **`seurat-scrna-v2`** — the prerequisite QC + clustering + initial
  `FindAllMarkers` recipe. This recipe is the deeper dive on the
  method-choice step that recipe ends on.
- **`seurat-integration`** — multi-sample integration with `IntegrateLayers`;
  remember to `JoinLayers()` before any DE on the integrated object.
- **`seurat-rna-atac-integration`** — DE on the joint WNN clusters of a 10x
  Multiome sample. The DE itself uses this recipe; the clusters come from
  there.
- **`deseq2-r`** / **`bulk_rnaseq_de`** — true bulk RNA-seq DE (one column
  per sample, no aggregation needed).
- **`seurat-index`** — the framework menu for picking the right Seurat
  recipe for your task.
