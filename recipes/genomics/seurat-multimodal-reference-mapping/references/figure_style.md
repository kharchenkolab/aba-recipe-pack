# Figure style — Seurat collection conventions

Plot styling shared by every recipe in the Seurat collection. Apply
verbatim. Load this when the agent needs to author a new plot, or
when a figure regression (palette swap, missing alpha-poke, missing
white background) needs to be diagnosed.

This file is duplicated VERBATIM across all sibling Seurat recipes
(per the meta-skill's `body_structure.md` §9 "Current default:
duplicate the shared section in each recipe"). Do not rewrite it per
recipe — if a convention drifts, update the foundation
(`seurat-scrna-v2`) first, then re-copy to the siblings.

## Boilerplate at the top of every R block that plots

```r
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(cowplot)
})
```

`library(Seurat)` does NOT attach `ggplot2`/`dplyr`/`cowplot` — load
them by name. `patchwork` is auto-loaded by Seurat for `&`-composed
plots, but attach it explicitly if you build patchworks manually.

## ggsave defaults

```r
ggsave(<filename>, plot = <p>,
       bg = "white", dpi = 120,
       width = <w>, height = <h>)
```

Conventions:

- `bg = "white"` — never leave transparent backgrounds; they confuse
  the figure-display layer.
- `dpi = 120` — readable in the chat UI without bloating files.
- `dpi = 180` for heatmaps (`DimHeatmap`, marker heatmaps) — they
  need finer rendering.
- File extension `.png` always, unless the user requests otherwise.
- **Bare relative filenames.** The kernel's cwd is the per-run
  scratch dir the artifact-harvester scans; absolute paths or
  `file.path(WORK_DIR, …)` may land outside it.
- Filename must match the `produces` frontmatter entry exactly.

For base/grid plots (e.g. `DimHeatmap` legacy mode):

```r
png(<filename>, width = <w>, height = <h>, units = "in",
    res = 120, bg = "white")
<plot call>
dev.off()
```

## theme_cowplot()

Apply `theme_cowplot()` on every ggsave figure (white background, no
panel border, clean axis lines). For patchwork compositions, broadcast
with `&`:

```r
p <- p1 | p2 | p3
p & theme_cowplot()
```

Optional grid additions on top of `theme_cowplot()`:

```r
theme(
  plot.title = element_text(size = 12, face = "bold"),
  panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
  panel.grid.minor = element_blank()
)
```

## Palettes

**Diverging (zero-centered).** Use for log-fold-change, modality
weights, z-scored expression, anything signed:

```r
scale_*_gradient2(low      = "#2166ac",
                  mid      = "grey90",
                  high     = "#b2182b",
                  midpoint = 0)
```

Substitute `midpoint = 0.5` and `limits = c(0, 1)` for modality
weights (RNA.weight / ATAC.weight) — the diverging axis runs 0 →
0.5 → 1.

**Sequential.** Use for unsigned magnitudes (gene expression,
nUMI, percent.mt):

```r
scale_colour_gradient(low = "grey85", high = "#b2182b")
```

For FeaturePlot, prefer the `cols = c("lightgrey", "<terminal>")`
argument; use `cols = c("lightgrey", "#1b7837")` (green) to mark
*imputed* / *predicted* values (predicted_ADT, predicted gene
activity) so they're visually distinct from measured RNA (red).

**Categorical fill.** For 2–8 categories use Set2; for >8, switch
to scale_fill_discrete (which uses the default ggplot palette and
doesn't recycle silently):

```r
scale_fill_brewer(palette = "Set2", guide = "none")    # ≤ 8 categories
scale_fill_discrete(guide = "none")                    # > 8 categories
```

**Kept / filtered (QC scatter).** For QC plots that color cells by
"will be kept" vs "will be removed":

```r
scale_colour_manual(values = c(kept = "black", filtered = "red"),
                    guide = guide_legend(override.aes = list(alpha = 1, size = 1.5)))
```

## Common plot conventions

**VlnPlot:**

```r
geom_violin(width = 0.85, colour = "black", linewidth = 0.4,
            scale = "width", trim = FALSE, alpha = 0.85)
```

**Threshold lines (QC):**

```r
geom_hline(yintercept = <thresh>, colour = "red",
           linetype = "dashed", linewidth = 0.5)
```

**DimPlot:**

```r
p <- DimPlot(obj, reduction = "umap",
             label = TRUE, repel = TRUE, pt.size = 0.4) +
  theme_cowplot() + coord_fixed() +
  NoLegend()         # AFTER theme_cowplot — cowplot 1.2.0 overrides legend.position
```

`coord_fixed()` keeps modality UMAPs visually comparable when stacked
side by side (e.g. `umap.rna` / `umap.atac` / `umap.wnn`).

**DotPlot:**

```r
DotPlot(obj, features = genes_to_show, cluster.idents = FALSE) +
  RotatedAxis()
# width = max(12, 0.18 * length(genes_to_show)); height = 6.5
```

**FeaturePlot:**

```r
FeaturePlot(obj, features = features, order = TRUE,
            pt.size = 0.3, ncol = 3)
# n_rows <- ceiling(length(features) / 3)
# width = 13; height = max(4, 4 * n_rows)
```

**PCA elbow:**

per-PC solid blue (`#1f77b4`) bars + cumulative `grey40` dashed line
on the SAME axis (not side-by-side panels), y-axis as percent of
total HVG-matrix variance.

**DimHeatmap:**

append `& theme(legend.position = "none")` to drop the per-panel
legend stack.

## Alpha-poke (DimPlot / FeaturePlot / VariableFeaturePlot)

Seurat does not expose `alpha` directly. After building the plot,
walk `p$layers` for `GeomPoint` and set `aes_params$alpha`:

```r
for (i in seq_along(p$layers)) {
  if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
    p$layers[[i]]$aes_params$alpha <- 0.6
  }
}
```

Values:

- `0.6` for UMAP / DimPlot / FeaturePlot.
- `0.35` for HVG MA-style scatter (`VariableFeaturePlot`).

For patchwork (FeaturePlot ncol > 1):

```r
for (i in seq_along(p)) {
  pl <- p[[i]]
  if (!is.null(pl$layers)) {
    for (k in seq_along(pl$layers)) {
      if (inherits(pl$layers[[k]]$geom, "GeomPoint")) {
        pl$layers[[k]]$aes_params$alpha <- 0.6
      }
    }
    p[[i]] <- pl
  }
}
```

## Patchwork composition (Seurat patchworks)

Seurat returns patchworks for multi-panel calls (FeaturePlot ncol>1,
some CoveragePlot panels). Broadcast theme/scale with **`&`**, not
`+`:

```r
p & theme_cowplot()                                  # right
p + theme_cowplot()                                  # wrong (applies to last panel only)

p & scale_colour_gradient2(low = "#2166ac",
                           mid = "grey90",
                           high = "#b2182b")
```

## Modality-specific exceptions

- **Signac CoveragePlot.** Don't overlay `theme_cowplot()` — Signac's
  default styling is genome-browser tuned (gene tracks, peak
  annotations, expression dots). Adding cowplot breaks the layout.
- **CITE-seq biaxial plots.** Use `theme_cowplot()` + the standard
  palette; manually log10-transform protein UMI counts before plotting
  (so the scale is comparable to RNA).
- **Spatial overlays.** Use `theme_cowplot()` + sequential palette
  (red); set `image.alpha = 0` if the underlying tissue image makes
  the colormap unreadable.

## Sizing cheatsheet

| Plot type | width × height (in) | dpi |
|---|---|---|
| Single UMAP / DimPlot | 7 × 6.5 | 120 |
| DimPlot with long titles | 8 × 6.5 | 120 |
| FeaturePlot 3-col | 13 × `4 × n_rows` | 120 |
| Violin (1 metric, many groups) | 8 × 4 | 120 |
| QC violin (4 metrics side by side) | 10 × 4.5 | 120 |
| DotPlot (marker panel) | `max(12, 0.18 × n_genes)` × 6.5 | 120 |
| Heatmap (DimHeatmap) | 12 × 8 | 180 |
| Heatmap (markers, ComplexHeatmap) | 13.5 × 8.1 | 120 |
| Coverage plot | 9 × 6.5 | 120 |
| Histogram (faceted) | 11 × 7 | 120 |
| Modality weights (FeaturePlot) | 7.5 × 6.5 | 120 |
