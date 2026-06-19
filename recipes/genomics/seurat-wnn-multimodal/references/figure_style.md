# Figure style — Seurat collection (shared across recipes)

This is the **shared figure-style cheatsheet** for every recipe in the
Seurat collection (`seurat-scrna-v2`, `seurat-integration`,
`seurat-sctransform`, `seurat-cite-seq`, `seurat-wnn-multimodal`,
`seurat-de-testing`, `seurat-reference-mapping`,
`seurat-multimodal-reference-mapping`, `seurat-rna-atac-integration`).
Duplicated verbatim in each recipe's `references/figure_style.md` so the
agent reading one recipe gets the full styling without needing to load a
second file.

Loaded when the agent is composing a plot — palette choice, theme
ordering quirks, `ggsave` conventions, the alpha-poke pattern. The
boilerplate is short; copy-paste discipline.

## Boilerplate at the top of every plotting block

```r
suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(cowplot)
})
```

If the recipe uses `patchwork` (most do), add `library(patchwork)` to
the block too.

## `theme_cowplot()` on every ggsave figure

Always append `theme_cowplot()` (NOT `theme_bw`, NOT `theme_minimal`).
The cowplot signature is what makes the collection's figures look
consistent — no panel border, clean axis lines, faint grey grid where
the recipe explicitly adds one.

```r
p <- DimPlot(obj, reduction = "umap", label = TRUE) +
  theme_cowplot() +
  theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank())
```

## `theme_cowplot()` and `NoLegend()` ordering — the cowplot 1.2.0 gotcha

cowplot 1.2.0 changed how `theme_cowplot()` handles `legend.position`.
**Put `NoLegend()` AFTER `theme_cowplot()`** or the legend stays
visible:

```r
# WRONG — legend reappears under cowplot 1.2.0+
p + NoLegend() + theme_cowplot()

# RIGHT — NoLegend wins because it's last
p + theme_cowplot() + NoLegend()
```

When broadcasting across a patchwork:

```r
# WRONG — '+' adds to the FIRST sub-panel only
p_patch + theme_cowplot() + NoLegend()

# RIGHT — '&' broadcasts to ALL sub-panels
p_patch & theme_cowplot() & NoLegend()
```

## `ggsave` conventions

Bare relative filename. White background. 120 dpi for most figures;
180 dpi for heatmaps. Explicit width × height in inches:

```r
ggsave("umap_clusters.png", plot = p, bg = "white", dpi = 120,
       width = 7, height = 6.5, units = "in")
```

NEVER use `file.path(WORK_DIR, …)` or any absolute path inside the
recipe — the kernel's cwd IS the thread's scratch dir, and the
artifact-harvester scans it. Absolute paths leak outside the harvester.

## Palettes — the named set

### Diverging (zero-centered) — log-fold-change, z-scores, residuals

```r
scale_colour_gradient2(low  = "#2166ac",      # ColorBrewer RdBu blue
                       mid  = "grey90",
                       high = "#b2182b",      # ColorBrewer RdBu red
                       midpoint = 0)
# For fills (heatmaps): scale_fill_gradient2(...) with the same colors.
```

Broadcast across a patchwork with `&`:

```r
p_patch & scale_colour_gradient2(low = "#2166ac", mid = "grey90",
                                  high = "#b2182b", midpoint = 0)
```

### Sequential — expression, counts, weights (always non-negative)

```r
scale_colour_gradient(low = "grey85", high = "#b2182b")
# For fills: scale_fill_gradient(low = "grey85", high = "#b2182b").
```

### Categorical fill (violins, ridges) — small N

```r
scale_fill_brewer(palette = "Set2", guide = "none")
# Use ONLY when N <= 8 — Set2 has 8 colors and recycles silently above that.
# For N > 8 categories (typical WNN clustering at resolution=2 gives 15-25):
scale_fill_discrete(guide = "none")    # ggplot's hue_pal — N distinct hues
```

### Kept/filtered (QC scatter) — categorical, two levels

```r
scale_colour_manual(
  values = c(kept = "black", filtered = "red"),
  guide  = guide_legend(override.aes = list(alpha = 1, size = 1.5))
)
```

## Plot-type conventions

### `VlnPlot`

```r
geom_violin(width = 0.85, colour = "black", linewidth = 0.4,
            scale = "width", trim = FALSE, alpha = 0.85)
# Always `pt.size = 0` to suppress overlaid jitter — the violins
# already encode density. Re-enable jitter only when the cluster
# sample sizes are small (<200 cells).
```

### Threshold lines (QC plots)

```r
geom_hline(yintercept = <thresh>, colour = "red",
           linetype = "dashed", linewidth = 0.5)
```

### `DimPlot`

```r
DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE,
        pt.size = 0.4) + NoLegend()
```

Cluster labels on the embedding (not in a side legend); `repel = TRUE`
to avoid label overlap.

### `DotPlot`

```r
DotPlot(obj, features = genes_to_show, cluster.idents = FALSE) +
  RotatedAxis()
# Sizing: width = max(12, 0.18 * length(genes_to_show)); height = 6.5
```

### `FeaturePlot`

```r
FeaturePlot(obj, features = feats,
            reduction = "umap", order = TRUE, pt.size = 0.3,
            ncol = 3)
# Sizing: height = max(4, 4 * n_rows); width = 13
# `order = TRUE` puts high-expression cells on top — non-default but
# essential for low-prevalence markers.
```

### `DimHeatmap`

Append `& theme(legend.position = "none")` to drop the per-panel legend
stack — heatmap legends crowd the figure.

### PCA elbow

Per-PC bars (solid `#1f77b4` blue) + cumulative `grey40` dashed line on
the **same** axis (not side-by-side panels). Y-axis as percent of total
HVG-matrix variance — `sdev^2 / sum(sdev^2) * 100`.

## Alpha-poke — DimPlot, FeaturePlot, VariableFeaturePlot

Seurat does NOT expose `alpha` directly on `DimPlot`, `FeaturePlot`, or
`VariableFeaturePlot`. Walk `p$layers` for `GeomPoint` and set
`aes_params$alpha` post-hoc:

```r
# 0.6 for UMAP/DimPlot/FeaturePlot — preserves cluster structure
# 0.35 for HVG MA-style scatter — denser overplotting
for (i in seq_along(p$layers)) {
  if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
    p$layers[[i]]$aes_params$alpha <- 0.6
  }
}
```

For a patchwork, nest the loop in `for (panel in seq_along(p_patch))`:

```r
for (panel in seq_along(p_patch)) {
  pl <- p_patch[[panel]]
  if (!is.null(pl$layers)) {
    for (k in seq_along(pl$layers)) {
      if (inherits(pl$layers[[k]]$geom, "GeomPoint")) {
        pl$layers[[k]]$aes_params$alpha <- 0.6
      }
    }
    p_patch[[panel]] <- pl
  }
}
```

## Patchwork: broadcast theme/scale with `&`, NOT `+`

Critical and easy to miss. `+` adds to the first sub-panel only; `&`
broadcasts to every sub-panel:

```r
# WRONG
(b1 | b2 | b3) + theme_cowplot() + NoLegend()

# RIGHT
(b1 | b2 | b3) & theme_cowplot() & NoLegend()
```

## Modality-specific notes

For modality-specific plots (CITE-seq biaxial, spatial overlays, ATAC
coverage), follow the same discipline: theme_cowplot, ggsave with
explicit bg/dpi, named palettes from the list above. If a Seurat
function returns a patchwork object, broadcast theme/scale with `&`.

## Pitfalls quick reference

- `NoLegend() + theme_cowplot()` — legend reappears under cowplot
  1.2.0+. Swap order.
- `+ theme_cowplot()` on a patchwork — applies only to the first
  sub-panel. Use `& theme_cowplot()`.
- `scale_fill_brewer(palette = "Set2")` on >8 levels — silently
  recycles colors. Use `scale_fill_discrete(...)` for N > 8.
- `FeaturePlot` without `order = TRUE` — low-prevalence markers buried
  under zeros.
- DimPlot/FeaturePlot without alpha-poke — overplotting hides
  cluster structure. Walk `p$layers` and set `aes_params$alpha`.
- Absolute path inside `ggsave` — figure lands outside the harvester's
  scratch dir. Bare relative filename only.
- Transparent background (`bg = NULL` default on some themes) —
  confuses the figure-display layer. Always `bg = "white"`.
