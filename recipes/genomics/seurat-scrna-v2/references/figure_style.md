# Figure style — Seurat collection conventions

The opinionated figure conventions every Seurat recipe in this collection
follows. Duplicated verbatim into each Seurat recipe so the agent reading
one recipe gets the full styling without loading a second file (per
`aba-skill-authoring/references/body_structure.md` §9).

Load this reference when:
- A figure regressed from the canonical style (legend in wrong slot, wrong
  palette, alpha=1 dots on UMAP) and you need to trace the convention.
- Authoring a new figure in a recipe and you want to copy the conventions
  rather than invent.

## File-level conventions

| Property | Value | Why |
|---|---|---|
| Resolution | `dpi = 120` for most plots, `dpi = 180` for high-detail heatmaps (DimHeatmap, complex marker heatmaps) | Readable in chat UI without bloating PNG sizes >1 MB |
| Background | `bg = "white"` ALWAYS | Transparent backgrounds confuse the figure-display layer and pixel-sampling sanity checks |
| Format | `.png` (not pdf, not svg, not jpeg) | Consistent across the collection; the artifact harvester is keyed off `.png` |
| Filename | Bare relative filename, no path prefix | The kernel's cwd IS the run's scratch dir — `ggsave("foo.png", ...)` writes where the harvester scans. `file.path(WORK_DIR, "foo.png")` is banned (see `translating_vignettes.md` Idiom 1) |
| Width / height | `width = 7, height = 6.5` for square-ish single-panel; `width = 12, height = 5.3` for a 2-panel row; `width = 13, height = 4` per row of FeaturePlot | Empirical for chat-UI readability |

## Theme — `theme_cowplot()` is the canonical theme

```r
library(cowplot)
# After every ggplot construction:
+ theme_cowplot()
```

`theme_cowplot()` strips the grey panel background, drops gridlines, adds
black axis lines, sets axis text to 12pt. It's the publication-ready Seurat
look. NEVER mix `theme_minimal()` or `theme_classic()` into the same recipe.

### Theme ordering trap — `NoLegend()` must come AFTER `theme_cowplot()`

```r
p + theme_cowplot() + NoLegend() + coord_fixed()  # CORRECT
p + NoLegend() + theme_cowplot() + coord_fixed()  # WRONG — legend reappears
```

`theme_cowplot()` (cowplot 1.2.0) RE-SETS `legend.position` in its `theme()`
call, overriding any earlier `NoLegend()` (which is itself a
`theme(legend.position = "none")` wrapper). The fix is to put `NoLegend()`
AFTER `theme_cowplot()`. Same applies to `theme(legend.position = ...)`
called manually.

This bug was caught live on the pilot recipe (Round 1 patch in
`test-seurat-scrna-v2/TEST_LOG.md`). The pattern silently produced a stray
legend on every UMAP DimPlot until the order was reversed.

### Title sizing trap — `coord_fixed()` + long titles

The DimPlot UMAP commonly carries an informative title like
`"Louvain res=0.5 . dims=1:30 . n=2698 cells . 8 clusters"`. With
`theme_cowplot()`'s default 14pt title and `coord_fixed()` constraining the
draw width, the title overflows at `width = 7` in. Two fixes (in order):

1. **Shrink the title:** `+ theme(plot.title = element_text(size = 12, face = "bold"))`.
   Apply this between `theme_cowplot()` and `NoLegend()`.
2. **Widen the figure:** bump `width = 7.5` or `width = 8`. Use when (1) is
   not enough (e.g. an even longer SCT title with 13 clusters).

The recipes apply both for any DimPlot that has more than `nClusters > 8`
or `coord_fixed()`.

## Palette — diverging blue / grey / red for signed values

For ANY plot showing a signed quantity (log-fold change, residuals, scaled
expression, signed correlation):

```r
+ scale_fill_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                      midpoint = 0)
# or for color aesthetic:
+ scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                        midpoint = 0)
```

Hex codes are RColorBrewer's `RdBu` endpoints (semi-canonical for biology
viz). NEVER use viridis (`scale_colour_viridis_c()`) for signed values — it
runs yellow → purple, which has no semantic mapping to "negative" vs
"positive".

For SEQUENTIAL quantities (gene expression on a UMAP, anything ≥ 0):

```r
+ scale_colour_gradient(low = "grey85", high = "#b2182b")
```

Grey to red. Low values are de-emphasized; high values pop. The default
Seurat FeaturePlot palette (blue → red) is acceptable but less consistent
with the rest of the collection — overlay the explicit scale.

## Per-plot details

### DimPlot (UMAP coloured by cluster)

```r
p <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE,
             pt.size = 0.4) +
  ggtitle("<algorithm>, res=<R>, dims=1:<D>, n=<N> cells, <K> clusters") +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold")) +
  NoLegend() +              # AFTER theme_cowplot — see ordering trap above
  coord_fixed()

# Alpha-poke — DimPlot has no `alpha` argument. Walk p$layers, find the
# GeomPoint, set aes_params$alpha. Public ggplot2 API; safe across versions.
for (k in seq_along(p$layers)) {
  if (inherits(p$layers[[k]]$geom, "GeomPoint")) {
    p$layers[[k]]$aes_params$alpha <- 0.6
  }
}

ggsave("umap_clusters.png", plot = p,
       bg = "white", dpi = 120, width = 7, height = 6.5)
```

`label = TRUE, repel = TRUE` puts cluster numbers ON the embedding (not in
a side legend); `repel = TRUE` uses ggrepel to keep them from overlapping.
`coord_fixed()` is required — UMAP coordinates aren't physical units, but
NOT fixing the aspect ratio makes the embedding warp under window resizes.

### FeaturePlot (UMAP coloured by gene expression)

```r
p <- FeaturePlot(obj, features = canonical, order = TRUE,
                 pt.size = 0.3, ncol = 3) &
  scale_colour_gradient(low = "grey85", high = "#b2182b") &
  theme_cowplot() &
  theme(legend.position = "right", legend.key.size = unit(0.4, "cm"))

# Alpha-poke each sub-panel
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

Note the `&` operator (NOT `+`) — FeaturePlot returns a `patchwork` object
when given multiple genes, and `&` applies a theme/scale to every panel.
Using `+` only modifies the top-level patchwork (which has no aesthetics).

`order = TRUE` plots high-expression cells on top of low ones — without it,
positive cells get visually buried under the grey background.

### DotPlot (cluster markers)

```r
p <- DotPlot(obj, features = genes_to_show, cluster.idents = FALSE) +
  RotatedAxis() +
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0, name = "avg expr") +
  ggtitle(...) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10))

ggsave("markers_dotplot.png", plot = p, bg = "white", dpi = 120,
       width = max(12, 0.18 * length(genes_to_show)), height = 6.5)
```

`cluster.idents = FALSE` keeps clusters in numeric order. Set to TRUE only
to group clusters hierarchically (helpful for cell-type validation, not for
diagnostic dot plots).

`width = max(12, 0.18 * length(genes_to_show))` scales width by gene count
so the rotated x-axis labels don't crowd.

### DimHeatmap (PC inspection)

```r
p <- DimHeatmap(obj, dims = 1:10, cells = 500,
                balanced = TRUE, fast = FALSE, combine = TRUE,
                ncol = 5) &
  scale_fill_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                       midpoint = 0) &
  theme(legend.position = "none")

ggsave("pca_heatmap.png", plot = p,
       bg = "white", dpi = 180, width = 20, height = 9)
```

`fast = FALSE` is the load-bearing argument — returns real ggplots so the
diverging `scale_fill_*` applies. `fast = TRUE` (the default) uses `image()`
under the hood and ignores ggplot scales, giving you Seurat's purple/yellow
default no matter what.

`combine = TRUE` + `ncol = 5` arranges into a 2×5 grid for 10 PCs. The
`& theme(legend.position = "none")` drops 10 identical per-panel legends.

`dpi = 180` because heatmaps need fine detail; `width = 20, height = 9` for
the 2×5 grid.

### Pre-filter QC violins (kept/filtered overlay)

The QC violin pattern uses a `pivot_longer` reshape so all four metrics
(`nFeature`, `nCount`, `percent.mt`, `percent.ribo`) live on one figure via
`facet_wrap(~ metric, ncol = 1, scales = "free_x")`. Threshold lines drawn
via `geom_hline` (red dashed, `linewidth = 0.5`). Cells coloured by a
`qc_kept` factor (black = kept, red = filtered). Fills via Set2 brewer
palette (categorical), one color per metric.

Full code in the SKILL.md Step 3 block; the convention is:
- `coord_flip()` so metric labels read horizontally on the left strip.
- `strip.placement = "outside"` + `strip.background = element_blank()`
  for a clean left-side metric label.
- `geom_jitter(width = 0.32, height = 0, size = 0.20, alpha = 0.10)` for
  the cell density underlay (10% opacity so the violin envelope stays
  legible).
- `legend.position = "bottom"` for the kept/filtered legend (the metric
  fill is `guide = "none"`).

## Pre-filter scatter (`patchwork` composition)

```r
p <- (s1 | s2) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "QC scatters, pre-filter",
                  theme = theme(plot.title = element_text(size = 13, face = "bold"))) &
  theme(legend.position = "bottom")
```

`|` puts two plots side by side. `plot_layout(guides = "collect")` merges
the two identical kept/filtered legends into one. `&` applies the
`legend.position = "bottom"` to BOTH sub-panels (necessary; the
`plot_annotation` `theme = ...` only affects the wrapper).

## PCA elbow — per-PC + cumulative on shared axis

The elbow plot fuses two lines on one axis:
- **Per-PC variance** (`prop`): solid blue (`#1f77b4`, matplotlib default).
- **Cumulative variance** (`cum`): dashed grey (`grey40`).

Both normalized to the TOTAL HVG-matrix variance (not the sum of the 50
PCs). This makes the cumulative curve plateau at the true fraction of
variance the PCs capture (~15–30% for a single 10x sample — never 100%).

`geom_vline(xintercept = DIMS_CHOSEN, colour = "red", linetype = "dashed")`
+ a small `annotate("text", ...)` callout shows the chosen dim cutoff.

The legend goes inline (`legend.position = c(0.98, 0.55)`) so the plot uses
its full width for data.

## Inline-display helper (used in test notebooks, not in production recipe)

For the Phase 2 validation notebooks, every `ggsave` is followed by:

```r
embed_png <- function(path, w = NULL, h = NULL) {
  bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  IRdisplay::display_png(data = bytes, width = w, height = h)
}
ggsave("foo.png", plot = p, ...)
embed_png("foo.png")
```

This is purely a notebook artifact for the rendered HTML — production
recipes don't include it (the agent's chat UI displays the saved PNGs
directly via the artifact harvester). See `validation_loop.md` §"Inline
display — required, not optional" for the rationale.

## What NOT to do

- ❌ `theme_minimal()`, `theme_classic()`, `theme_bw()` — not the
  collection's convention.
- ❌ `scale_colour_viridis_c()` for signed quantities — viridis has no
  midpoint semantics; use `scale_colour_gradient2`.
- ❌ `pt.size = 1` on UMAP — opaque dense dots obscure structure. Default
  `0.4` + alpha-poke 0.6.
- ❌ `NoLegend() + theme_cowplot()` ordering — the legend reappears. Put
  `NoLegend()` LAST.
- ❌ `file.path(WORK_DIR, ...)` for output paths — the harvester only
  scans the cwd. Bare filename always.
- ❌ Absolute paths in `ggsave` — write outside the harvester scan and
  the artifact won't be registered.
- ❌ `dpi = 300` on every figure — bloats PNG sizes >2 MB and the chat UI
  doesn't need that resolution. Reserve for actual print figures.
