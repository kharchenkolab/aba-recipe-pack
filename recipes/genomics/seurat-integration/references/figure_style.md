# Seurat collection — shared figure style

When to load this: producing any plot from a Seurat-collection recipe;
fixing a figure whose theme/palette/legend looks off; bringing a new
plot into the collection's style.

This reference is **duplicated verbatim** across every recipe in the
seurat-* collection. The collection-level decision (see
`aba-skill-authoring/references/body_structure.md` §9) is to keep
shared conventions in each recipe rather than centralize them — zero
coupling between siblings, at the cost of light drift if someone
forgets to copy the latest version. When you update this file, copy
the same change into the other seurat-* sibling recipes.

## Backbone — `theme_cowplot()`

Every plot in the collection uses `theme_cowplot()` as its base
theme. It gives a clean white background, no panel border, faint
axis lines, no major-grid by default (we add a subtle one back
manually).

```r
ggplot_call + theme_cowplot()
```

Apply theme modifications AFTER `theme_cowplot()`. Order matters
because cowplot ≥1.2.0 sets `legend.position` and other slots that a
later `theme()` will override.

## Critical ordering — `NoLegend()` AFTER `theme_cowplot()`

For DimPlots with `label = TRUE`, suppress the side legend so labels
on the embedding are the only legend. Place `NoLegend()` AFTER
`theme_cowplot()`:

```r
p <- DimPlot(obj, ..., label = TRUE, repel = TRUE, pt.size = 0.4) +
  ggtitle(...) +
  theme_cowplot() +          # 1. cowplot base
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.4)) +
  NoLegend() +               # 2. AFTER theme_cowplot — cowplot 1.2.0 overrides legend.position
  coord_fixed()
```

Reversed order (`+ NoLegend() + theme_cowplot()`) silently restores
the side legend on some cowplot versions. The pilot caught this
across ≥5 sibling recipes.

## Alpha-poke for DimPlot / FeaturePlot

`DimPlot` and `FeaturePlot` do NOT expose an `alpha` argument
directly. To get semi-transparent points (which read better at any
cell count), walk the plot object's layers and set `alpha` on
each `GeomPoint` layer's `aes_params`:

```r
for (k in seq_along(p$layers)) {
  if (inherits(p$layers[[k]]$geom, "GeomPoint")) {
    p$layers[[k]]$aes_params$alpha <- 0.6
  }
}
```

Apply AFTER all other plot composition (after `+ NoLegend() +
coord_fixed()`). The walk is idempotent — re-applying does no harm.

## ggsave defaults

```r
ggsave("<filename>.png", p,
       width = 7, height = 6, units = "in", dpi = 120, bg = "white")
```

Conventions:

- **`bg = "white"`** — always. Transparent backgrounds confuse the
  figure-display layer.
- **`dpi = 120`** — readable in chat UI, files stay reasonably small.
- **`.png`** extension always, unless the user explicitly requests
  another.
- **Bare relative filename** — the kernel's cwd is the per-run scratch
  dir the harvester scans; NEVER write to `/tmp` or use
  `file.path(WORK_DIR, ...)`.
- **`width = 7, height = 6` for single DimPlots, `coord_fixed()`** —
  the square panel ratio matches `coord_fixed()`. Side-by-side
  patchwork: `width = 14, height = 6`. ElbowPlot / line charts:
  `width = 7, height = 4.8`.

## Palette conventions

- **Categorical (cell types, clusters, samples)** — Seurat's default
  hue palette is fine for ≤12 categories. For >12, switch to
  `scale_color_brewer(palette = "Set1")` or a hand-defined vector.
- **Diverging (log fold change, score signed)** —
  `scale_color_gradient2(low = "#2166ac", mid = "grey90", high =
  "#b2182b")`. NEVER use viridis for diverging scales (yellow-purple
  doesn't read as positive vs negative).
- **Sequential (expression intensity, density)** — near-white →
  saturated red works:
  `scale_color_gradient(low = "grey90", high = "#b2182b")`.

## ElbowPlot / line charts — dual-axis convention

For elbow plots showing per-PC + cumulative on one panel:

- per-PC: solid line, colour `#1f77b4` (blue)
- cumulative: dashed line, colour `grey40`
- chosen-dims marker: red dashed vline + small annotation

```r
ggplot(df_long, aes(PC, value, colour = kind, linetype = kind)) +
  geom_vline(xintercept = DIMS_CHOSEN, colour = "red",
             linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = DIMS_CHOSEN - 0.7, y = max(df_long$value) * 0.95,
           label = sprintf("dims = 1:%d (chosen)", DIMS_CHOSEN),
           colour = "red", hjust = 1, size = 3.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
  scale_colour_manual(values = c("per-PC" = "#1f77b4",
                                 "cumulative" = "grey40")) +
  scale_linetype_manual(values = c("per-PC" = "solid",
                                   "cumulative" = "dashed")) +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        legend.position = c(0.98, 0.55),
        legend.justification = c(1, 0.5),
        legend.background = element_rect(fill = alpha("white", 0.5), colour = NA),
        legend.title = element_blank())
```

## Base graphics / pheatmap / ComplexHeatmap — wrap with `png()`

For non-ggplot output (heatmaps, base R plots):

```r
png("<filename>.png",
    width = 13.5, height = 8.1, units = "in", res = 120, bg = "white")
# draw call here
dev.off()
```

Note `res = 120` (base graphics naming) vs `dpi = 120` (ggsave).

## patchwork side-by-side

```r
library(patchwork)
p_side <- (p_left | p_right) +
  plot_annotation(title = "<title>",
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))
ggsave("<filename>.png", p_side,
       width = 14, height = 6, units = "in", dpi = 120, bg = "white")
```

When both panels use `coord_fixed()`, the side-by-side layout
preserves equal aspect ratios.

## Title style

```r
ggtitle(sprintf("%s — %d cells, %d clusters", method_label, ncol(obj),
                length(unique(obj$cluster)))) +
  ...
  theme(plot.title = element_text(size = 12, face = "bold"), ...)
```

Long sprintf titles can clip at the right edge when combined with
`coord_fixed()` + side legend at `width = 7in`. Either widen to
`width = 8.5` or shorten the title.

## What NOT to do

- **Transparent backgrounds.** Always set `bg = "white"` on both
  ggsave and png.
- **`file.path(WORK_DIR, ...)`** — that's a vignette idiom that leaks
  variable names into tool calls. Use bare filenames.
- **Viridis for diverging.** Yellow-to-purple doesn't read as
  positive-vs-negative.
- **`+ NoLegend() + theme_cowplot()`** — wrong order; cowplot 1.2.0
  restores the side legend.
- **`pt.size = 0` or `pt.size = 1.5`+** — both unread. Stick with
  `pt.size = 0.4` for single-panel UMAPs, `pt.size = 0.6` for
  side-by-side.
