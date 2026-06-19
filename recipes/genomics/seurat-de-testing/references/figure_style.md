# Figure style

The styling conventions shared across the Seurat recipes — palette, theme,
dpi, alpha-poke patterns. Duplicated verbatim into each Seurat-collection
recipe so the agent reading one recipe gets the full styling without
needing to load a second file. If you change this section, propagate
to siblings.

Load this reference when:
- You're authoring or modifying a figure-producing step in this recipe.
- A reviewer flags a palette / theme / dpi regression.
- The validation loop's "Per-figure checks" (palette, background, theme)
  fails and you need the spec it's checking against.

## ggsave / png — common arguments

```r
# ggplot — single-panel, white background, 120 dpi. Bare filename: the
# kernel's cwd IS the thread's scratch dir, and the artifact-harvester
# scans cwd for outputs to register. No file.path(WORK_DIR, ...).
ggsave("de_volcano.png", p_volc,
       width = 7, height = 5.5, units = "in", dpi = 120, bg = "white")

# Base/grid (heatmaps etc.) — wrap with png()/dev.off()
png("marker_heatmap.png",
    width = 13.5, height = 8.1, units = "in", res = 120, bg = "white")
# ... plot calls ...
dev.off()
```

Conventions:

- **`bg = "white"`** — never leave transparent backgrounds; they confuse
  the figure-display layer (the chat UI's PNG renderer doesn't composite
  transparent PNGs against a known background).
- **`dpi = 120`** — readable in the chat UI without bloating files.
- **`.png` always** — unless the user explicitly requests otherwise.
- **Bare relative filenames.** The kernel's cwd is the per-run scratch dir
  the harvester scans; absolute paths or `file.path(WORK_DIR, …)`
  constructs may land outside it (and `WORK_DIR` is itself a banned
  vignette idiom). Never write to `/tmp` or any other path the harvester
  doesn't see.
- **Filename matches `produces` frontmatter exactly.** The recipe's
  `produces` lists `de_volcano.png`, `de_topgenes_dotplot.png`,
  `cluster_markers.csv`, `pseudobulk_de.csv` — every figure-saving call
  uses one of those names.

## Theme

`theme_cowplot()` on every ggplot figure. It gives clean axis lines, no
panel border, and a flat black-on-white look that composites well into
the chat UI. Add the recipe's standard overlay tweaks on top:

```r
+ theme_cowplot()
+ theme(plot.title       = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank())
```

For DotPlots specifically, rotate the x-axis labels:

```r
+ RotatedAxis()
+ theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        axis.text.y = element_text(size = 10),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank())
```

## Palette — diverging blue→grey→red, centered on zero

For any plot whose colour axis is **signed** (log fold change, z-scores,
correlations), use the diverging Brewer-RdBu-derived palette:

```r
+ scale_colour_gradient2(low   = "#2166ac",
                         mid   = "grey90",
                         high  = "#b2182b",
                         midpoint = 0,
                         name  = "log2FC")    # or "avg expr", "z-score", etc.
```

- `#2166ac` — Brewer-RdBu blue, the "down" end.
- `grey90` — neutral mid, slightly off-white so it doesn't disappear against
  the white background.
- `#b2182b` — Brewer-RdBu red, the "up" end.
- `midpoint = 0` — center the palette on biological "no change".

The same scale applies for fill-based plots: swap `scale_colour_gradient2`
for `scale_fill_gradient2` with the same args.

**Do NOT use viridis on signed axes.** viridis (yellow→purple) is a
sequential palette designed for unsigned magnitude. Using it on
`avg_log2FC` confuses the reader because positive and negative shifts are
not visually symmetric. This is the most common figure-style regression
the validation-loop sweep catches.

For unsigned magnitudes (counts, fractions, p-value strengths), use
sequential — typically `scale_*_gradient(low = "grey95", high = "#b2182b")`
or `scale_*_distiller(palette = "Reds")`. Not viridis there either, but
that's a softer rule.

## Threshold lines on volcano / scatter plots

```r
+ geom_hline(yintercept = -log10(0.05), colour = "red",
             linetype = "dashed", linewidth = 0.5)
+ geom_vline(xintercept = c(-0.5, 0.5), colour = "red",
             linetype = "dashed", linewidth = 0.5)
```

- **Solid `colour = "red"`** for the cutoff lines so they read as
  threshold markers, not data.
- **`linetype = "dashed"`** distinguishes them from any solid trend or
  fit lines.
- **`linewidth = 0.5`** — visible but not dominant.

## Point styling on volcano / scatter plots

```r
+ geom_point(alpha = 0.5, size = 0.9)
```

- `alpha = 0.5` — semi-transparent so overplotting reveals density.
  Opaque dots (`alpha = 1`) hide the data shape; the validation-loop
  "alpha-poke" check flags this.
- `size = 0.9` — small enough that ~10k points fit without saturating;
  bump to 1.2 for smaller plots or to 0.6 if the plot is very dense.

For DimPlot (not used in this recipe, but for consistency with siblings):

```r
DimPlot(obj, ..., pt.size = 0.4, label = TRUE, repel = TRUE) +
  NoLegend() +
  theme_cowplot()
# Walk p$layers and set alpha = 0.6 on the points (Seurat doesn't expose alpha):
# p$layers[[1]]$aes_params$alpha <- 0.6
```

## Dimensions and aspect ratios

Recipe-specific defaults that work across most fixtures:

| Plot type | Default size | Notes |
|---|---|---|
| Volcano | `width = 7, height = 5.5` | Slightly wider than tall — reads as horizontal scatter |
| DotPlot (top markers) | `width = max(12, 0.18 * n_genes), height = 6.5` | Width scales with gene count |
| Marker heatmap | `width = 13.5, height = 8.1` | Wider; rows = genes, cols = cells |
| DimPlot (UMAP) | `width = 7, height = 6.5` | Roughly square (with cluster label space) |

All at `dpi = 120, bg = "white"`.

## Sources

- Brewer palettes — `RColorBrewer::brewer.pal(11, "RdBu")` for the
  blue/grey/red triple (Brewer 2003); the `#2166ac` / `#b2182b` hex
  values are entries 2 and 10 of that palette.
- cowplot — Wilke 2020, ggplot2-themes wiki (github.com/wilkelab/cowplot)
- Seurat v5 plotting vignette — satijalab.org/seurat/articles/visualization_vignette
- The validation loop's per-figure checklist —
  `aba-skill-authoring/references/validation_loop.md` §"Per-figure checks"
