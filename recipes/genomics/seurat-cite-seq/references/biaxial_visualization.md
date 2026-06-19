# Biaxial visualization — FeatureScatter, FeaturePlot, and the `adt_`/`rna_` keys

Deep detail behind Step 5 (FeatureScatter biaxials) and Step 6 (FeaturePlot
overlays + ADT-vs-RNA pairs) of `seurat-cite-seq`. Loaded when the agent
needs to switch FeatureScatter axes between CLR and raw counts, adapt
the panel to a different organism, or read protein > RNA / RNA > protein
discordance off the cross-modality plot.

## The `adt_<protein>` / `rna_<gene>` feature-key convention

Seurat exposes every assay through a key that overrides `DefaultAssay()`:

```r
Key(obj[["RNA"]])    # "rna_"
Key(obj[["ADT"]])    # "adt_"
```

Any Seurat plotting function that takes a `features =` argument accepts
those prefixed names and routes the query to the right assay regardless
of `DefaultAssay()`. This is the **explicit assay** pattern the SKILL.md
prefers over `DefaultAssay(obj) <- "ADT"; FeatureScatter(...)`:

```r
# Without keys — relies on DefaultAssay being ADT at this call site
DefaultAssay(obj) <- "ADT"
FeatureScatter(obj, feature1 = "CD4", feature2 = "CD8")
# Fragile — a later step that sets DefaultAssay back to RNA breaks this.

# With keys — assay is explicit in the feature name
FeatureScatter(obj, feature1 = "adt_CD4", feature2 = "adt_CD8")
# Robust — independent of DefaultAssay; safe to mix in a patchwork.
```

The agent should **always use the prefixed form** in plots. Even when
DefaultAssay is already what you want, the prefix makes the assay
contract visible to a human reader of the recipe.

### Cross-modality plots — protein vs RNA on the same call

The prefixed-key convention is what makes Step 6's ADT-vs-RNA pair plot
work as a single `FeaturePlot` call:

```r
FeaturePlot(obj, features = c("adt_CD3", "rna_CD3E"), reduction = "umap")
```

Without keys you'd have to switch DefaultAssay between the two
features — which `FeaturePlot` doesn't support in one call. The key
syntax is non-negotiable for cross-modality work.

## `FeatureScatter` — CLR-normalized vs raw counts

Default behavior reads from the data layer (CLR-normalized, after
Step 3):

```r
FeatureScatter(obj, feature1 = "adt_CD4", feature2 = "adt_CD8")
```

To switch to raw counts (useful when comparing to a flow-cytometry FCS
file with the same panel):

```r
# Seurat v4: slot = "counts"
# Seurat v5: layer = "counts" (slot still accepted as alias)
FeatureScatter(obj, feature1 = "adt_CD4", feature2 = "adt_CD8",
               slot = "counts")
```

The raw axes are unbounded and log-scaled visually; CLR axes are
~0-centered and roughly Gaussian per protein. For diagnostic biaxials
the CLR view is preferred — quadrants are stable across cells. For
"does this match the flow plot we ran last month" the raw view is
what you want.

## Panel adaptation — three canonical biaxials per study type

The SKILL.md's Step 5 shows three PBMC biaxials. Adapt to other study
types:

### PBMC / blood (the SKILL default)

```r
b1 <- FeatureScatter(obj, feature1 = "adt_CD4",  feature2 = "adt_CD8")   # T-cell subsets
b2 <- FeatureScatter(obj, feature1 = "adt_CD3",  feature2 = "adt_CD19")  # T vs B
b3 <- FeatureScatter(obj, feature1 = "adt_CD14", feature2 = "adt_CD16")  # monocyte subsets
```

### Bone marrow (bmcite-style)

Lowercase `CD8a` — bmcite's panel name; verify with `rownames(obj[["ADT"]])`.

```r
b1 <- FeatureScatter(obj, feature1 = "adt_CD34", feature2 = "adt_CD38")  # HSPC subsets
b2 <- FeatureScatter(obj, feature1 = "adt_CD4",  feature2 = "adt_CD8a")  # T-cell subsets
b3 <- FeatureScatter(obj, feature1 = "adt_CD14", feature2 = "adt_CD16")  # monocytes
```

### Tonsil / secondary lymphoid

```r
b1 <- FeatureScatter(obj, feature1 = "adt_CD3",   feature2 = "adt_CD19")
b2 <- FeatureScatter(obj, feature1 = "adt_CD27",  feature2 = "adt_IgD")  # memory vs naive B
b3 <- FeatureScatter(obj, feature1 = "adt_CXCR5", feature2 = "adt_PD-1") # Tfh
```

In every case: verify the marker exists with
`feat %in% rownames(obj[["ADT"]])` before plotting — a typo or
panel-mismatch gives an obscure ggplot error mid-patchwork.

## `FeatureScatter` returns a Pearson r in the title

By default the panel title is the Pearson correlation between the two
features (e.g. `"-0.55"` for CD4 vs CD8a — anticorrelated by lineage).
If you want a custom title:

```r
b1 <- FeatureScatter(obj, "adt_CD4", "adt_CD8") +
  ggtitle("CD4 vs CD8 (CLR-normalized)")    # overrides the r-in-title
```

Reading the default r-titles is the cheap interpretation lever:
strong negative r on a lineage pair (CD4 vs CD8, CD3 vs CD19) is the
sanity check that the panel and clustering are doing the right thing.

## `FeaturePlot` overlay on the RNA UMAP — sequential red gradient

For ADT overlays on the RNA UMAP (Step 6), the recipe sets a
sequential `grey85 → #b2182b` gradient via the patchwork `&` operator:

```r
p_adt_fp <- FeaturePlot(obj,
                        features = paste0("adt_", panel_show),
                        reduction = "umap",
                        order = TRUE,
                        pt.size = 0.3,
                        ncol = 3) &
  scale_colour_gradient(low = "grey85", high = "#b2182b")
```

Three details to keep right:

- **`&` not `+`.** `FeaturePlot` with multiple features returns a
  patchwork object; `+` adds layers to ONE panel, `&` broadcasts to
  ALL. The recipe uses `&` correctly.
- **`order = TRUE`.** Cells with high expression are drawn on TOP so
  the positive signal is visible. Default `FALSE` paints in observation
  order; positives are buried under negatives for any low-prevalence
  marker.
- **The "Scale for colour is already present" warning** is cosmetic —
  `FeaturePlot` installs its own default scale per sub-panel and the
  `&` then overrides it. Suppress with `suppressMessages(ggsave(...))`
  if it's distracting; do NOT switch off `&` to silence it.

## ADT-vs-RNA pair plot — reading protein/RNA discordance

The Step 6 pair plot lays out each lineage marker as protein on the
left, RNA on the right (2-column ncol). Three patterns to read:

| Pattern | Meaning | Action |
|---|---|---|
| **Both gradients track the same cluster** | Antibody + transcript agree; cluster is real | Good — annotate cluster with this lineage |
| **Protein high in cluster, RNA flat** | Translation lag — surface-stable protein remains while the transcript turns off (memory T cells, naive B cells, plasma cells) | Biology, not error. Mark the cluster as "post-transcriptional" / "memory" |
| **RNA high in cluster, protein flat** | Recently-induced transcript before the surface protein has accumulated (activated cells), OR cleaved/internalized protein | Possible activation state; corroborate with another activation marker |
| **Both flat across all clusters** | Either the antibody failed (re-check panel QC) or the gene isn't expressed (wrong organism, wrong symbol) | Audit symbol mapping; see `references/adt_qc_and_panel.md` "Common symbol mismatches" |
| **Anticorrelated within one cluster** | Two cell states within the cluster — protein-high and RNA-high subsets | Consider sub-clustering, or step up to WNN (`seurat-wnn-multimodal`) |

## Alpha-poke is mandatory for biaxials too

`FeatureScatter` and `FeaturePlot` don't expose `alpha` directly.
Without alpha, overplotting at dense quadrant corners hides density
structure. Walk `p$layers` for `GeomPoint` and set `aes_params$alpha`
(0.6 for both — same as DimPlot):

```r
for (i in seq_along(p$layers)) {
  if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
    p$layers[[i]]$aes_params$alpha <- 0.6
  }
}
```

For a patchwork the loop sits inside an outer `for (i in
seq_along(patchwork))` so it broadcasts to each sub-panel —
the SKILL.md's Step 5 code already does this.

## Quick reference — feature-key plot calls

```r
# Single-panel biaxial
FeatureScatter(obj, "adt_CD4", "adt_CD8")

# Patchwork biaxial — three panels with shared theme + alpha-poke
(b1 | b2 | b3) & theme_cowplot() & NoLegend()

# ADT overlay on RNA UMAP
FeaturePlot(obj, features = "adt_CD3", reduction = "umap", order = TRUE)

# Cross-modality pair on one call
FeaturePlot(obj, features = c("adt_CD3", "rna_CD3E"), reduction = "umap",
            order = TRUE, ncol = 2)

# Raw counts on a biaxial (flow-style)
FeatureScatter(obj, "adt_CD4", "adt_CD8", slot = "counts")    # v5 accepts slot=
```

## References

- Seurat multimodal vignette
  <https://satijalab.org/seurat/articles/multimodal_vignette> — source
  of the `slot = "counts"` raw-axes pattern and the `adt_`/`rna_` key
  convention.
- `?FeatureScatter` and `?FeaturePlot` (Seurat 5.5.0) — `order=TRUE`
  default semantics and the `slot`/`layer` argument behavior.
