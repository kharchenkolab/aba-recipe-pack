# Clustering choices — `FindClusters` algorithm, resolution, and `dims`

The decisions that shape `FindNeighbors` → `FindClusters` → `RunUMAP`:
which algorithm, what resolution, how many PCs (`dims`), and how to tell
when each is wrong.

Load this reference when:
- The cluster count looks wrong (too many, too few, fragmented populations).
- The user asks "why Louvain vs Leiden?"
- A UMAP looks weird (over-mixed populations, split clusters, artifactual
  bridges).
- You're considering switching `algorithm` or tuning `resolution`.

## The three steps and their default semantics

Source-verified against Seurat 5.5.0 (`args(Seurat:::FindNeighbors.Seurat)`,
`args(Seurat:::FindClusters.Seurat)`, `args(Seurat:::RunUMAP.Seurat)`):

```r
FindNeighbors(obj,
              reduction = "pca", dims = 1:10,        # NOTE the default
              k.param = 20, return.neighbor = FALSE,
              compute.SNN = TRUE,
              prune.SNN = 1/15,
              nn.method = "annoy", annoy.metric = "euclidean",
              ...)

FindClusters(obj,
             modularity.fxn = 1,
             resolution = 0.8,                       # NOTE the default
             algorithm = 1,                          # 1 = Louvain
             method = NULL,
             leiden_method = c("leidenbase", "igraph"),
             leiden_objective_function = c("modularity", "CPM"),
             n.start = 10, n.iter = 10, random.seed = 0,
             ...)

RunUMAP(obj,
        dims = NULL, reduction = "pca",
        umap.method = "uwot",                        # uwot, not umap-learn
        n.neighbors = 30L, n.components = 2L,
        metric = "cosine",                           # was "correlation" in v4
        min.dist = 0.3, spread = 1,
        seed.use = 42L,
        ...)
```

Two things change vs the pre-v5 defaults the PBMC3k tutorial documents:

- **`FindNeighbors` default `dims = 1:10`** is the Seurat default but
  the recipe overrides with `dims = 1:DIMS_CHOSEN` (typically 30) for
  larger samples. The vignette's `1:10` is correct for ~3k PBMCs; v5
  routine practice has converged on ~30 because most modern samples
  carry more cell-type heterogeneity than PBMC3k.
- **`RunUMAP` defaults to `umap.method = "uwot"` + `metric = "cosine"`** in
  v5. The pre-v5 default was `umap.method = "umap-learn"` (Python) + `metric
  = "correlation"`. UMAP layouts are NOT identical across major versions; for
  cross-session reproducibility on a pre-v5 analysis, pass
  `umap.method = "umap-learn", metric = "correlation"` explicitly.

## `dims` — how many PCs feed the neighbor graph

The PCs you pass to `FindNeighbors` and `RunUMAP` determine which axes of
biological variation the clustering sees. Too few PCs: structure compressed,
distinct populations bleed together. Too many PCs: noise PCs dominate
late-PC distance contributions, clusters fragment.

### Picking `DIMS_CHOSEN` from the elbow plot

The recipe's PCA elbow plot shows per-PC + cumulative variance. Three
heuristics, in order of robustness:

1. **The cumulative curve plateau.** If cumulative variance levels off
   sharply by PC 15, use `dims = 1:15`. If it climbs steadily through PC 50,
   use `dims = 1:30` (don't go to `1:50` — the late PCs are mostly noise).
2. **The per-PC step.** A visible "shoulder" — where the per-PC variance
   drops by an order of magnitude — separates signal from noise PCs.
3. **Sanity check against the heatmap.** `DimHeatmap(dims = 1:10, ...)`
   should show coherent gene loadings on early PCs. If PC 8's heatmap is
   random noise, the structure stops by PC 7.

PBMC3k vignette uses `dims = 1:10`; v5 routine practice uses `dims = 1:30`
on most non-PBMC samples (because broader cell-type heterogeneity demands
more axes). The recipe defaults to 30 and overrides off the elbow.

### Excluding a technical PC

If PC 1 is dominated by `MT-*` genes (recipe Step 5 prints loadings to
catch this), the options are:

1. **Tighten QC** — go back to Step 3 with a lower `percent.mt` ceiling.
2. **Regress in `ScaleData`** — `ScaleData(obj, vars.to.regress = "percent.mt")`
   removes the linear contribution. Slow but clean.
3. **Skip the PC** — `dims = 2:30` excludes PC 1. Quick but blunt.

The recipe defaults to (1); (2) and (3) are surfaced when (1) isn't an
option (e.g. PC1 is cell cycle and the biology is in cycling cells).

## `FindClusters(algorithm = …)` — four options

Source-verified table:

| Value | Method | Speed | Dependencies | When to use |
|---|---|---|---|---|
| `1` (default) | Original Louvain (Blondel 2008) | Fast | None (built-in C++) | Default. Modular communities on SNN graph |
| `2` | Louvain with multilevel refinement | Slightly slower than 1 | None | Marginal quality improvement; rarely worth the slowdown for first-pass |
| `3` | SLM (Smart Local Moving, Waltman 2013) | Comparable to 1 | None | Some benchmarks suggest tighter clusters than Louvain; experimental in most workflows |
| `4` | Leiden (Traag 2019) | Slower than 1 | Python `leidenalg` via `reticulate` (`reticulate::py_install("leidenalg")`) OR R `leidenbase` package | Higher-quality communities by guaranteeing well-connectedness. Slow to wire up but the gold standard for publication-quality clustering |

`leiden_method = c("leidenbase", "igraph")` controls which Leiden
implementation is used when `algorithm = 4` — `leidenbase` (the R port)
when present, falling back to `igraph` otherwise. Pass
`leiden_method = "leidenbase"` explicitly to skip the Python detour.

For the recipe's default Louvain (`algorithm = 1`), no extra setup is
needed. If the user explicitly requests Leiden, install `leidenbase` or
`leidenalg` first.

`modularity.fxn = 1` (default) is standard modularity; `modularity.fxn = 2`
is an alternative formulation rarely needed for default scRNA workflows.

## `resolution` — controls cluster COUNT, not identity

The single most-tuned parameter. Higher resolution → more, smaller
clusters; lower → fewer, larger. The recipe defaults to `0.5` (moderate);
the vignette uses `0.5`; SCT-based workflows often use `0.8` because SCT
recovers finer structure than LogNormalize.

Rough mapping for PBMC3k:
- `resolution = 0.1`: 4–5 clusters (broad lineages: T, B, myeloid, NK)
- `resolution = 0.3`: 7–8 clusters (lineage + subset distinction)
- `resolution = 0.5`: 8–10 clusters (standard PBMC granularity)
- `resolution = 0.8`: 11–14 clusters (sub-lineage detail)
- `resolution = 1.5`: 18–22 clusters (often over-clustered)

### Crucial gotcha — resolution does NOT preserve cluster identity

Cluster 3 at `resolution = 0.5` is NOT "cluster 3" at `resolution = 0.8`.
The community-finding algorithm restarts from scratch at each call; the
numeric labels are arbitrary outputs of the algorithm. Implications:

- **Pick one resolution per analysis.** Don't toggle and assume the same
  cluster persists.
- **Marker calling, cluster annotation, downstream sub-clustering all
  depend on a fixed resolution.** Changing it invalidates the marker table.
- **`RenameIdents` is keyed by cluster label.** If you label cluster 3 as
  "CD4 mem T" then re-run `FindClusters`, the label vanishes.

### Sanity-check via cluster size distribution

After `FindClusters`, the recipe prints `table(Idents(obj))`. Heuristics:

- **A cluster < 1% of total cells (e.g. <30 cells in PBMC3k):** often a
  doublet bridge or a genuine rare population. Inspect the markers — if
  they mix two parent populations, it's a doublet bridge; if coherent,
  it's a rare cell type to flag.
- **One mega-cluster > 50% of cells:** resolution may be too low. Bump to
  0.8 / 1.0.
- **Most clusters < 50 cells:** resolution too high. Drop to 0.3.

## Multi-resolution exploration

Seurat stores every `FindClusters` result as a metadata column named
`<assay>_snn_res.<R>`:

```r
obj <- FindClusters(obj, resolution = c(0.3, 0.5, 0.8))
# Creates obj$RNA_snn_res.0.3, RNA_snn_res.0.5, RNA_snn_res.0.8
```

Useful for the `clustree` package — visualizes how cluster membership shifts
across resolutions, helping pick a stable point. Optional; the recipe stays
at one resolution by default.

## UMAP — visualization, NOT clustering

The clusters come from `FindClusters` (operating on the SNN graph). The UMAP
is a 2D projection FOR VISUALIZATION. Common confusions:

- **Read cluster identity from `obj$seurat_clusters` (the active idents),
  NOT from drawing polygons on the UMAP.** A cluster might be split into
  two visual islands on UMAP and still be one cluster (and vice versa).
- **UMAP is non-deterministic without a seed.** `seed.use = 42L` is the
  default, but if you `library(uwot)` and pull `uwot::umap` directly, set
  the seed. Different runs of the SAME object can produce mirror-image or
  rotated layouts — clusters preserve identity but visual positions shift.
- **UMAP distances are NOT meaningful.** Two clusters far apart on UMAP
  aren't necessarily more dissimilar than two adjacent ones. Use Euclidean
  distance in PCA space for that, not UMAP coordinates.

## `min.dist` and `spread` — UMAP layout tuning

When clusters look too tightly packed or smeared:

| Parameter | Default | Effect of increasing |
|---|---|---|
| `min.dist` | `0.3` | Looser packing, more empty space between clusters |
| `spread` | `1.0` | Cluster diameter grows; useful when clusters overlap visually |
| `n.neighbors` | `30L` | Larger n_neighbors smooths local structure (more global) |

Don't tune these to make a story look better — they're aesthetic, not
algorithmic. Tune `dims` and `resolution` first; touch UMAP params only
when the visualization is genuinely hard to read.

## When to skip clustering and jump to reference mapping

For tissues with a high-quality reference atlas (PBMC → Azimuth's
pbmcref; brain → Allen atlas; gut → HCA gut atlas), `SingleR` or
`Azimuth::RunAzimuth` annotates cells directly without clustering. Faster
and less subjective than cluster-then-annotate.

The recipe doesn't take this branch by default because (a) most users want
to see the clustering structure of their own data first, and (b) reference
mapping has its own recipe (`seurat-reference-mapping`,
`seurat-multimodal-reference-mapping`).
