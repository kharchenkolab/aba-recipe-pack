# WNN internals — FindMultiModalNeighbors argument reference + artifact layout

Deep detail behind Step 2 of `seurat-wnn-multimodal`. Loaded when the agent
needs to defend a `dims.list` choice, tune `k.nn` for a small dataset,
diagnose a missing `RNA.weight` column after a CLR/SCT normalization, or
explain what each of `weighted.nn` / `wknn` / `wsnn` is for.

## Function signature (Seurat 5.5.0)

`formals(FindMultiModalNeighbors)` returns:

```r
function (object,
          reduction.list,                          # required: list of reduction names
          dims.list,                               # required: list of dim ranges, same length
          k.nn                  = 20,              # per-modality k-nearest neighbors
          l2.norm               = TRUE,            # L2-normalize each reduction before NN
          knn.graph.name        = "wknn",          # output kNN graph name
          snn.graph.name        = "wsnn",          # output SNN graph name
          weighted.nn.name      = "weighted.nn",   # output Neighbor object name
          modality.weight.name  = NULL,            # see "The single-string gotcha" below
          knn.range             = 200,             # candidate neighbor range for prediction
          prune.SNN             = 1/15,            # SNN graph pruning threshold
          sd.scale              = 1,               # bandwidth for the kernel
          cross.contant.list    = NULL,            # (sic) advanced — leave default
          smooth                = FALSE,           # smooth weights across modalities
          return.intermediate   = FALSE,           # for debugging weight learning
          modality.weight       = NULL,            # advanced — pre-computed ModalityWeight
          verbose               = TRUE)
```

Three required + named: `object`, `reduction.list`, `dims.list`. Everything
else has a sensible default for typical CITE-seq / multiome data.

## The single-string `modality.weight.name` gotcha (READ THIS)

Source-level (Seurat 5.5.0 R/integration_methods.R):

```r
modality.weight.name <- modality.weight.name %||%
                        paste0(modality.weight@modality.assay, ".weight")
if (length(modality.weight.name) != length(reduction.list)) {
  warning("The number of provided modality.weight.name is not equal to the number of modalities. ", ...)
  modality.weight.name <- paste0(modality.assay, ".weight")
}
```

What this means in practice:

| What you pass | What lands on the object |
|---|---|
| `modality.weight.name = NULL` | Per-assay `<assay>.weight` from the underlying assays — fragile if `DefaultAssay` is `"SCT"` (you get `SCT.weight`, not `RNA.weight`) |
| `modality.weight.name = "RNA.weight"` (single string) | **WARNS, then falls back** to per-assay names. The vignette uses this form and works only because the assays happen to be `"RNA"` / `"ADT"` |
| `modality.weight.name = c("RNA.weight", "ADT.weight")` | **Correct** — explicit per-modality names land on the object regardless of `DefaultAssay` |

The Step 2 code in the SKILL.md passes the 2-element form for this
reason. The vignette's single-string form *appears* to work but is
silently equivalent to passing nothing in many setups. **Always pass
the 2-element vector** — it's the explicit, defensive form.

When the second modality is ATAC, name it `ATAC.weight`:

```r
modality.weight.name = c("RNA.weight", "ATAC.weight")
```

The recipe's `## Step 6 / Final response checklist` references
`obj$RNA.weight` by name — if the user later inspects this column and
finds it missing, the diagnosis is almost always a single-string
`modality.weight.name` argument combined with `DefaultAssay = "SCT"`
(which lands `SCT.weight` on the object instead).

## `reduction.list` and `dims.list` — what to pass

`reduction.list` and `dims.list` are **paired R lists** of equal length.
Each entry of `reduction.list` is the string name of a reduction already
present on the object; each entry of `dims.list` is the integer range
of dimensions to keep from that reduction.

| Use case | `reduction.list` | `dims.list` |
|---|---|---|
| CITE-seq (RNA + ADT) — bmcite, vignette default | `list("pca", "apca")` | `list(1:30, 1:18)` |
| 10x Multiome (RNA + ATAC) — Signac convention | `list("pca", "lsi")` | `list(1:30, 2:30)` |
| Spatial + RNA (rare) | `list("pca", "ipca")` | `list(1:30, 1:50)` |
| Three modalities (RNA + ADT + ATAC, experimental) | `list("pca", "apca", "lsi")` | `list(1:30, 1:18, 2:30)` |

**Why LSI starts at component 2.** Component 1 of LSI captures
sequencing depth, not biology — a Signac/scATAC convention. Dropping
it is required (`dims.list = list(1:30, 2:30)` for multiome), not
optional.

**Why ADT often stops at 18.** A 25-protein panel has at most 24
informative PCs (rank constraint); the SLM clusterer doesn't benefit
from including PC 19–24 because they encode panel-specific noise. The
vignette empirically settled on 18 for bmcite. For larger panels
(40–50 markers) raise the upper bound (`1:25` or `1:30`); for smaller
(<15 markers) lower it (`1:10`).

**The dims-tuning loop.** If you don't know the right cap:

```r
ElbowPlot(obj, ndims = ncol(Embeddings(obj, "apca")), reduction = "apca")
# Read the knee — typically where stdev drops below 1% of total.
# The first knee is the dim cap; below it is signal, above is noise.
```

The recipe's defaults (`1:30` RNA, `1:18` ADT) work for the bmcite-shape
panel out of the box; tune only if the elbow disagrees.

## `k.nn` and `knn.range` — typically leave default

`k.nn = 20` is the per-modality nearest-neighbor count used for the
weight-learning step. The default fits datasets of 5k–100k cells. For
very small datasets (<2k cells), drop to `k.nn = 10` so the
neighborhoods don't engulf the whole cluster.

`knn.range = 200` is the candidate pool size for the prediction step —
how far to look for matching cells when computing the modality weight.
Lower means stricter local fit; higher means more global smoothing.
**Don't change** unless a methods-paper reviewer asks for a
sensitivity analysis.

## Artifact placement — where each output lives

Three output structures end up in three different slots:

```r
obj@neighbors$weighted.nn   # Neighbor object — used by RunUMAP(nn.name=)
obj@graphs$wknn             # k-NN graph    — adjacency, sparse Matrix
obj@graphs$wsnn             # SNN graph     — used by FindClusters(graph.name=)
obj$RNA.weight              # per-cell weight metadata column (0..1)
obj$ADT.weight              # per-cell weight for the OTHER modality (1 - RNA.weight if 2 modalities)
```

Access patterns:

```r
# Verify all three before downstream calls — STEP 2 sanity-check.
stopifnot("weighted.nn"   %in% names(obj@neighbors))
stopifnot(all(c("wknn", "wsnn") %in% names(obj@graphs)))
stopifnot("RNA.weight"     %in% colnames(obj@meta.data))

# Read the Neighbor object directly if needed (rare).
nn <- obj[["weighted.nn"]]              # IDX matrix + DIST matrix, k columns each

# Inspect graph density.
sum(obj[["wknn"]]@x > 0) / length(obj[["wknn"]]@x)   # fraction of non-zero entries
```

## `l2.norm = TRUE` — when to flip

`l2.norm` scales each cell's reduction vector to unit length before the
neighbor search. This is the **right default for PCA-vs-PCA / PCA-vs-LSI**
where the two reductions have different scales. Keep `TRUE`.

The only time to set `FALSE` is when ALL reductions are already on a
comparable scale (e.g. both are L2-normalized SCT residuals, or both
are spectral embeddings on the same scale). In practice every recipe
uses the default.

## `prune.SNN = 1/15` — the SNN density tuner

`prune.SNN` is the Jaccard-coefficient threshold for keeping an edge
in the SNN graph. The default `1/15 ≈ 0.067` is right for most
datasets. Two failure modes:

- **Too dense graph** (everything is connected to everything) →
  `FindClusters` produces 2–3 mega-clusters. Bump `prune.SNN` to
  `1/10 = 0.1`.
- **Too sparse graph** (clusters fragment) → `FindClusters` produces
  >50 clusters. Drop to `1/20 = 0.05`.

In WNN specifically, `prune.SNN` interacts with the per-cell weight:
high-weight cells contribute thicker edges to one modality's
neighborhood, which can make the joint graph denser than either
single-modality graph. The default usually handles this fine.

## `sd.scale = 1` — the kernel bandwidth

`sd.scale` controls the bandwidth of the kernel that converts
neighborhood-prediction distances into a weight. Lower = sharper
weights (closer to 0/1 per cell); higher = softer (closer to 0.5 per
cell).

Symptom-driven adjustment:

- **All cells `RNA.weight ≈ 0.5`** — the weights aren't separating
  populations. Try `sd.scale = 0.5` for sharper weights.
- **Bimodal weights with clear modality dominance per cluster but the
  joint UMAP looks no different from the better single-modality UMAP**
  — weights are too sharp. Try `sd.scale = 2`.

Default works for CITE-seq + multiome out of the box.

## `smooth = FALSE` — leave alone unless you have a reason

`smooth = TRUE` averages each cell's weight with its neighbors'
weights — an extra step the vignette doesn't use. Keeps `FALSE`
unless you observe high-frequency speckle in the modality_weight
plot that doesn't track biology.

## What `weighted.nn` is — the structure WNN exports

After `FindMultiModalNeighbors`, the `weighted.nn` Neighbor object
carries, for each cell, the indices of its `k.nn = 20` nearest
neighbors in the **weighted** distance — where the weight is
applied per cell to combine the two reductions. This is what
`RunUMAP(nn.name = "weighted.nn", ...)` consumes to compute the joint
UMAP. You can also feed it back into other downstream tools:

```r
# Use the WNN graph for label-transfer-style operations (advanced):
# Find cells in cluster X with the closest WNN neighbor structure.
weighted_nn_idx <- obj[["weighted.nn"]]@nn.idx
```

This is a power-user pattern; the SKILL.md's normal path is just
`RunUMAP` + `FindClusters` consuming `weighted.nn` and `wsnn`
respectively.

## References

- Hao Y., et al. (2021). *Integrated analysis of multimodal single-cell
  data.* Cell 184, 3573–3587. doi:10.1016/j.cell.2021.04.048. The WNN
  paper; methods section is the canonical description of the
  weight-learning algorithm.
- Seurat WNN vignette
  <https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis>
  — bmcite worked example with the `dims.list = list(1:30, 1:18)`
  default.
- `?FindMultiModalNeighbors` (Seurat 5.5.0) — argument documentation;
  `formals(FindMultiModalNeighbors)` is the authoritative source for
  default values used here.
