# Mapping internals — what FindTransferAnchors / TransferData / MapQuery actually do

When to load this: you need to tune `k.anchor` / `k.filter` / `npcs`;
the anchor count or median score looks off; you want to understand
`pcaproject` vs `cca` vs `lsiproject`; you need to debug a low-confidence
mapping; you need to interpret `prediction.score.max` vs
`mapping.score`.

The SKILL.md body shows the canonical pcaproject + LogNormalize +
1:30 dims pattern with `MapQuery` wrapping the three sub-calls. This
reference adds: what each sub-call does, the inner args' semantics
(from `?FindTransferAnchors`, `?TransferData`, `?MapQuery`, Seurat
5.5.0), score interpretation, and patches for low-anchor / low-score
outcomes.

## The three-stage anchor → transfer → project pipeline

```
FindTransferAnchors(reference, query, …)
  └─ outputs: AnchorSet (anchor cell pairs + per-anchor scores)

TransferData(anchorset, refdata, …)
  └─ outputs: data.frame of per-query-cell predictions (predicted.id +
              prediction.score.<label> columns + prediction.score.max)

IntegrateEmbeddings(anchorset, reference, query, …)
  └─ outputs: query cells embedded in the reference's PCA space (ref.pca)

ProjectUMAP(query[["ref.pca"]], reference[["umap"]]@misc$model)
  └─ outputs: query cells projected onto the reference's UMAP (ref.umap)
```

`MapQuery` is the convenience wrapper that runs all four; the SKILL.md
Step 4 uses it. Step 3 uses only `TransferData` (Path A).

## `FindTransferAnchors` — arg semantics

Per `?FindTransferAnchors` (Seurat 5.5.0):

```
FindTransferAnchors(reference, query,
  normalization.method = "LogNormalize",   # or "SCT"
  recompute.residuals = TRUE,              # SCT-only
  reference.assay = NULL,
  reference.neighbors = NULL,
  query.assay = NULL,
  reduction = "pcaproject",                # see "Reduction choice" below
  reference.reduction = NULL,              # named reduction on reference
  project.query = FALSE,
  features = NULL,
  scale = TRUE,
  npcs = 30,
  l2.norm = TRUE,
  dims = 1:30,
  k.anchor = 5,
  k.filter = NA,
  k.score = 30,
  max.features = 200,
  nn.method = "annoy",
  n.trees = 50,
  eps = 0,
  approx.pca = TRUE,
  mapping.score.k = NULL,
  verbose = TRUE)
```

### Reduction choice (`reduction` arg)

`reduction` controls the basis the anchor search runs IN. Per `?FindTransferAnchors`:

| Value | When to use |
|---|---|
| `"pcaproject"` (default) | scRNA-seq query against scRNA-seq reference. Projects the reference's PCA onto the query, then finds anchors in that projected PCA. |
| `"lsiproject"` | scATAC-seq query against scATAC-seq reference. LSI must already be on the reference; shared features (peaks / bins) must match. See `seurat-rna-atac-integration` recipe. |
| `"rpca"` | Reciprocal PCA projection — reference→query AND query→reference. More conservative than pcaproject; useful when reference and query are technology-different. |
| `"cca"` | Run CCA on the joint reference+query. Slower; can find anchors when pcaproject can't (e.g. very different protocols). |

The SKILL.md body uses the default `"pcaproject"` for the standard
scRNA-seq case.

### The PC-related args

`reference.reduction` — names which reduction on the reference holds
the PCA to project. Must exist on the reference object. The recipe
uses `"pca"`; if you have an integrated reference, use
`"integrated.cca"` / `"integrated.rpca"` / `"harmony"` as appropriate.

`npcs = 30` — when `reference.reduction = NULL`, FindTransferAnchors
computes a fresh PCA on the reference with this many PCs. **When you
pass `reference.reduction = "pca"`, `npcs` is IGNORED** — the existing
reduction is used as-is and its dim count is whatever
`reference[["pca"]]@cell.embeddings` has.

`dims = 1:30` — which PCs feed the anchor finder. Match the reference's
chosen dim count (the same number you used to build the reference's
UMAP).

`approx.pca = TRUE` — uses irlba for the PCA when one is computed; set
to FALSE for exact PCA on small references (<5000 cells).

### Anchor-shape args (the knobs to tune)

| Arg | Default | What it controls |
|---|---|---|
| `k.anchor` | 5 | k-nearest-neighbors for the initial anchor candidates. Higher = more anchors, looser correspondence. |
| `k.filter` | NA (disabled) | If set, filters anchors whose query partner isn't in the reference's k-nearest neighbors. Tightens anchor quality. |
| `k.score` | 30 | k for the local anchor scoring. |
| `max.features` | 200 | Cap on number of features used in the scoring step. |
| `n.trees` | 50 | Annoy index trees for kNN. More = better recall, slower. |
| `mapping.score.k` | NULL | If set, computes per-query-cell mapping.score (similarity to the local reference neighborhood). |

For a typical 3000-cell reference vs 3000-cell query, defaults yield
roughly **500–2000 anchors**.

### Low-anchor diagnoses

| Anchor count | Likely cause | Patch |
|---|---|---|
| <50 | Gene-name mismatch (symbol vs Ensembl) | Run `length(intersect(rownames(reference), rownames(query)))` — should be >5000 |
| 50–200 | Very different cell-type compositions, or one is technology-different | Try `reduction = "cca"` or `reduction = "rpca"` |
| 200–500 | Small reference + small query, working as expected | Proceed; expect noisier predictions |
| 500+ | Healthy | Standard outcome |

If `reduction = "pcaproject"` returns <100 anchors, switching to
`reduction = "cca"` usually rescues — at higher cost.

## `TransferData` — arg semantics

Per `?TransferData` (Seurat 5.5.0):

```
TransferData(anchorset, refdata,
  reference = NULL, query = NULL, query.assay = NULL,
  weight.reduction = "pcaproject",       # the reduction used for weighting
  l2.norm = FALSE,
  dims = NULL,                           # NULL = reuse anchor dims
  k.weight = 50,                         # neighbors for anchor weighting
  sd.weight = 1,                         # Gaussian kernel bandwidth
  eps = 0, n.trees = 50,
  verbose = TRUE, slot = "data",
  prediction.assay = FALSE,
  only.weights = FALSE,
  store.weights = TRUE)
```

`refdata` has two valid shapes:

1. **A vector** — names are reference cell barcodes, values are labels
   (or numeric). The SKILL.md uses
   `refdata = reference[[LABEL_COL, drop = TRUE]]`.
2. **A matrix** — column names are reference cell barcodes; rows are
   features. Transfers continuous expression rather than categorical
   labels. Useful for projecting reference gene-module scores onto the
   query.

`k.weight = 50` — for each query cell, how many anchors weigh in to its
predicted label. Higher = smoother predictions, lower confidence range;
lower = sharper predictions but noisier on cells far from any anchor.

`prediction.assay = TRUE` — instead of a data.frame, returns the
predictions as a new Assay on the query (one row per reference label).
Lets you `FeaturePlot` the per-cell score for any label.

## `MapQuery` — wraps three sub-calls

Per `?MapQuery` (Seurat 5.5.0):

```
MapQuery(anchorset, query, reference,
  refdata = NULL,                        # named list or vector
  new.reduction.name = NULL,
  reference.reduction = NULL,
  reference.dims = NULL,
  query.dims = NULL,
  store.weights = FALSE,
  reduction.model = NULL,                # the reference's UMAP DimReduc (with model)
  transferdata.args = list(),            # passthrough to TransferData
  integrateembeddings.args = list(),     # passthrough to IntegrateEmbeddings
  projectumap.args = list(),             # passthrough to ProjectUMAP
  verbose = TRUE)
```

`refdata` as a **named list** (the SKILL.md's pattern):

```r
refdata = list(celltype = LABEL_COL)
```

- The **list NAME** (`celltype`) becomes the column-name suffix on the
  query: `query$predicted.celltype` + `query$predicted.celltype.score`.
- The **list VALUE** (`LABEL_COL`, which is `"celltype"`) names the
  reference column to transfer FROM.
- These are independent — `refdata = list(myname = "celltype")` produces
  `query$predicted.myname`.

Multi-column transfer:

```r
# Transfer two columns in one call. Output: predicted.l1 + predicted.l2 + their .score
query <- MapQuery(
  anchorset = anchors, reference = reference, query = query,
  refdata = list(
    l1 = "celltype.l1",
    l2 = "celltype.l2"
  ),
  reference.reduction = "pca", reduction.model = "umap",
  verbose = FALSE
)
```

`reduction.model` — names the UMAP DimReduc on the **reference** that
has `return.model = TRUE`. The default is `"umap"`; pass the actual
name if the reference's UMAP was saved under a different reduction key
(e.g. `"umap.integrated"`).

### Tuning `MapQuery` via the args lists

```r
query <- MapQuery(
  anchorset           = anchors,
  reference           = reference,
  query               = query,
  refdata             = list(celltype = LABEL_COL),
  reference.reduction = "pca",
  reduction.model     = "umap",
  transferdata.args   = list(k.weight = 30),         # tighter than default 50
  projectumap.args    = list(n.neighbors = 15),      # passes to uwot::umap_transform
  verbose             = FALSE
)
```

## Score interpretation

### `prediction.score.max` (TransferData path)

- Range 0–1.
- A query cell's max score across all reference labels.
- **High (≥0.8)** — confident assignment. Most query cells in a
  well-matched reference fall here.
- **Mid (0.4–0.8)** — ambiguous. Likely the cell sits between two
  reference labels (e.g. a transitioning population) or the reference
  doesn't have the right label for this cell.
- **Low (<0.4)** — query cell type is poorly represented in the
  reference. Either the reference's label vocabulary is too coarse, or
  the cell is something the reference doesn't contain (e.g. tumor cells
  in a healthy-tissue reference).

A healthy mapping has a **bimodal** score distribution: a high-confidence
peak near 1.0 + a tail of mid/low-confidence cells. A flat distribution
centered ~0.3 means the reference doesn't cover the query.

### `predicted.celltype.score` (MapQuery path)

Same idea, but suffixed with the list-name from `refdata`. Same
interpretation.

### `mapping.score` (when `mapping.score.k` is set in FindTransferAnchors)

Per-query-cell similarity to its local reference neighborhood — a
**different** signal from `prediction.score.max`. Low `mapping.score`
means the query cell lands in a region of the reference's embedding
space that has few reference cells nearby (extrapolation), regardless
of the label confidence.

Use both together:

| `prediction.score.max` | `mapping.score` | Interpretation |
|---|---|---|
| High | High | Strong assignment in a well-populated reference region. Trust. |
| High | Low | Query cell lands in a sparsely-populated region of the reference UMAP; the label may be confidently assigned but the cell is at the edge of what the reference knows. Inspect. |
| Low | High | Query cell is in a reference region but doesn't strongly match any one label — ambiguous biology. |
| Low | Low | Query cell type is alien to the reference. Don't trust either signal. |

## Diagnosing a low-confidence mapping

1. **Check shared features.** `length(intersect(rownames(reference),
   rownames(query)))` — should be >5000. <2000 means a name mismatch.
2. **Check anchor count.** `nrow(anchors@anchors)` — see the table above.
3. **Inspect the score histogram.**
   ```r
   hist(query$predicted.celltype.score, breaks = 50,
        main = "mapping confidence", xlab = "predicted.celltype.score")
   ```
   - Bimodal with a peak near 1.0 → good.
   - Unimodal centered ~0.3 → reference doesn't cover query.
4. **Look at the low-confidence subset.**
   ```r
   low <- query[, query$predicted.celltype.score < 0.5]
   table(predicted = low$predicted.celltype)
   ```
   If low-confidence cells cluster into one or two predicted types, those
   types are likely the reference's nearest-but-not-right labels.
5. **Switch `reduction = "cca"` in `FindTransferAnchors`** for a
   second pass when anchor count is low.
6. **Confirm the reference is the right tissue / organism.** Mapping
   a tumor query against a healthy-tissue reference is technically
   possible but the predicted labels will be misleading.

## When to use Path A (TransferData only) vs Path B (MapQuery)

| Situation | Path |
|---|---|
| Reference has saved UMAP model (`return.model = TRUE`) | B — MapQuery gives you A's output + UMAP projection |
| Reference UMAP lacks saved model | A — TransferData only |
| You want continuous reference expression projected (matrix `refdata`) | A — set `refdata = GetAssayData(reference[["RNA"]])` |
| You want to transfer multiple metadata columns in one call | Either — A takes a named list; B takes a named list directly |
| You only need labels for downstream filtering, no visualization | A — lighter weight |
| You want side-by-side visualization on the same UMAP | B — produces `ref.umap` on the query |

Path B is the standard. Path A exists for cases where the reference
isn't UMAP-model-saved.

## Updating an existing reference with `return.model = TRUE`

If the reference UMAP was built without `return.model = TRUE`:

```r
# Reload the reference and rebuild its UMAP with model saved.
reference <- RunUMAP(reference,
                     dims          = 1:30,
                     reduction     = "pca",
                     return.model  = TRUE,
                     verbose       = FALSE)
saveRDS(reference, "/path/to/labeled_reference.rds")
```

The UMAP coordinates change slightly run-to-run (UMAP is
non-deterministic without a seed), so this re-builds rather than
restoring the prior embedding. Set `seed.use = 42` (or whatever the
original used) if you need a reproducible UMAP across re-runs.
