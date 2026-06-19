# Azimuth references — what exists and how to get one

Catalogue of pre-built multimodal references, the `Azimuth` package wrapper,
and how to obtain reference RDS files when the package install is not
available. Load this when the user has a query but is unsure which
reference to map against, or when the recipe's default
`pbmc_multimodal_2023.rds` cannot be obtained.

Sources consulted: the Azimuth app at `azimuth.hubmapconsortium.org`, the
`Azimuth` R package source at `github.com/satijalab/azimuth`, the per-
reference vignettes linked from the app's "Reference Files" tab, and
Hao et al. *Cell* 184, 3573–3587 (2021) for the PBMC 162k reference's
construction.

## The Azimuth reference catalogue

Each reference has a tissue scope, an organism, a gene-symbol convention,
and a label hierarchy. Mapping a query against the wrong reference yields
high-confidence wrong labels — these scopes matter.

| Reference ID (Azimuth) | Tissue | Organism | Modalities | Labels | Size |
|---|---|---|---|---|---|
| `pbmcref` | Peripheral blood (CITE-seq) | Human | RNA + ADT (228 panel) | `celltype.l1` (8), `celltype.l2` (30), `celltype.l3` (57) | 162k cells |
| `bonemarrowref` | Bone marrow (CITE-seq) | Human | RNA + ADT (25 panel) | `celltype.l1` (5), `celltype.l2` (22) | 30k cells (bmcite) |
| `fetusref` | Fetal whole-body | Human | RNA only | tissue + cell type | 380k cells |
| `kidneyref` | Adult kidney | Human | RNA only | `annotation.l1`, `annotation.l2`, `annotation.l3` | 65k cells |
| `lungref` | Adult lung (HLCA) | Human | RNA only | `ann_level_*` (5 levels) | 580k cells |
| `heartref` | Adult heart | Human | RNA only | `celltype` | 50k cells |
| `motorcortexref` | Adult motor cortex (BICCN) | Human | RNA only | `subclass` and `cluster` | 70k cells |
| `pancreasref` | Adult pancreas | Human | RNA only | `annotation` | 30k cells |
| `tonsilref` | Tonsil (CITE-seq + ATAC) | Human | RNA + ADT | `subset`, `cell.types` | 200k cells |
| `adiposeref` | Adult adipose | Human | RNA only | `cell.types` | 30k cells |

Notes:
- Only the references with **RNA + ADT** modalities have a useful
  `predicted_ADT` imputation; on RNA-only references, the `refdata` list
  should drop `predicted_ADT`.
- Only references with both `spca` and `wnn.umap` reductions support
  this recipe. RNA-only references generally have `pca` + `umap` only —
  use the `seurat-reference-mapping` recipe (RNA-only mapping) instead.
- The list above is the set of references the Azimuth app exposed as of
  the Hao et al. paper + the additions through ~2024; new ones may have
  appeared. Check `Azimuth::AvailableData()` or
  `azimuth.hubmapconsortium.org` for the current catalogue.

## The `Azimuth::RunAzimuth` wrapper

`RunAzimuth(query, reference = "<id>")` collapses this entire recipe into
one call — fetches the reference if not already cached, runs
`FindTransferAnchors` + `MapQuery` with the per-reference defaults, and
returns the query with `predicted.celltype.l1`/`l2`/`l3` populated.

```r
# install once
if (!requireNamespace("Azimuth", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("satijalab/azimuth")
}

library(Azimuth)
query <- RunAzimuth(query, reference = "pbmcref")
# query now has predicted.celltype.l1, predicted.celltype.l2,
# predicted.celltype.l3 + .score columns, the predicted_ADT assay, and
# the ref.umap reduction.
```

When the wrapper is right vs when to use this recipe:

| Goal | Wrapper or this recipe? |
|---|---|
| Just want labels, no parameter inspection | Wrapper |
| User specifies "use Azimuth" | Wrapper |
| Need to inspect or filter anchors before transfer | This recipe |
| Need to tune `dims` away from the wrapper's default | This recipe |
| Reference is custom-built (not in `AvailableData()`) | This recipe |
| Need only some refdata entries, e.g. l1 + ADT only | This recipe |

## Obtaining a reference RDS file (when Azimuth install fails)

When the `Azimuth` R package can't install (common in sandbox environments
without compiler tools, or with stale BiocManager mirrors), download the
reference RDS directly:

1. **From the Azimuth app website.** Navigate to
   `azimuth.hubmapconsortium.org/references`, pick the tissue, click
   "Reference Files" — the page exposes a Zenodo / cloud link for the
   reference `.rds` (and for the older `.h5seurat` form).
2. **From SeuratData.** `SeuratData::InstallData("bmcite")` installs
   the bone-marrow CITE-seq dataset (the bmcite dataset Hao et al. used
   to derive `bonemarrowref`). NOTE: this is the *raw* dataset, not the
   pre-built reference; build the reference yourself with the
   `seurat-wnn-multimodal` recipe + `RunSPCA(assay = "RNA", graph = "wsnn")`.
3. **From the upstream paper's Zenodo.** Hao et al. 2021 deposited the
   PBMC 162k reference at the Zenodo DOI cited in the paper's methods
   section; download `pbmc_multimodal_2023.rds` (or the older
   `pbmc_multimodal.h5seurat`).

The exact public URL for the standalone `.rds` is not documented in the
vignette body — the most stable path is `Azimuth::AzimuthReference("pbmcref")`,
which both fetches and caches.

## Reference-specific ADT naming conventions

The PBMC 162k reference's ADT panel uses **dash-suffixed names** (e.g.
`"CD3-1"`, `"CD56-1"`, `"CD45RA"`) — a convention from the CITE-seq
panel that includes multiple clones for some targets. Other references
use plain names (`"CD3"`, `"CD56"`):

| Reference | ADT naming convention |
|---|---|
| `pbmcref` (PBMC 162k) | dash-suffixed (`CD3-1`, `CD4-1`, `CD8a`, `CD56-1`) |
| `bonemarrowref` (BM CITE) | plain (`CD3`, `CD4`, `CD8a`, `CD56`, `HLA.DR`) |
| `tonsilref` | mixed — check `rownames(reference[["ADT"]])` |

The recipe's Step 6 uses an `intersect(...)` fallback that handles
mismatches gracefully (drops names not present in the panel), but for
the cleanest result, inspect `rownames(reference[["ADT"]])` after
loading and adapt the `proteins_show` list to the actual names.

## Custom (non-Azimuth) references

A "multimodal reference" in this recipe's sense is any Seurat object
with `spca` + `wnn.umap` reductions. You can build one yourself from
a CITE-seq dataset:

```r
# Build a custom reference from a CITE-seq Seurat object (sketch):
# 1. WNN: ref <- FindMultiModalNeighbors(ref, list("pca", "apca"), list(1:30, 1:18))
# 2. UMAP: ref <- RunUMAP(ref, nn.name = "weighted.nn", reduction.name = "wnn.umap",
#                          return.model = TRUE)
# 3. Clusters (so spca has a graph to optimize against):
#    ref <- FindClusters(ref, graph.name = "wsnn", resolution = 2)
# 4. spca: ref <- RunSPCA(ref, assay = "RNA", graph = "wsnn")
# 5. saveRDS(ref, "my_reference.rds")
```

See the `seurat-wnn-multimodal` recipe for the full WNN construction
recipe and `?RunSPCA` for the supervised PCA call. The key invariant:
the saved reference must carry `spca` AND `wnn.umap` (the latter saved
with `return.model = TRUE`), or `MapQuery` will fail with a "no UMAP
model found" error.

## When the reference is wrong for the query

The recipe's "Decisions" #1 calls out the failure mode: a
tissue-mismatched reference yields confident wrong labels. Cross-check
before you trust the output:

- **Feature overlap.** If `intersect(rownames(query), rownames(reference)) /
  length(rownames(reference)) < 0.8`, the query is using a different
  gene-symbol convention (Ensembl IDs instead of symbols, mouse genes
  instead of human, etc.) — fix the input.
- **Anchor count.** Per
  `references/multimodal_mapping_internals.md`'s anchor heuristic,
  anchors / query cells < 10% is a red flag.
- **Score distribution.** A reference that's a poor match yields
  predicted scores skewed toward 0.5 (or with a long low-confidence
  tail). See `references/prediction_qc.md`.
- **Population presence.** A tissue with cell types the reference
  doesn't have will see those cells confidently mapped to the nearest
  reference type with low score. Inspect the
  `predicted.celltype.l1.score` distribution on `ref.umap` — clustered
  low-score regions are the candidates.
