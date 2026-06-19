# Installation, dependencies, and input I/O for Seurat v5 scRNA-seq

Reader functions for the three common scRNA-seq input shapes (10x triplet,
10x .h5, generic MTX), what each reader expects on disk, and how the load
filter (`min.cells`, `min.features`) interacts with downstream QC.

Load this reference when:
- The user's input doesn't match the SKILL.md's default 10x CellRanger
  directory shape (renamed flat triplet, GEO download, custom MTX, AnnData
  bridge).
- A `Read10X` / `Read10X_h5` / `ReadMtx` call errored and you need to
  trace why.
- You need to confirm a transitive dependency (`Matrix`, `hdf5r`, `dplyr`).

## Install — what each package is for

Verified against the env at validation time (Seurat 5.5.0, sctransform 0.4.3,
ggplot2 + dplyr pulled by Seurat, cowplot 1.2.0, patchwork 1.3.2, ggrepel
0.9.8, presto 1.0.0, MAST 1.32.0, DESeq2 1.46.0).

| Package | Purpose | When `library()` actually attaches it |
|---|---|---|
| `Seurat` | Core S4 classes (Seurat, Assay5, SCTAssay), `Read10X*`, `CreateSeuratObject`, `NormalizeData`, `FindVariableFeatures`, `ScaleData`, `RunPCA`, `FindNeighbors`, `FindClusters`, `RunUMAP`, `FindAllMarkers`, `DimPlot`, `DotPlot`, `FeaturePlot` | always |
| `ggplot2` | Underlying plotting layer; `ggsave`, `theme`, `scale_*` | NOT auto-attached by Seurat — `library(ggplot2)` explicitly |
| `dplyr` | `%>%`, `group_by`, `slice_max`, `select` for the marker post-processing block | NOT auto-attached — `library(dplyr)` explicitly |
| `cowplot` | `theme_cowplot()` — the canonical figure theme across this collection | NOT a Seurat dep — `library(cowplot)` explicitly |
| `patchwork` | `|`, `/`, `+`, `plot_layout`, `plot_annotation` operators for composing the QC scatter pair | Imported by Seurat but not attached — `library(patchwork)` explicitly |
| `tidyr` | `pivot_longer` for reshaping the QC metric metadata into long form for `facet_wrap` | Attached by `library(Seurat)` indirectly via `library(dplyr)`'s deps, but safer to attach explicitly |
| `ggrepel` | `geom_text_repel`, `LabelPoints` non-overlapping labels for HVG plot. Transitive dep of Seurat's `LabelPoints` | NOT a hard Seurat dep — guard with `requireNamespace` before the HVG plot |
| `presto` | C++-accelerated Wilcoxon via `wilcoxauc()`. Seurat v5 auto-detects it for `FindAllMarkers(test.use="wilcox")` and switches transparently (~10× faster, identical results) | NOT attached; never required (Wilcoxon falls back to native R if absent) |
| `Matrix` | Sparse matrix backing (`dgCMatrix`) — the type returned by every reader | Attached as a Seurat dep |
| `hdf5r` | Required by `Read10X_h5` ONLY. Bioconductor-flavored HDF5 binding (NOT `rhdf5` — those are sibling packages with different APIs) | Imported on demand by `Read10X_h5` |
| `qs` | OPTIONAL — `qs::qsave` / `qs::qread` for large-object serialization (2–3× faster than `saveRDS`, smaller files). Pure quality-of-life for downstream sessions | NOT a recipe dep |

The recipe's `library()` block deliberately attaches only the names the body
calls directly, so a missing attach surfaces immediately at load time rather
than 30 cells later at first use.

## Reader function signatures (source-verified against Seurat 5.5.0)

```r
Read10X(data.dir, gene.column = 2, cell.column = 1,
        unique.features = TRUE, strip.suffix = FALSE)

Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)

ReadMtx(mtx, cells, features,
        cell.column = 1, feature.column = 2,
        cell.sep = "\t", feature.sep = "\t",
        skip.cell = 0, skip.feature = 0,
        mtx.transpose = FALSE, unique.features = TRUE, strip.suffix = FALSE)
```

### `Read10X(data.dir = …)` — 10x CellRanger directory triplet

Expects a directory with:
- `barcodes.tsv[.gz]`
- `features.tsv[.gz]` (CellRanger v3+) OR `genes.tsv[.gz]` (CellRanger v2)
- `matrix.mtx[.gz]`

Auto-detects gzip from the `.gz` extension and v2 vs v3 from the feature
filename. Returns a `dgCMatrix` (genes × cells) when the features file has
≤2 columns, or a NAMED LIST of matrices (one per modality: `"Gene Expression"`,
`"Antibody Capture"`, `"Peaks"`, …) when v3 features.tsv has a 3rd column
naming the modality. For an RNA-only run, the matrix is bare; for multimodal
runs, subset with `counts[["Gene Expression"]]`.

**`gene.column = 2`** is the column of `features.tsv` to use as rownames.
Column 1 is Ensembl ID; column 2 is HGNC/symbol; column 3 (when present) is
the modality label. Use `gene.column = 1` only when you need Ensembl IDs
downstream (e.g. for matching to a reference annotation that's Ensembl-keyed).

**`cell.column = 1`** — the barcodes file has one column; rarely changed.

**`strip.suffix = FALSE`** — by default the trailing `-1` on 10x barcodes is
kept. Set `TRUE` to drop it (rarely useful; the suffix is meaningful in
multi-sample contexts).

### `Read10X_h5(filename = …)` — single HDF5 file

10x's `filtered_feature_bc_matrix.h5` (or `raw_feature_bc_matrix.h5`).
`use.names = TRUE` returns gene symbols (the `feature_names` HDF5 dataset);
`FALSE` returns Ensembl IDs.

Like `Read10X`, returns a NAMED LIST for multimodal runs. The `.h5` is faster
to load than a triplet and smaller on disk; preferred when both are available.

### `ReadMtx(mtx = …, cells = …, features = …)` — renamed / GEO-style triplets

Use when files are NOT named the canonical `barcodes.tsv` / `features.tsv` /
`matrix.mtx` (GEO uploads often prepend the GSM accession, e.g.
`GSM5746259_matrix.mtx.gz`). Pass each path explicitly.

**`feature.column = 2`** parallels `Read10X`'s `gene.column = 2` (HGNC by
default). If your features file has only one column (rare CellRanger v1
shape), pass `feature.column = 1`.

**`cell.column = 1`** — barcodes file column. **`cell.sep` and `feature.sep`
default to `"\t"`.** For comma-delimited files, pass `cell.sep = ","`.

**`skip.cell` / `skip.feature` = 0** — number of header lines to drop. Set
these only if the files carry headers (10x's don't).

### Reading from `.h5ad` or `.h5Seurat` — interop bridges

Native Seurat readers do NOT cover AnnData (`.h5ad`) directly. The bridges:

```r
# h5Seurat (legacy SeuratDisk format)
SeuratDisk::LoadH5Seurat("/path/to/sample.h5Seurat")

# h5ad via SingleCellExperiment + zellkonverter (Bioconductor)
sce <- zellkonverter::readH5AD("/path/to/sample.h5ad")
obj <- Seurat::as.Seurat(sce, counts = "counts", data = "logcounts")
```

The roundtrip is lossy — `obsm` reductions, `varm`, `obsp` graphs, etc.
sometimes don't survive. Document which slots you keep.

## Empty-droplet load filter — `min.cells` and `min.features`

`CreateSeuratObject` takes two early-stage filters that are NOT QC and
should NOT be confused with the `nFeature_RNA` / `percent.mt` thresholds:

| Argument | Default | What it removes | Why it's NOT QC |
|---|---|---|---|
| `min.cells = 3` | `0` | Genes detected in fewer than 3 cells | Pure sparsity — you cannot compute variance on a gene seen 1× across the dataset. Has no biological signal either way. |
| `min.features = 200` | `0` | Cells with fewer than 200 unique genes | Empty-droplet floor — 10x runs include ~5–20% empty/ambient droplets that the cell calling missed. Keeping them inflates QC distributions and pulls thresholds toward noise. |

The recipe uses `min.cells = 3, min.features = 200` to prune obvious junk
BEFORE Step 2 computes QC metrics. The real QC (Step 3's `subset()`) reads
the post-load distribution and applies tissue-specific cutoffs.

Setting `min.cells = 0, min.features = 0` makes EVERY gene and EVERY barcode
available — fine for `Read10X_h5` outputs from CellRanger filtered matrices
(which already drop empties), but adds noise to raw matrices.

## Organism-specific MT prefix — case sensitivity matters

`PercentageFeatureSet(obj, pattern = "<regex>", col.name = "percent.mt")`
greps gene names against `pattern`. The pattern IS case-sensitive.

| Organism | MT prefix | Ribosomal pattern |
|---|---|---|
| Human (HGNC) | `^MT-` | `^RP[SL]` |
| Mouse (MGI) | `^mt-` | `^Rp[sl]` |
| Rat (RGD) | `^Mt-` | `^Rp[sl]` |
| Zebrafish (ZFIN) | `^mt-` | `^rp[sl]` |
| Drosophila (FlyBase) | `^mt:` | `^Rp[SL]` |
| C. elegans (WormBase) | `^MTCE\.` (cluster prefix) | `^rpl-` / `^rps-` |

Sanity-check at load time:
```r
n_mt <- sum(grepl("^MT-", rownames(obj)))
stopifnot(n_mt > 0)  # if 0, pattern is wrong; percent.mt would be uniform 0
```

If `n_mt == 0` and `percent.mt` is computed anyway, the metric is uniformly
zero and dying cells slip through QC. Always confirm the prefix matches
before computing.

## Gene name normalization gotchas

Seurat names gene rows by what's in the features file (column 2 if
`gene.column = 2`). Common issues:

- **Underscores → dashes** — `CreateSeuratObject` emits a warning and
  replaces `_` with `-` in gene names because Seurat's downstream parsing
  uses `_` as a delimiter. The original symbol is lost from the object;
  if you need it, save the pre-rename rownames before `CreateSeuratObject`.
- **Duplicate symbols** — `unique.features = TRUE` (default) appends `.1`,
  `.2`, … to collisions. Some Ensembl annotations have multiple Ensembl
  IDs mapping to the same HGNC symbol; this produces silent duplicates.
- **Ensembl vs symbol mismatch with reference data** — if you used
  `gene.column = 1` (Ensembl) for the loader but your marker panel /
  reference is in HGNC, downstream subset / FeaturePlot calls return
  empty. Pick one ID system and stick with it through the workflow.

## Saving / loading the processed object

```r
saveRDS(obj, "seurat_processed.rds")
obj <- readRDS("seurat_processed.rds")
```

`.rds` is the canonical Seurat serialization. Carries all assays (RNA,
SCT if present), all reductions (pca, umap, tsne), all graphs (RNA_nn /
RNA_snn), every metadata column, every layer in every assay. Fully
round-trip-safe.

`.qs` (the `qs` package) is a drop-in 2–3× faster alternative; useful when
the object is >500 MB.
