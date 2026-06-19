# ATAC preprocessing — TF-IDF, SVD/LSI, QC, depth correlation

What `RunTFIDF` + `RunSVD` produce, why LSI component 1 correlates with
sequencing depth, how `FindTopFeatures` filters peaks,
`NucleosomeSignal` / `TSSEnrichment` semantics, and how to inspect the
depth-correlation diagnostic. Load this when the agent needs to defend
the `dims = 2:40` choice, debug a depth-dominated UMAP, or interpret
ATAC QC distributions.

Sources consulted: Cusanovich et al. *Science* 348 (2015) (the original
LSI-on-chromatin paper); Stuart et al. *Nat Methods* 18 (2021) (Signac
design, including the depth-correlation discovery); Signac source
`R/preprocessing.R` (`RunTFIDF`, `RunSVD`); Signac vignettes at
`stuartlab.org/signac/articles/pbmc_vignette` and
`pbmc_multiomic`.

## TF-IDF — what it does and why

ATAC count matrices are dominated by sparsity and uneven depth: most
cells have ≤1 fragment per peak, and per-cell total counts vary by
~10x. Plain log-normalization (good for RNA) doesn't help here because
the cell-by-cell variance is dominated by total depth.

TF-IDF (Term Frequency × Inverse Document Frequency, borrowed from
text-document analysis) normalizes BOTH dimensions:

- **TF (per cell):** divide each peak's count by the cell's total
  count → fraction of cell's fragments in that peak.
- **IDF (per peak):** weight peaks by `log(n_cells / n_cells_with_peak)`
  → up-weight rare peaks, down-weight ubiquitous peaks.
- **Output:** `TFIDF[i,j] = log(1 + TF[i,j] * IDF[i] * scale_factor)`.

In Signac:

```r
obj <- RunTFIDF(obj,
                method       = 1,            # Stuart et al. 2021 default
                scale.factor = 1e4)
```

Methods (Signac source `R/preprocessing.R`):

- `method = 1` (default) — `log(1 + TF * IDF * scale_factor)`. Used by
  the Signac vignettes; matches Cusanovich et al. 2015.
- `method = 2` — `log(1 + TF) * log(1 + IDF * scale_factor)`. Variant
  that handles very rare peaks slightly differently.
- `method = 3` — `log(TF * IDF * scale_factor)` (no `1+`); avoid
  (numerical issues at zeros).
- `method = 4` — `TF * IDF` (no log); deprecated.

Stick with `method = 1` unless you have a paper-replication reason to
pick otherwise. The TFIDF result is stored in `obj[["ATAC"]]@data`.

## SVD → LSI — the ATAC equivalent of PCA

`RunSVD()` runs truncated singular value decomposition on the TFIDF
matrix. The result (`obj@reductions$lsi`) is conceptually like PCA but
the input is TFIDF, not log-normalized counts, and the math is SVD,
not eigendecomposition of the covariance.

```r
obj <- RunSVD(obj,
              assay     = "ATAC",
              reduction = "lsi",
              n         = 50,            # number of components to retain
              scale.embeddings = TRUE)
```

The `lsi` reduction has the same shape as `pca`: cells × components.
Use it as `reduction = "lsi"` everywhere you'd use `"pca"` for RNA.

## Why LSI 1 correlates with sequencing depth

Empirical discovery (Stuart et al. 2021, supplementary fig.): on
nearly every ATAC dataset, the FIRST LSI component is strongly
correlated with total per-cell fragment count. The interpretation:

- TFIDF normalization handles per-cell depth on average, but residual
  signal (from technical variation in capture efficiency, library
  complexity, etc.) still loads heavily on the dominant direction of
  variation.
- That dominant direction IS LSI 1. So LSI 1 ≈ depth.

If you include LSI 1 in downstream dimensionality reduction (UMAP,
clustering), the embedding will separate cells by depth rather than
cell identity. This is the silent-failure mode the recipe's Step 4
warns about.

The fix: drop LSI 1. Use `dims = 2:40` (or `2:N` for some N).

To VERIFY this on your own dataset, use Signac's `DepthCor`:

```r
DepthCor(obj, n = 30)   # plots correlation of each LSI component with total counts
```

The expected plot:
- Component 1: |correlation| > 0.7 (often > 0.9).
- Components 2+: |correlation| < 0.2.

If a later component (e.g. LSI 5) ALSO shows high depth correlation,
that's unusual — drop it too, or investigate.

## `FindTopFeatures` — peak filtering

ATAC peak sets from Cell Ranger ARC are large (50–100k peaks).
Most peaks are accessible in only a handful of cells; LSI on the full
peak set wastes compute on noise. `FindTopFeatures` filters:

```r
obj <- FindTopFeatures(obj, min.cutoff = 5)     # keep peaks accessible in >= 5 cells
obj <- FindTopFeatures(obj, min.cutoff = "q0")  # keep all (no filter)
obj <- FindTopFeatures(obj, min.cutoff = "q5")  # keep top 95% by total counts
```

Options:

- **Integer N** — keep peaks accessible in ≥N cells.
- **`"qX"`** — quantile cutoff on total peak count. `"q0"` = no
  filter, `"q5"` = drop bottom 5%, `"q25"` = drop bottom 25%.

Defaults from Signac vignettes:
- **PBMC multiome vignette (10x):** `min.cutoff = 5`. Small datasets
  (~3k cells); keeps ~70k peaks.
- **PBMC scATAC-only vignette:** `min.cutoff = "q0"`. Different
  rationale; uses depth-correlation as the filter instead.

Pick once and report the choice. The choice affects:
- LSI components (different peak set → different SVD).
- GeneActivity (uses the filtered peak set's overlap with annotation).
- CoveragePlot (only filtered peaks render; the underlying fragments
  file is queried directly so coverage is correct, but per-peak
  annotations come from the filtered set).

## NucleosomeSignal — what it measures

`NucleosomeSignal()` computes a per-cell ratio: mononucleosomal
fragments (147–294 bp) / nucleosome-free fragments (<147 bp). Biology:
nucleosome-free regions are open chromatin; mononucleosomal fragments
are intermediate. The ratio is a quality metric:

- **Healthy cells:** ratio ~0.5–2 (some nucleosome-bound DNA, mostly
  open). Per Signac vignette defaults: filter `nucleosome_signal < 2`.
- **Low-quality cells:** ratio > 4. Often these are cells with
  damaged/fragmented chromatin or low-quality libraries.

Stored in `obj$nucleosome_signal`. Compute with
`obj <- NucleosomeSignal(obj)` — call requires fragments file (it
reads fragment sizes from the tabix-indexed `.tsv.gz`).

## TSSEnrichment — what it measures

`TSSEnrichment()` computes per-cell signal at transcription start
sites relative to flanking regions:

```
TSS.enrichment = mean(coverage in TSS ± 500bp) / mean(coverage in flanks 1-2kb away)
```

Biology: high-quality ATAC cells have strong TSS accessibility (most
genes have open TSS regions). Low TSS enrichment means weak chromatin
structure or low signal-to-noise.

- **Healthy cells:** TSS.enrichment > 2 (typically 5–15).
- **Low-quality cells:** TSS.enrichment < 1.

Stored in `obj$TSS.enrichment`. Compute with
`obj <- TSSEnrichment(obj)` — also requires fragments file.

The per Signac vignette default filter: `TSS.enrichment > 1` (the
recipe uses this). For higher-quality datasets, tighten to > 2 or > 3.

## The QC quadrant — interpreting the four metrics together

The recipe's Step 2 filters cells on FIVE thresholds:
- RNA: `nCount_RNA > 1000 & < 25000`, `percent.mt < 20`
- ATAC: `nCount_ATAC > 1000 & < 100000`, `nucleosome_signal < 2`,
  `TSS.enrichment > 1`

Inspect all five distributions BEFORE setting thresholds — what's
typical for one tissue / depth combo doesn't transfer.

```r
# Per-cell summary
qc_df <- obj@meta.data[, c("nCount_RNA", "nFeature_RNA", "percent.mt",
                            "nCount_ATAC", "nFeature_ATAC",
                            "nucleosome_signal", "TSS.enrichment")]
apply(qc_df, 2, quantile, probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99))
```

The typical 10x Multiome pattern (PBMC unsorted 3k as a calibration
point):

| Metric | 1st pct | 25th | Median | 75th | 99th |
|---|---|---|---|---|---|
| nCount_RNA | 1500 | 4000 | 7000 | 12000 | 30000 |
| percent.mt | 1 | 4 | 8 | 14 | 35 |
| nCount_ATAC | 2000 | 8000 | 18000 | 40000 | 120000 |
| nucleosome_signal | 0.3 | 0.6 | 0.8 | 1.2 | 4 |
| TSS.enrichment | 1.5 | 4 | 6 | 9 | 16 |

The "lose ~10–25% on ATAC even though RNA looks fine" pattern the
recipe documents reflects the higher technical failure rate on the
ATAC side — open-chromatin assays are more sensitive to fixation /
permeabilization than RNA.

## Putting it together — the canonical Step 4 in detail

```r
DefaultAssay(obj) <- "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = 5)    # filter peaks
obj <- RunTFIDF(obj)                           # normalize
obj <- RunSVD(obj)                             # reduce to lsi
# Verify depth correlation:
print(DepthCor(obj, n = 10))                    # confirm LSI 1 is depth-correlated
# Drop LSI 1:
obj <- RunUMAP(obj, reduction = "lsi", dims = 2:40,
               reduction.name = "umap.atac",
               reduction.key  = "atacUMAP_")
```

The five lines above are the entire ATAC dimensionality-reduction
pipeline; everything downstream (clustering, WNN, GeneActivity) reads
from `obj[["ATAC"]]@data` (TFIDF) or `obj@reductions$lsi`.
