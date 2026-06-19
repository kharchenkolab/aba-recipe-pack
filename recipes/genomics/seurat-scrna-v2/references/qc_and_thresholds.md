# QC metrics, thresholds, and failure modes

QC metric definitions, where to read thresholds from the data (quantile
tables, NOT prior datasets' defaults), tissue-specific cutoff conventions,
and the doublet-detection options the recipe punts to.

Load this reference when:
- QC removed too many cells (>20%) or too few (<2%) and you need to defend
  the threshold choice.
- The user asks "why these cutoffs?"
- You suspect the threshold is intercepting the wrong population
  (e.g. dropping high-RNA biology mistaken for doublets).
- The user asks about doublet detection.

## QC metric definitions

Each metric is computed BEFORE QC filtering (Step 2 of the recipe), then
the recipe filters on the trio in one `subset()` call.

| Metric | Formula | Failure mode it detects |
|---|---|---|
| `nFeature_RNA` | Number of unique genes detected in the cell (`colSums(counts > 0)`) | LOW: empty / ambient droplet, dead cell; HIGH: doublet |
| `nCount_RNA` | Total UMI count in the cell (`colSums(counts)`) | LOW: undersequenced cell, ambient; HIGH: doublet, hyper-active rare cell |
| `percent.mt` | `100 * colSums(counts[mt_genes, ]) / colSums(counts)` | HIGH: dying / lysed cell with cytosolic mRNA leaked but mitochondrial intact |
| `percent.ribo` | Same formula with `^RP[SL]` (human) genes | Variable — high values often biological (T cells, plasma cells); rarely useful as a filter |

`nFeature_RNA` and `nCount_RNA` are computed automatically by
`CreateSeuratObject` and live in `obj$nFeature_RNA` / `obj$nCount_RNA`.
The MT / ribo metrics are added explicitly via `PercentageFeatureSet`.

## Reading thresholds off the quantile tables (do this, not "use defaults")

The recipe prints two quantile tables per metric:
- `qs_lo = c(0.005, 0.01, 0.025, 0.05, 0.10)` — for picking floors.
- `qs_hi = c(0.50, 0.75, 0.90, 0.95, 0.975, 0.99)` — for picking ceilings.

### `nFeature_RNA` floor

Pick from the 1st–5th percentile of `qs_lo`. That's where the empty-droplet
shoulder sits — cells with hundreds of genes that are clearly distinct from
the bulk distribution centered in the thousands.

- **PBMC, 10x v3 chemistry**: floor ~200 (vignette default); the 5th
  percentile usually lands around 600–800 — anything below 200 is almost
  certainly debris.
- **Neutrophils, eosinophils, low-RNA cells**: floor ~100; biology is
  genuinely sparse. A 200 floor erases an entire population.
- **Erythrocytes**: floor ~50; even lower. These cells barely have nuclei.
- **Neurons, hepatocytes, high-RNA cells**: floor 500–1000.

### `nFeature_RNA` ceiling

Pick from the 99th percentile of `qs_hi`. That's where the doublet shoulder
begins. Homotypic doublets (two T cells fused) look like 2× a single T cell;
heterotypic doublets (T + monocyte) have intermediate counts but the gene
count is roughly additive.

- **PBMC, 10x v3, ~3k cells loaded**: ceiling 2500–5000.
- **PBMC, ~10k cells loaded (>5% doublet rate)**: pair the ceiling with a
  proper doublet detector (scDblFinder / DoubletFinder, see below).
- **Neurons, plasma cells, high-RNA cells**: ceiling 8000–12000.

The ceiling is the noisier choice. If the 99th percentile is much higher
than 2× the median of `nFeature_RNA`, you have either real high-RNA biology
or a heavy doublet load. The pre-filter scatter (`qc_scatters_pre.png`)
shows which.

### `percent.mt` ceiling

Pick from `qs_hi`. Healthy cells have <5% MT reads; dying cells have ~20–80%.
Tissue-specific guidance:

| Tissue / sample | Typical `percent.mt` ceiling |
|---|---|
| PBMC, fresh blood | 5–10 |
| Solid tumor, dissociated | 15–25 |
| Brain (cortical neurons) | 15–20 |
| Single-nucleus (snRNA-seq) | 5 (high MT = nucleus contaminated with cytosol) |
| Frozen tissue (variable) | 20–30 |

The PBMC3k vignette uses `< 5`; the recipe defaults to `15` because PBMC3k
is unusual in having extremely clean MT distribution and the recipe is
calibrated for the broader genomics-recipe population.

### `percent.ribo` — usually skip the filter

Ribosomal protein gene fraction is bimodal in T cells, plasma cells, and
some other populations — it reflects translation activity, NOT cell death.
Filtering on it removes biology. Compute it for inspection, do NOT include
it in the `subset()` call unless you see a clearly bimodal failure-mode
shoulder (rare).

## The single-pass `subset()` filter

```r
obj <- subset(obj, subset =
  nFeature_RNA > THRESH$nFeature_low  &
  nFeature_RNA < THRESH$nFeature_high &
  percent.mt   < THRESH$mt_high)
```

Three conjoined conditions, single `subset()` call. Sequential `subset()`
calls work but are slower (each rebuilds the metadata frame). The strict
`>` / `<` (not `>=` / `<=`) matches the vignette convention.

## Diagnostic — % cells removed

Track:
```r
n_before <- ncol(obj); n_after <- ncol(obj_post)
pct_removed <- 100 * (n_before - n_after) / n_before
```

Heuristics:

- **<2% removed**: thresholds may be too loose. Sanity-check by tightening
  one (`percent.mt < 10`) and re-running; if cluster structure improves,
  the original was indeed loose.
- **2–15% removed**: typical for a CellRanger-filtered 10x sample.
- **15–25% removed**: thresholds tight but defensible. Re-read the quantile
  table — make sure each filter targets a different population (otherwise
  you're double-counting one failure mode).
- **>25% removed**: thresholds are wrong, or the dataset has a serious
  quality problem upstream. Return to Step 2.

## Failure-mode diagnostic from the post-filter violin

After filtering, the post-filter violin (`qc_violins_post.png`) should show
**no truncated peak at the threshold**. A flat peak butting up against the
cutoff means the threshold is intercepting a real population (often the
high tail of an active subpopulation). Loosen it.

A smooth distribution that tapers naturally below the cutoff means the
threshold caught the failure mode it was meant to catch.

## Doublet detection — when QC isn't enough

QC catches:
- **Empty droplets** (low `nFeature_RNA` / `nCount_RNA`).
- **Dying cells** (high `percent.mt`).
- **Obvious homotypic doublets** (extreme high `nFeature_RNA`).

QC does NOT catch:
- **Heterotypic doublets** — two different cell types fused. Their
  `nFeature_RNA` is intermediate; their transcriptome is a mix of two
  signatures. Look like an "intermediate" cluster between two populations
  on the UMAP. The standard detectors:

| Tool | Mechanism | When |
|---|---|---|
| **scDblFinder** (Bioconductor) | Simulates synthetic doublets from random pairs, trains a classifier per cell. Sensitive + reasonably specific | First-line. Single command: `scDblFinder(sce)` |
| **DoubletFinder** | Same idea, original implementation | Lower-throughput than scDblFinder; older codebase |
| **scrublet** (Python) | Same simulate-and-classify idea, scanpy ecosystem | Use the Python sibling recipe if you're on Python |

**Order of operations:** run the doublet detector AFTER Step 3's QC, BEFORE
Step 4's normalization. Doublet detectors need raw counts and a roughly
clean cell population — running before QC confuses the simulated-vs-real
boundary; running after normalization throws away the count structure the
detector needs.

The recipe does NOT include doublet detection in the canonical flow because
it adds a Bioconductor dependency and is only needed for >5% expected
doublet rate (10x runs loading >7k cells per channel). Add it explicitly
when:
- The chemistry was high-load (>7k cells/channel).
- You see "intermediate" clusters on the UMAP that share markers from two
  major populations.
- The downstream analysis is sensitive to spurious cell-type intermediates
  (trajectory inference, RNA velocity, cell-cell communication).

## Cell-cycle scoring — optional and usually after QC

If a downstream PC is dominated by cell cycle (`MKI67`, `TOP2A`, `MCM*`,
`PCNA`), regress it via `CellCycleScoring(obj)` then add `"S.Score"` /
`"G2M.Score"` to `vars.to.regress` in `ScaleData`. The reference panels are
human gene names in Seurat's `cc.genes.updated.2019`; orthologs for other
organisms need explicit translation.

The recipe doesn't include cell-cycle regression by default because (a) most
samples don't need it and (b) it can erase real cycling biology (e.g.
tumor proliferation gradients). Surface as a Step 4 / Step 5 option only
when the PC loadings show the confound.
