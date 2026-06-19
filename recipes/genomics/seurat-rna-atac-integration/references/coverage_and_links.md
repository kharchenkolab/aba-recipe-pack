# GeneActivity, CoveragePlot, LinkPeaks

ATAC-derived gene-activity matrices, genome-browser-style coverage
plots, and peak-to-gene regulatory link inference. Load this when the
agent needs to interpret a GeneActivity matrix, customize a
CoveragePlot layout, or understand what LinkPeaks measures.

Sources consulted: Stuart et al. *Nat Methods* 18 (2021) (the
Signac design paper, which introduced the gene-activity and
peak-to-gene-link machinery); Pliner et al. *Mol Cell* 71 (2018) (the
Cicero paper that informed peak-to-gene correlation methodology);
Signac source `R/objects.R` (GeneActivity), `R/visualization.R`
(CoveragePlot), `R/links.R` (LinkPeaks); Signac vignettes at
`stuartlab.org/signac/articles/pbmc_multiomic` and
`peak_calling`.

## GeneActivity — ATAC fragments → pseudo-RNA matrix

`GeneActivity()` builds a gene × cell matrix from ATAC fragments by
summing fragments that overlap each gene's "activity window": the
gene body PLUS a small upstream/downstream extension. The result
behaves like a coarse RNA count matrix, but reflects chromatin
accessibility (a leading indicator of expression) rather than mRNA
abundance.

```r
DefaultAssay(obj) <- "ATAC"
gene.activities <- GeneActivity(obj,
                                features          = NULL,        # all genes (default)
                                extend.upstream   = 2000,
                                extend.downstream = 0,
                                biotypes          = "protein_coding",
                                max.width         = 5e5,
                                gene.id           = FALSE)         # use gene symbols, not IDs
```

Argument semantics (from Signac source + `?GeneActivity`):

- **`features`** — character vector of gene symbols to compute
  activity for. `NULL` = all genes in the annotation that pass
  `biotypes` + `max.width` filters.
- **`extend.upstream`** — bp upstream of the gene's 5'-most TSS to
  include in the window. Default 2000 (capture promoter signal).
- **`extend.downstream`** — bp downstream of the gene's 3'-most
  position to include. Default 0 (don't extend; the gene body is
  large enough).
- **`biotypes`** — filter to gene biotypes from the EnsDb annotation.
  Default `"protein_coding"`. Set to NULL to include lncRNA, miRNA,
  etc. (the activity matrix becomes much larger).
- **`max.width`** — drop genes wider than this many bp (avoids
  enormously long genes dominating, e.g. DLG2, RBFOX1 at >1Mb).
  Default 5e5 (500kb).
- **`gene.id`** — `FALSE` (default) uses gene symbols as row names;
  `TRUE` uses Ensembl IDs.

The return is a `dgCMatrix` (gene × cell, integer counts).

### Why the matrix is "pseudo-RNA"

The activity score for a gene is the fragment count in its window. If
the window is open (chromatin accessible at the promoter + gene body),
the count is high — and on average, accessible promoters precede or
accompany transcription. So gene activity tracks RNA expression
imperfectly: positively correlated, lagging slightly in time, noisier.

The recipe normalizes it with `LogNormalize` to mimic the RNA
log-normalization workflow:

```r
obj[["GeneActivity"]] <- CreateAssayObject(counts = gene.activities)
obj <- NormalizeData(obj, assay = "GeneActivity",
                     normalization.method = "LogNormalize",
                     scale.factor = median(obj$nCount_GeneActivity))
```

After this, `FeaturePlot(obj, features = c("CD3D", "MS4A1"),
reduction = "umap.wnn", assay = "GeneActivity")` will overlay
gene-activity gradients on the joint UMAP — useful as a sanity check
that the RNA + ATAC modalities agree on canonical lineage markers.

### Comparing GeneActivity to RNA

Divergence between RNA and GeneActivity for the same gene is
informative:

| RNA high, GeneActivity high | Active transcription — chromatin open AND gene expressed |
| RNA low, GeneActivity high | Primed but not transcribed — chromatin open, transcription not engaged (often progenitor populations) |
| RNA high, GeneActivity low | Silenced locus with continued mRNA — rare, suggests mRNA stability or alternative regulation |
| RNA low, GeneActivity low | Inactive locus |

A cluster where many lineage markers fall in the "primed but not
transcribed" quadrant is often a stem / progenitor population — worth
flagging for downstream annotation.

## CoveragePlot — anatomy of the genome-browser panel

`CoveragePlot()` is Signac's flagship visualization: a per-group
ATAC fragment pileup track at a genomic locus, optionally with a gene
model, peak track, RNA expression overlay, and link arcs.

```r
DefaultAssay(obj) <- "ATAC"

p_cov <- CoveragePlot(
  object             = obj,
  region             = "MS4A1",              # gene symbol OR "chr11:60450000-60470000"
  features           = "MS4A1",              # gene to overlay expression for (optional)
  expression.assay   = "SCT",                # which assay for the expression overlay
  group.by           = "wnn_clusters",       # NULL -> Idents(obj)
  idents             = NULL,                 # subset of groups (NULL = all)
  extend.upstream    = 2000,
  extend.downstream  = 2000,
  annotation         = TRUE,                 # show gene model below the coverage
  peaks              = TRUE,                 # show peak track below the gene model
  peaks.group.by     = NULL,                 # group peaks by metadata column
  links              = TRUE,                 # show LinkPeaks arcs (if LinkPeaks was run)
  ymax               = NULL,                 # custom y-axis cap per panel
  window             = 100                   # smoothing window for the pileup (bp)
)
```

The return is a **patchwork** of stacked ggplot panels:

1. **Coverage tracks** (one per group) — fragment density across the
   region.
2. **Gene model track** — exons (boxes), introns (lines), strand
   (arrows). Pulled from `Annotation(obj[["ATAC"]])`.
3. **Peak track** — the called peaks overlapping the region.
4. **Link arcs** (optional) — connect peaks to gene TSS, colored by
   correlation strength. Requires `LinkPeaks` first.
5. **Expression panel** (optional) — violin/box of `features` from
   the `expression.assay` per group, side by side with coverage.

### Customizing CoveragePlot

The recipe's Step 8 says "don't add cowplot theming" — Signac's
default styling is genome-browser tuned. If you DO need to customize:

```r
# Modify individual panels in the patchwork
p_cov[[1]] <- p_cov[[1]] + theme_minimal()    # the first coverage panel
ggsave("coverage_custom.png", p_cov, width = 9, height = 6.5, dpi = 120, bg = "white")
```

Common adjustments:

- **`ymax = N`** — cap the per-panel y-axis. Use when one cluster has
  much higher coverage than others, drowning out others.
- **`window = 50`** (lower) — finer smoothing; useful for narrow
  promoters.
- **`peaks.group.by = "wnn_clusters"`** — when LinkPeaks has been run,
  color peaks by which group they're linked to.

### Common failure modes

- **Blank coverage tracks.** Almost always seqlevels-style mismatch
  (`chr1` vs `1`) — see `references/signac_chromatinassay.md`.
- **"Cannot find region"** — the gene symbol doesn't exist in
  `Annotation(obj[["ATAC"]])`. Check `head(Annotation(obj[["ATAC"]])$gene_name)`.
- **All panels look identical** — group.by didn't split the cells
  (e.g. `wnn_clusters` was numeric and got coerced wrong). Check
  `table(obj$wnn_clusters)`.
- **Expression panel missing** — `expression.assay` doesn't carry the
  feature; check `rownames(obj[["SCT"]])`.

## LinkPeaks — peak-to-gene regulatory links

`LinkPeaks()` scores each peak's regulatory association with nearby
genes by correlating per-cell peak accessibility with per-cell gene
expression, controlling for GC content + peak width (matched-control
distribution).

```r
# Pre-req: compute GC content + width per peak. Needs the BSgenome.
obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# Then score links
obj <- LinkPeaks(
  object           = obj,
  peak.assay       = "ATAC",
  expression.assay = "SCT",
  genes.use        = c("MS4A1", "CD3D"),     # only score these — defaults to all variable genes
  distance         = 5e5,                      # peaks within ±500kb of TSS
  min.cells        = 10,                       # peaks accessible in ≥10 cells
  score_cutoff     = 0.05,                     # Pearson r threshold (after matched-control)
  pvalue_cutoff    = 0.05
)
```

Argument semantics:

- **`peak.assay`** — typically `"ATAC"`.
- **`expression.assay`** — typically `"SCT"` (NOT `"RNA"`, because
  the SCT residuals are what the variance modeling uses).
- **`genes.use`** — restrict to a subset of genes (default = all
  variable features from the expression assay). On full datasets,
  scoring all genes is slow (~1 hour for 3k cells × 80k peaks).
- **`distance`** — max bp between peak and TSS to consider. Default
  5e5 (500kb).
- **`min.cells`** — drop sparse peaks before scoring.
- **`score_cutoff` / `pvalue_cutoff`** — significance thresholds.

The result is stored in `Links(obj[["ATAC"]])` as a GRanges with
columns `score` (the correlation), `pvalue`, `gene`, `peak`. Render
with `CoveragePlot(..., links = TRUE)`.

### Interpreting a link

A positive link (high correlation): the peak's accessibility tracks
the gene's expression across cells. Often interpreted as a
candidate enhancer for that gene — but correlation isn't causation.
Validate with TF-motif analysis or perturbation.

A negative link: accessibility and expression anti-correlate. Less
common; may indicate a silencer-like regulatory element or a peak
shared with a different gene's repression context.

The `score` is a partial correlation controlling for GC + width via
matched-control resampling (see Stuart et al. 2021 Methods §"Peak-gene
linkage").

### When to skip LinkPeaks

The recipe's Step 8 makes it optional because:

- It's slow on full peak sets (1+ hour for 3k cells, longer at scale).
- It requires the `BSgenome` (an extra ~700MB install for human).
- It's only useful when the downstream interpretation needs
  regulatory hypotheses; for cluster annotation alone, GeneActivity +
  CoveragePlot is enough.

Skip unless the user names "enhancer" / "regulatory" / "peak-gene
link" / "LinkPeaks" in their question.
