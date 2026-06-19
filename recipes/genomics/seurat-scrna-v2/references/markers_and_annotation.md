# Cluster markers, DE tests, and annotation

`FindAllMarkers` test choices, the math behind each, marker interpretation
heuristics, and the two annotation paths (manual / reference-mapped).

Load this reference when:
- The user asks "why Wilcoxon?" or "should I use MAST/DESeq2 for markers?"
- A cluster has zero significant markers and you need to diagnose.
- The user asks to annotate clusters and you need to choose manual vs
  automated.
- A marker plot looks unusual (lots of MT genes, all clusters share top
  markers).

## `FindAllMarkers` — per-cluster vs.-rest DE

Source-verified args (Seurat 5.5.0, `args(Seurat:::FindAllMarkers)`):

```r
FindAllMarkers(object,
               assay = NULL,                      # uses DefaultAssay
               features = NULL,                   # NULL = all
               group.by = NULL,                   # NULL = active idents
               logfc.threshold = 0.1,             # recipe overrides to 0.25
               test.use = "wilcox",
               slot = "data",                     # uses log-normalized data
               min.pct = 0.01,                    # recipe overrides to 0.25
               min.diff.pct = -Inf,
               only.pos = FALSE,                  # recipe sets TRUE
               max.cells.per.ident = Inf,
               random.seed = 1,
               latent.vars = NULL,
               min.cells.feature = 3,
               min.cells.group = 3,
               mean.fxn = NULL,
               fc.name = NULL,
               base = 2,                          # log2 fold change
               return.thresh = 0.01,              # post-test p-value floor
               densify = FALSE,
               ...)
```

The recipe overrides three defaults:

| Arg | Default | Recipe value | Why |
|---|---|---|---|
| `logfc.threshold` | `0.1` | `0.25` | Drops weak markers; signal-to-noise gain >> the tail of small-effect genes you lose |
| `min.pct` | `0.01` | `0.25` | Drops genes detected in <25% of either group; markers should be majority-detected |
| `only.pos` | `FALSE` | `TRUE` | Returns only genes UP in the cluster vs. rest. Marker discovery is the asymmetric question "what defines this cluster" |

### `slot = "data"` — what test runs on what

`slot = "data"` uses the log-normalized expression matrix
(`NormalizeData`'s output). For Wilcoxon, t-test, ROC, bimod, LR — fine,
they're rank/distribution-based and log-normalization is appropriate.

For `negbinom` / `poisson`, the matrix needs to be raw counts (model raw
UMIs). Switch to `slot = "counts"` for these tests, or they error.

## `test.use` — six DE tests, when to pick each

Source: Satija lab `FindAllMarkers` reference page + sctransform paper for
the SCT-specific gotcha.

| `test.use` | Statistical model | Required package | When to use |
|---|---|---|---|
| `"wilcox"` (default) | Wilcoxon rank-sum (non-parametric) | none; auto-detects `presto` for ~10× speedup | Default for marker discovery. Robust, no distributional assumptions, fast |
| `"wilcox_limma"` | Wilcoxon via limma | `limma` | For reproducing pre-v4 results — limma's Wilcoxon implementation. Rarely needed |
| `"bimod"` | Likelihood-ratio test of bimodal expression (McDavid 2013) | built-in | Sensitive to genes that are bimodal (off in some cells, on in others) — useful for binary markers |
| `"roc"` | ROC analysis returning AUC per gene | built-in | When you want a single "cleanest marker" per cluster. AUC = 1 means perfect separator. Slow for large datasets |
| `"t"` | Student's t-test | built-in | Parametric, rarely better than Wilcoxon for sparse count data |
| `"negbinom"` | Negative binomial GLM (UMI counts) | built-in | UMI-aware; **requires `slot = "counts"`**. Slower than Wilcoxon |
| `"poisson"` | Poisson GLM (UMI counts) | built-in | UMI-aware, simpler than NB; **requires `slot = "counts"`**. Usually worse than NB on overdispersed data |
| `"LR"` | Logistic regression with likelihood ratio | built-in | Useful when adjusting for `latent.vars` (e.g. batch, donor) — those are added as additional terms |
| `"MAST"` | Hurdle model (zero-inflated + continuous part) | `MAST` (Bioc) | When zero-inflation is suspected and you want to model it explicitly. Slower than Wilcoxon, often slightly more sensitive |
| `"DESeq2"` | DESeq2's negative binomial Wald test | `DESeq2` (Bioc) | Designed for bulk RNA-seq; on single-cell data slow and over-conservative. Use for pseudobulk DE (`seurat-de-testing`), NOT direct per-cell |

### Default = Wilcoxon — why

- No distributional assumptions (counts are sparse + zero-inflated; that
  doesn't fit Gaussian, often doesn't fit NB exactly either).
- Robust to outliers.
- Fast with `presto` (~10× over the native R `wilcox.test`).

### When to pick `MAST` instead

Zero-inflated genes (high % zeros + bimodal continuous distribution) — a
common case in scRNA. MAST models the two parts (proportion detected,
expression conditional on detection) separately. Slightly more sensitive
than Wilcoxon when zero-inflation is strong, similar otherwise.

### When NEVER to use `DESeq2` / `pydeseq2` for per-cell scRNA DE

DESeq2 is built for BULK RNA-seq, where each "sample" is a tissue. On
direct per-cell scRNA-seq, it's:
- Slow (minutes vs. seconds for Wilcoxon).
- Over-conservative (tens to thousands of cells per group break its
  dispersion estimation).
- Not the right tool for the question.

The Seurat docs allow `test.use = "DESeq2"`, but the recommended use of
DESeq2 in scRNA is **pseudobulk DE** — aggregate cells per (cluster, donor)
into pseudo-bulk samples and run DESeq2 on those. See the
`seurat-de-testing` recipe for the pseudobulk path; never call
`FindAllMarkers(test.use = "DESeq2")` directly on real per-cell data.

## SCT assay markers — the `PrepSCTFindMarkers` gate

`FindAllMarkers` on an `SCT` assay requires `PrepSCTFindMarkers()` first.
The gate corrects SCT counts to a shared median UMI across multiple SCT
models (merged objects). For a SINGLE-sample SCT object with one model the
correction is a no-op, but the call is REQUIRED — Seurat errors otherwise.

```r
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)
markers <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE, ...)
```

See the `seurat-sctransform` recipe and `sct_de_gotchas.md` for the full
SCT-vs-RNA DE choice and why some workflows prefer the RNA assay for DE
even when clustering used SCT.

## Marker interpretation heuristics

Per cluster, after `FindAllMarkers`, the recipe prints the top 5 by
`avg_log2FC` (`group_by(cluster) %>% slice_max(avg_log2FC, n = 5)`). Read:

- **`gene`** — symbol. Verify against a curated panel (PanglaoDB,
  CellMarker, a published atlas) for the tissue/organism.
- **`avg_log2FC`** — log2 fold change of the cluster mean vs. rest mean.
  >1.5 is strong; 0.5–1.5 is moderate; <0.5 is weak.
- **`pct.1`** — fraction of cluster cells detecting the gene. Higher is
  better; <50% means the marker isn't broadly expressed in the cluster.
- **`pct.2`** — fraction of OTHER cells detecting the gene. Lower is
  better; <20% means the marker is specific. Genes with `pct.1 ≈ pct.2`
  are not informative even if `p_val_adj` is significant.
- **`p_val_adj`** — Bonferroni-adjusted p-value. The recipe's defaults
  filter at `p_val_adj < 0.01` implicitly via `return.thresh`. Very small
  values (`< 1e-100`) are common for clean lineage markers in PBMC.

A marker with `avg_log2FC = 3, pct.1 = 0.95, pct.2 = 0.05, p_val_adj < 1e-50`
is unambiguous. A marker with `avg_log2FC = 0.4, pct.1 = 0.8, pct.2 = 0.7,
p_val_adj = 0.005` is statistically significant but biologically near-noise.

## Why a cluster might have zero significant markers

If a cluster has 0 genes passing `min.pct = 0.25`, `logfc.threshold = 0.25`,
`p_val_adj < 0.01`:

1. **Under-clustering** — the cluster is heterogeneous, mixing two
   sub-populations whose markers cancel. Bump `resolution` and re-cluster.
2. **Doublet bridge** — the cluster is doublets of two parent populations,
   so markers from both parents partially light up but no clean signature
   dominates. Run a doublet detector.
3. **Cluster too small** — Bonferroni adjustment over ~13k tests is
   punitive at low cell count. Check cluster size; <30 cells often loses
   significance even on strong markers.
4. **Algorithm fragmentation** — Louvain occasionally produces tiny
   fragment clusters from graph topology that don't correspond to biology.
   Visible as a few-cell "cluster" floating off the main UMAP island.

## Annotation — biology, not numbers

Cluster labels come from BIOLOGICAL evidence, not just markers. Two paths:

### Manual annotation

1. Take top 5–10 markers per cluster.
2. Look each up in a curated reference for the tissue:
   - **PBMC / immune**: PanglaoDB, CellMarker, Hao 2021 (Azimuth pbmcref).
   - **Brain**: Allen Brain Atlas, Tabula Muris (mouse), HCA brain.
   - **Gut**: Smillie 2019 colon, HCA gut atlas.
   - **Tumor**: ANY published study on the same tumor type — they're
     never directly comparable but the gross lineage markers transfer.
3. Match each cluster to a canonical cell type. Be conservative — if you
   can't identify a cluster, label it "unknown" rather than guessing.
4. `obj <- RenameIdents(obj, "0" = "CD14+ mono", "1" = "CD4 mem T", ...)`.

DO NOT invent a label from a single marker. "FCGR3A+ monocytes" requires
FCGR3A *plus* the rest of the non-classical-monocyte signature (e.g.
`MS4A7`, `CDKN1C`, `LST1`, low `LYZ`).

### Automated reference mapping

Faster and less subjective when a high-quality reference exists:

| Tool | Mechanism | Strengths |
|---|---|---|
| **SingleR** (Bioc) | Per-cell correlation with reference profiles | Fast, multiple built-in references (HumanPrimaryCellAtlasData, etc.) |
| **Azimuth** | Seurat's reference mapping framework (anchors + label transfer) | Best for PBMC / lung / kidney / heart / motor cortex (Hao 2021 references) |
| **scArches** | Conditional VAE for query-to-reference projection | Best when the reference is itself an integrated atlas (HCA, scvi-tools models) |

These run AFTER clustering (often the cluster labels and the
reference-mapped labels are cross-checked) or even instead of clustering
for atlases where the reference cell-type tree is the canonical
identity system.

See `seurat-reference-mapping` and `seurat-multimodal-reference-mapping`
recipes for the workflow.

## Marker visualization — pick TWO figures, not five

The recipe produces two marker figures by convention:

1. **Dot plot of top 5 per cluster** — broad overview of expression and
   specificity. Each row is a cluster, each column a marker; dot size =
   fraction expressed, dot color = average expression (diverging palette).
2. **FeaturePlot of 4–8 hand-picked canonical lineage markers** — for the
   tissue / organism. PBMC default: `CD3D`, `CD8A`, `MS4A1`, `GNLY`, `LYZ`,
   `PPBP`. ALWAYS REPLACE for non-PBMC samples.

What NOT to do:

- ❌ A wall of 40 FeaturePlots over the full top-5-per-cluster list.
  Overwhelms more than informs.
- ❌ Heatmap of every marker. Use the dot plot — heatmap is too noisy at
  scale and the dot plot already shows expression + specificity.
- ❌ VlnPlot per marker. The dot plot subsumes it; a per-marker violin is
  warranted only when you're investigating ONE marker in detail.

## Latent variable adjustment — `latent.vars`

For `test.use = "LR"` / `"negbinom"` / `"poisson"` / `"MAST"`, you can pass
`latent.vars = c("percent.mt", "donor_id")` to adjust the test for those
covariates. Useful when:
- Differential expression results are confounded by donor effect (multi-
  sample DE without integration).
- Differential expression results are confounded by cell-cycle (cycling
  cells dominate one cluster).

The recipe doesn't use `latent.vars` by default because (a) it doesn't
apply to Wilcoxon (the default test), (b) the recipe is single-sample so
donor adjustment doesn't apply, and (c) latent variables for cluster-
markers are usually a sign that pseudobulk DE is the better approach.
