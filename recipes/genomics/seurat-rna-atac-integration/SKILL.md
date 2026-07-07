---
name: seurat-rna-atac-integration
description: Joint scRNA + scATAC analysis on a 10x Multiome dataset with R/Seurat (v5) + Signac — build ONE Seurat object holding an RNA assay AND a Signac ChromatinAssay, preprocess each modality separately (RNA→SCTransform→PCA; ATAC→TF-IDF→SVD/LSI), fuse them with `FindMultiModalNeighbors` (WNN), cluster + UMAP on the weighted graph, and add a gene-activity matrix + CoveragePlot diagnostics. The RNA+ATAC specialization of the general WNN pattern.
when_to_use: ONE 10x Multiome sample (paired RNA + ATAC from the same cells) where you want a joint embedding + clusters that use both modalities — input is the Cell Ranger ARC `filtered_feature_bc_matrix.h5` plus the `atac_fragments.tsv.gz` (with its `.tbi`). Use THIS when the user names "10x Multiome", "scRNA + scATAC", "joint RNA+ATAC", "Signac WNN", "ChromatinAssay". For the general multimodal WNN pattern (RNA + ADT / CITE-seq / arbitrary pairs) see `seurat-wnn-multimodal`. For ATAC-only analysis without RNA, or peak-calling, see the Signac single-modality vignettes. Multi-sample multiome → run this per sample first, then integrate.
invocation: interactive+batch
requires_tools: [run_r]
capabilities_needed: [Seurat, Signac]
keywords: [multiome, 10x Multiome, scATAC, scATAC-seq, ATAC-seq, scRNA-seq, joint RNA ATAC, multimodal, WNN, weighted nearest neighbors, FindMultiModalNeighbors, ChromatinAssay, CreateChromatinAssay, Signac, RunTFIDF, RunSVD, LSI, GeneActivity, CoveragePlot, LinkPeaks, fragments, EnsDb, BSgenome, hg38, Seurat v5, R]
produces: [umap_rna.png, umap_atac.png, umap_wnn.png, modality_weights.png, coverage_plot.png, multiome_processed.rds, multiome_processed.lstar.zarr]
domain: genomics
source: "Signac 'Joint RNA and ATAC analysis: 10x multiomic' vignette — stuartlab.org/signac/articles/pbmc_multiomic (Signac 1.17.0, Seurat 5.4.0). Stuart et al. Nat Methods 18, 1333–1341 (2021) for ChromatinAssay + LSI design."
---

# Joint scRNA + scATAC analysis on a 10x Multiome sample (Seurat + Signac, WNN)

Single-sample recipe for **10x Multiome**: paired RNA + ATAC measured from the
same nuclei/cells. We build ONE Seurat object that holds an `RNA` assay AND a
Signac `ChromatinAssay`, preprocess each modality on its own track
(RNA→SCTransform→PCA; ATAC→TF-IDF→SVD/LSI), then fuse them with
**Weighted Nearest Neighbors** (`FindMultiModalNeighbors`). The fused graph
drives a joint UMAP + clustering, and we add a gene-activity (ATAC-derived
"pseudo-RNA") matrix and a coverage plot for sanity.

Use this when the user names a multiome dataset (Cell Ranger ARC output) and
wants joint clusters. For the general WNN pattern (RNA + ADT, RNA + any other
modality), the framing is `seurat-wnn-multimodal`. For ATAC alone, see the
Signac single-modality vignettes.

## Bundled references — load on demand

This recipe is self-contained for the standard workflow. For deeper detail
on any aspect, load the matching reference file with `read_file` ONLY when
the task needs it — don't pre-load everything:

- `references/signac_chromatinassay.md` — what a Signac `ChromatinAssay`
  is, how `CreateChromatinAssay` consumes peaks + fragments + annotation,
  the GRanges interplay, the soft-dep on `biovizBase`,
  multimodal-H5 layout from Cell Ranger ARC.
- `references/atac_processing.md` — TF-IDF semantics, SVD / LSI
  construction, why LSI 1 correlates with depth and must be dropped,
  `FindTopFeatures` peak filtering, `NucleosomeSignal` /
  `TSSEnrichment` QC, depth-correlation diagnostic.
- `references/coverage_and_links.md` — `GeneActivity` matrix
  construction + interpretation (ATAC-derived pseudo-RNA), `CoveragePlot`
  layout + interpretation, `LinkPeaks` for peak-to-gene regulatory
  links (and the `RegionStats` / `BSgenome` prerequisites).
- `references/genome_pinning.md` — picking the right EnsDb / BSgenome
  combination per organism + Cell Ranger ARC reference build (hg38 vs
  hg19 vs mm10), seqlevels Style mapping (Ensembl `1` vs UCSC `chr1`)
  and why a mismatch silently breaks downstream calls.
- `references/figure_style.md` — Seurat plot styling shared with the
  rest of the collection (theme, palettes, alpha-poke, save dimensions,
  CoveragePlot exception).

## Install

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("Seurat",      quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("Signac",      quietly = TRUE)) install.packages("Signac")
if (!requireNamespace("ggplot2",     quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr",       quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("cowplot",     quietly = TRUE)) install.packages("cowplot")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# Annotation + sequence packages — Bioconductor. Pick the pair matching the
# Cell Ranger ARC reference (hg38 default for human; see references/genome_pinning.md).
for (pkg in c("EnsDb.Hsapiens.v86",
              "BSgenome.Hsapiens.UCSC.hg38",
              "biovizBase",            # SOFT-dep of Signac — required by GetGRangesFromEnsDb
              "GenomicRanges")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}
stopifnot(packageVersion("Seurat") >= "5.0.0")
stopifnot(packageVersion("Signac") >= "1.13.0")
```

> **Why `biovizBase` is in the list.** `GetGRangesFromEnsDb()` (used in
> Step 1) calls `biovizBase::crunch` — biovizBase is a Signac `Suggests`
> (not `Imports`), so it does NOT auto-pull. Without it, Step 1 errors
> with "Please install biovizBase". Install once.

Attach once per session:

```r
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(ggplot2); library(dplyr); library(cowplot)
  library(GenomicRanges)
  library(EnsDb.Hsapiens.v86)              # human v86 / hg38 — swap for mouse / other build
  library(BSgenome.Hsapiens.UCSC.hg38)     # only needed for RegionStats + LinkPeaks
})
```

`library(Seurat)` does NOT attach `ggplot2`/`dplyr`/`cowplot` — load them by
name. For the rationale behind the EnsDb / BSgenome pinning, see
`references/genome_pinning.md`.

## Decisions to surface up front

Tell the user these are the analysis-defining decisions:

1. **Genome / annotation source** — `EnsDb.Hsapiens.v86` +
   `BSgenome.Hsapiens.UCSC.hg38` for human hg38, or the matching mouse /
   other-build pair. MUST match the Cell Ranger ARC reference the FASTQs
   were aligned to, or peak coordinates won't resolve. See
   `references/genome_pinning.md`.
2. **ATAC peak filter (`FindTopFeatures` cutoff)** — `min.cutoff = 5`
   keeps peaks seen in ≥5 cells (Signac vignette default).
   `min.cutoff = "q0"` keeps all peaks. Tighten on small samples, loosen
   on rare-cell datasets. See `references/atac_processing.md`.
3. **PCA / LSI dimensions for WNN** — RNA uses `1:50`; ATAC uses
   `2:40` because **LSI component 1 correlates with sequencing depth**
   and is dropped from the embedding. Do NOT include LSI 1.
4. **Clustering resolution on the WNN graph** —
   `FindClusters(graph.name = "wsnn", resolution = ...)`. 0.5 is
   moderate; lower → coarser joint clusters, higher → finer.
5. **Mode** — INTERACTIVE (default) when this is the primary analysis;
   BATCH when invoked as one of many sub-runs in a larger plan. The
   orchestrator passes `args="batch"`; the agent declares the mode in
   its `present_plan`.

Figures to show as the analysis proceeds:
- `umap_rna.png` — RNA-only UMAP (PCA-driven)
- `umap_atac.png` — ATAC-only UMAP (LSI-driven)
- `umap_wnn.png` — joint UMAP from the weighted graph
- `modality_weights.png` — per-cell RNA/ATAC weights (which modality dominates where)
- `coverage_plot.png` — genome-browser sanity check at one marker

---

## Step 1 — Load the multiome H5 + build the Seurat + ChromatinAssay object

Read the Cell Ranger ARC `filtered_feature_bc_matrix.h5` (a list with both
`Gene Expression` and `Peaks` slots), then add the ATAC fragments file path
into a `ChromatinAssay`. The fragments file lives next to the H5 (Cell Ranger
ARC writes it as `atac_fragments.tsv.gz` with a tabix `.tbi` index).

```r
# --- Load the paired counts (returns a list: $`Gene Expression` and $Peaks) ---
counts <- Read10X_h5("/path/to/filtered_feature_bc_matrix.h5")

# --- Annotation: EnsDb -> GRanges, with UCSC-style "chr" seqlevels ---
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"            # adds "chr" prefix; matches Cell Ranger ARC bams
genome(annotation) <- "hg38"

# --- Seurat object: RNA first ---
obj <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA",
                          project = "<sample_id>")

# --- ChromatinAssay: peaks x cells + fragments file path + annotation ---
# sep = c(":", "-") parses the rownames "chr1:100-200" into a GRanges.
# fragments is the PATH to the bgzipped fragments .tsv.gz (with .tbi alongside).
obj[["ATAC"]] <- CreateChromatinAssay(
  counts     = counts$Peaks,
  sep        = c(":", "-"),
  fragments  = "/path/to/atac_fragments.tsv.gz",
  annotation = annotation,
  genome     = "hg38"
)
```

**Pitfalls inline:**

- `Read10X_h5` returns a **list** when the H5 has multiple modalities —
  index by name (`counts$`Gene Expression`` / `counts$Peaks`), not by
  position. Multiome H5s name the ATAC matrix `Peaks` (not `ATAC`).
- The fragments file path must point to the **.tsv.gz** (NOT the index).
  The tabix index `.tbi` must live next to it; Signac reads it implicitly.
- `seqlevelsStyle(annotation) <- "UCSC"` is REQUIRED — Cell Ranger ARC
  peaks are `chr1` / UCSC style; raw EnsDb is `1` / Ensembl style.
  Mismatched seqlevels yield empty `GeneActivity` matrices and silent
  `CoveragePlot` failures.

**Report:** RNA assay dims (genes × cells), ATAC assay dims (peaks × cells),
and the number of peaks the ChromatinAssay parsed into `granges(obj[["ATAC"]])`.

```r
cat(sprintf("RNA:  %d genes x %d cells\n", nrow(obj[["RNA"]]), ncol(obj)))
cat(sprintf("ATAC: %d peaks x %d cells\n", nrow(obj[["ATAC"]]), ncol(obj)))
cat(sprintf("ChromatinAssay GRanges: %d intervals (genome=%s)\n",
            length(granges(obj[["ATAC"]])), unique(genome(obj[["ATAC"]]))))
```

If ATAC peaks are 0 or RNA cells != ATAC cells, the H5 was not a multiome
output — stop and re-check the input.

For ChromatinAssay internals (fragment lookup, GRanges semantics, why the
`.tbi` is required), read `references/signac_chromatinassay.md`. For
EnsDb vs BSgenome pinning details, read `references/genome_pinning.md`.

---

## Step 2 — ATAC QC: NucleosomeSignal + TSSEnrichment, then filter

ATAC has modality-specific QC that has no RNA analog. The two metrics that
matter per cell:

- **`NucleosomeSignal`** — ratio of mononucleosomal (147–294bp) to
  nucleosome-free (<147bp) fragments. High values flag low-quality cells.
- **`TSSEnrichment`** — signal at transcription start sites relative to
  flanks. Low values flag low signal-to-noise ATAC.

```r
DefaultAssay(obj) <- "ATAC"
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj)
# Per-cell ATAC fragment counts in peaks (added by CreateChromatinAssay):
#   obj$nCount_ATAC, obj$nFeature_ATAC are already populated.

# RNA QC metrics (familiar from seurat-scrna-v2)
DefaultAssay(obj) <- "RNA"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")   # "^mt-" for mouse
```

Inspect the joint distribution before filtering — multiome QC must be passed
on BOTH modalities, so the right view is per-modality violins side by side.
Read off thresholds from the quantile tables (as in `seurat-scrna-v2` Step 2b)
and apply them with `subset()`:

```r
# Per-dataset thresholds — read from the quantiles, don't hard-code these.
obj <- subset(obj,
              subset = nCount_RNA      > 1000  & nCount_RNA      < 25000 &
                      nCount_ATAC      > 1000  & nCount_ATAC     < 100000 &
                      percent.mt       < 20    &
                      nucleosome_signal < 2    &
                      TSS.enrichment   > 1)
cat(sprintf("Post-QC: %d cells\n", ncol(obj)))
```

**Report:** pre-QC cells, post-QC cells, fraction dropped, and which metric
removed which fraction (it's common for ~10–25% of barcodes in a 10x Multiome
output to fail ATAC QC even though their RNA looks fine — that's expected).
If you lose >40%, the thresholds are too tight or the sample is genuinely
low-quality.

For NucleosomeSignal / TSSEnrichment semantics and threshold rationale, read
`references/atac_processing.md`.

---

## Step 3 — RNA preprocessing: SCTransform → PCA → RNA-only UMAP

Standard `SCTransform`. `SCT` handles normalization + HVG + scale in one call;
the resulting assay (named `SCT`) is what PCA runs on.

```r
DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)

# RNA-only UMAP — give it an explicit reduction name so the WNN UMAP later
# doesn't clobber it.
obj <- RunUMAP(obj, reduction = "pca", dims = 1:50,
               reduction.name = "umap.rna", reduction.key = "rnaUMAP_",
               verbose = FALSE)

p_rna <- DimPlot(obj, reduction = "umap.rna", label = FALSE, pt.size = 0.4) +
  ggtitle("RNA-only UMAP (SCT -> PCA)") +
  theme_cowplot() + coord_fixed() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank())
for (k in seq_along(p_rna$layers)) {
  if (inherits(p_rna$layers[[k]]$geom, "GeomPoint")) {
    p_rna$layers[[k]]$aes_params$alpha <- 0.6
  }
}
ggsave("umap_rna.png", p_rna, width = 7, height = 6.5, dpi = 120, bg = "white")
```

> **SCT note for multiome.** `SCTransform` makes a new `SCT` assay and sets
> `DefaultAssay(obj) <- "SCT"`. We rely on this in Step 7 (`LinkPeaks` uses
> `expression.assay = "SCT"`). Don't switch back to `RNA` mid-recipe.

**Report:** number of PCs computed, fraction of variance captured by the first
50, top loadings on PC1 (a quick "is PC1 biology vs MT/ribo" sanity check).

---

## Step 4 — ATAC preprocessing: TF-IDF → SVD (LSI) → ATAC-only UMAP

ATAC counts are sparse and binary-ish. The canonical reduction is **LSI**
(Latent Semantic Indexing) = TF-IDF normalization + truncated SVD. Then UMAP
on `lsi`, dropping component 1 because it correlates with sequencing depth.

```r
DefaultAssay(obj) <- "ATAC"

# Filter to peaks above a minimum cell count (the Signac default for multiome
# vignettes is min.cutoff = 5; pass "q0" to keep all peaks if you have a
# small sample and can't afford to lose any).
obj <- FindTopFeatures(obj, min.cutoff = 5)
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)        # writes the "lsi" reduction

# DROP LSI component 1 — it correlates with depth, not biology. Use 2:40.
obj <- RunUMAP(obj, reduction = "lsi", dims = 2:40,
               reduction.name = "umap.atac", reduction.key = "atacUMAP_",
               verbose = FALSE)

p_atac <- DimPlot(obj, reduction = "umap.atac", label = FALSE, pt.size = 0.4) +
  ggtitle("ATAC-only UMAP (TF-IDF -> SVD/LSI, dims 2:40)") +
  theme_cowplot() + coord_fixed() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank())
for (k in seq_along(p_atac$layers)) {
  if (inherits(p_atac$layers[[k]]$geom, "GeomPoint")) {
    p_atac$layers[[k]]$aes_params$alpha <- 0.6
  }
}
ggsave("umap_atac.png", p_atac, width = 7, height = 6.5, dpi = 120, bg = "white")
```

**Pitfalls inline:**

- `FindTopFeatures(obj, min.cutoff = 5)` — integer 5 (cells), or string `"q0"`
  to keep everything. Both are accepted; pick once and report it.
- `dims = 2:40` for ATAC LSI — NEVER include component 1. If you accidentally
  do, the UMAP separates by depth/library complexity, not cell type.
- LSI is computed on the **filtered** peak set from `FindTopFeatures` — re-run
  `RunSVD` if you re-filter peaks.

**Report:** number of peaks retained after `FindTopFeatures`, LSI computed,
and the depth-correlation diagnostic (Signac's `DepthCor(obj)` plots
correlation of each LSI component with total counts — confirm component 1
has high correlation and 2+ have low).

For TF-IDF / SVD / LSI theory, why LSI 1 correlates with depth, and how to
inspect `DepthCor`, read `references/atac_processing.md`.

---

## Step 5 — Fuse the modalities with WNN (`FindMultiModalNeighbors`)

The WNN step learns a per-cell weight balancing RNA vs ATAC and builds a
fused k-nearest-neighbors graph stored under `wknn` / `wsnn`.

```r
# Same dims as the per-modality UMAPs above: RNA 1:50, ATAC 2:40 (drop LSI 1).
obj <- FindMultiModalNeighbors(
  object               = obj,
  reduction.list       = list("pca", "lsi"),
  dims.list            = list(1:50, 2:40),
  modality.weight.name = c("RNA.weight", "ATAC.weight"),  # one name per modality
  verbose              = FALSE
)
```

After this call:

- `obj@neighbors$weighted.nn` holds the fused KNN.
- `obj@graphs$wknn` (KNN) and `obj@graphs$wsnn` (shared-NN, for clustering)
  are populated.
- `obj$RNA.weight` and `obj$ATAC.weight` are per-cell weights (sum to 1).

**Pitfalls inline:**

- `reduction.list` takes the **string names** of the reductions stored on
  the object (`"pca"` for RNA, `"lsi"` for ATAC) — not the assay names.
- `modality.weight.name` MUST be a length-`length(reduction.list)` vector
  (one name per modality). Passing a single string triggers an internal
  warning and Seurat silently overrides the names to
  `paste0(assay, ".weight")` — which in a multiome where `DefaultAssay(obj)`
  is `"SCT"` after `SCTransform()` yields `SCT.weight` + `ATAC.weight`,
  NOT the `RNA.weight` you expected. Always pass a length-2 character vector.

---

## Step 6 — Joint UMAP + clustering on the WNN graph

Run UMAP on the `weighted.nn` neighbor object (NOT on a stored reduction).
Cluster using `graph.name = "wsnn"`.

```r
obj <- RunUMAP(obj, nn.name = "weighted.nn",
               reduction.name = "umap.wnn", reduction.key = "wnnUMAP_",
               assay = "RNA", verbose = FALSE)

# Cluster on the WNN shared-NN graph (Louvain default; algorithm = 4 for Leiden
# if leidenalg is configured via reticulate).
obj <- FindClusters(obj, graph.name = "wsnn", resolution = 0.5,
                    algorithm = 1, verbose = FALSE,
                    cluster.name = "wnn_clusters")

p_wnn <- DimPlot(obj, reduction = "umap.wnn", group.by = "wnn_clusters",
                 label = TRUE, repel = TRUE, pt.size = 0.4) +
  ggtitle(sprintf("WNN joint UMAP - Louvain res=0.5 - n=%d cells - %d clusters",
                  ncol(obj), length(unique(obj$wnn_clusters)))) +
  theme_cowplot() + coord_fixed() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank()) +
  NoLegend()              # AFTER theme_cowplot — cowplot 1.2.0 overrides legend.position
for (k in seq_along(p_wnn$layers)) {
  if (inherits(p_wnn$layers[[k]]$geom, "GeomPoint")) {
    p_wnn$layers[[k]]$aes_params$alpha <- 0.6
  }
}
ggsave("umap_wnn.png", p_wnn, width = 7, height = 6.5, dpi = 120, bg = "white")
```

**Diagnostic — per-cell modality weight.** Plot `RNA.weight` on the WNN UMAP:
high RNA.weight cells are RNA-dominated (e.g. transcriptionally distinct),
high ATAC.weight cells are ATAC-dominated. Strong spatial structure in this
plot tells you WHICH lineages each modality is informing.

```r
p_w <- FeaturePlot(obj, features = "RNA.weight", reduction = "umap.wnn",
                   pt.size = 0.4, order = TRUE) &
  scale_colour_gradient2(low = "#2166ac", mid = "grey90", high = "#b2182b",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "RNA.weight") &
  theme_cowplot() &
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank()) &
  coord_fixed()
for (i in seq_along(p_w)) {
  pl <- p_w[[i]]
  if (!is.null(pl$layers)) {
    for (k in seq_along(pl$layers)) {
      if (inherits(pl$layers[[k]]$geom, "GeomPoint")) {
        pl$layers[[k]]$aes_params$alpha <- 0.6
      }
    }
    p_w[[i]] <- pl
  }
}
ggsave("modality_weights.png", p_w, width = 7.5, height = 6.5, dpi = 120, bg = "white")
```

**Report:** number of joint clusters, mean RNA.weight (cells <0.5 are
ATAC-dominated), whether the three UMAPs (RNA / ATAC / WNN) show the same
broad structure or whether one modality reveals populations the other misses.

---

## Step 7 — Gene activity matrix (ATAC → pseudo-RNA per gene)

`GeneActivity` sums ATAC fragments overlapping each gene body + 2 kb upstream
into a gene × cell matrix. Useful for (a) checking whether RNA and ATAC agree
on a known marker, and (b) annotating clusters from canonical RNA markers
when looking at ATAC-only populations.

```r
DefaultAssay(obj) <- "ATAC"
gene.activities <- GeneActivity(obj)

# Add it as a new assay alongside RNA + ATAC, then normalize.
obj[["GeneActivity"]] <- CreateAssayObject(counts = gene.activities)
obj <- NormalizeData(obj, assay = "GeneActivity",
                     normalization.method = "LogNormalize",
                     scale.factor = median(obj$nCount_GeneActivity))

cat(sprintf("GeneActivity: %d genes x %d cells\n",
            nrow(obj[["GeneActivity"]]), ncol(obj)))
```

`GeneActivity` returns a `dgCMatrix` (gene × cell, integer counts). Adding it
as a separate assay lets you do `DefaultAssay(obj) <- "GeneActivity"` and then
`FeaturePlot(obj, features = c("CD3D", "MS4A1"))` to compare gene activity
against RNA expression for the same markers — divergence is biology (silenced
locus, primed but not expressed, etc.).

For GeneActivity semantics (window definition, normalization rationale), read
`references/coverage_and_links.md`.

---

## Step 8 — CoveragePlot at one canonical marker (sanity check)

A genome-browser-style ATAC + RNA + gene-track view confirms peak calls land
on real promoters/enhancers and that the chosen clusters show differential
accessibility.

```r
# Pick a marker whose RNA expression is clear in the joint clusters.
# Substitute a real marker for YOUR tissue — examples below are PBMC.
DefaultAssay(obj) <- "ATAC"

p_cov <- CoveragePlot(
  object           = obj,
  region           = "MS4A1",                     # gene symbol (B-cell marker)
  features         = "MS4A1",                     # also overlay RNA expression
  expression.assay = "SCT",
  extend.upstream  = 2000,
  extend.downstream = 2000
)
ggsave("coverage_plot.png", p_cov, width = 9, height = 6.5, dpi = 120, bg = "white")
```

Signac's `CoveragePlot` returns a patchwork of stacked panels (per-cluster
coverage + gene-model track + optional expression dot). Don't add cowplot
theming — Signac's defaults are genome-browser-tuned.

**Report:** which marker you plotted, which cluster shows peak accessibility,
whether RNA expression in the same cluster agrees with ATAC openness.

> **Optional — LinkPeaks for peak-to-gene links.** If you want to score
> regulatory peaks for a gene, compute sequence GC content first
> (`RegionStats` needs the `BSgenome`), then call `LinkPeaks`:
> ```r
> obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)
> obj <- LinkPeaks(obj, peak.assay = "ATAC", expression.assay = "SCT",
>                  genes.use = c("MS4A1", "CD3D"))    # the genes you care about
> ```
> Then `CoveragePlot(obj, region = "MS4A1", features = "MS4A1",
> expression.assay = "SCT", peaks.group.by = "wnn_clusters")` will overlay
> link arcs. This is optional — `LinkPeaks` is slow on the full peak set.

For CoveragePlot layout customization, LinkPeaks semantics + interpretation,
and gene-window choices in GeneActivity, read
`references/coverage_and_links.md`.

---

## Step 9 — Save the joint object

```r
saveRDS(obj, "multiome_processed.rds")
cat(sprintf("Wrote multiome_processed.rds (%.1f MB)\n",
            file.info("multiome_processed.rds")$size / 1e6))
```

The saved `.rds` carries the RNA + ATAC + GeneActivity assays, the `pca` /
`lsi` reductions, the three UMAPs (`umap.rna`, `umap.atac`, `umap.wnn`), the
WNN graphs (`wknn`, `wsnn`), per-cell modality weights, and the joint
clusters. A follow-up session can resume at marker calling / annotation /
trajectory work without re-running 1–8.

---

## Batch variant — use INSTEAD of Steps 1–9 when invoked with args="batch"

Branch on `$ARGUMENTS` at the top of the body. In batch mode:

- Skip the per-step UMAP figures, the modality-weights plot, and the coverage
  plot — orchestrator handles cross-sample comparative plotting.
- Skip per-step "Report" blocks; print ONE final summary line.
- Still save `multiome_processed.rds` (the artifact the orchestrator consumes).

```r
counts <- Read10X_h5("/path/to/filtered_feature_bc_matrix.h5")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"; genome(annotation) <- "hg38"

obj <- CreateSeuratObject(counts$`Gene Expression`, assay = "RNA",
                          project = "<sample_id>")
obj[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks, sep = c(":", "-"),
  fragments = "/path/to/atac_fragments.tsv.gz",
  annotation = annotation, genome = "hg38")

DefaultAssay(obj) <- "ATAC"
obj <- NucleosomeSignal(obj); obj <- TSSEnrichment(obj)
DefaultAssay(obj) <- "RNA"
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Per-dataset QC thresholds — the orchestrator can override per sample.
obj <- subset(obj, subset =
  nCount_RNA  > 1000 & nCount_RNA  < 25000 &
  nCount_ATAC > 1000 & nCount_ATAC < 100000 &
  percent.mt  < 20 & nucleosome_signal < 2 & TSS.enrichment > 1)

DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj, verbose = FALSE)
obj <- RunPCA(obj, verbose = FALSE)

DefaultAssay(obj) <- "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = 5)
obj <- RunTFIDF(obj); obj <- RunSVD(obj)

obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"),
                               dims.list = list(1:50, 2:40),
                               modality.weight.name = c("RNA.weight", "ATAC.weight"),
                               verbose = FALSE)
obj <- RunUMAP(obj, nn.name = "weighted.nn",
               reduction.name = "umap.wnn", reduction.key = "wnnUMAP_",
               assay = "RNA", verbose = FALSE)
obj <- FindClusters(obj, graph.name = "wsnn", resolution = 0.5,
                    algorithm = 1, verbose = FALSE,
                    cluster.name = "wnn_clusters")

ga <- GeneActivity(obj)
obj[["GeneActivity"]] <- CreateAssayObject(counts = ga)
obj <- NormalizeData(obj, assay = "GeneActivity",
                     normalization.method = "LogNormalize",
                     scale.factor = median(obj$nCount_GeneActivity))

saveRDS(obj, "multiome_processed.rds")
cat(sprintf("batch ok | %d cells | %d RNA genes | %d ATAC peaks | %d clusters\n",
            ncol(obj), nrow(obj[["RNA"]]), nrow(obj[["ATAC"]]),
            length(unique(obj$wnn_clusters))))
```

---

## Final response checklist

Summarize:
- input H5 + fragments path, sample ID, and genome / annotation pinning
- RNA cells × genes, ATAC cells × peaks, GeneActivity cells × genes
- QC failure rate and which metric (RNA, ATAC, MT, nucleosome, TSS) drove most drops
- RNA processing: SCT used, PCs computed, top PC1 loadings
- ATAC processing: peak filter (`min.cutoff`), LSI dims used (2:40 — explain WHY LSI 1 dropped)
- WNN: RNA.weight median / range, joint cluster count and sizes
- whether the three UMAPs (RNA / ATAC / WNN) agree or whether one reveals
  structure the other misses
- coverage plot interpretation: marker chosen, cluster showing accessibility,
  whether RNA and ATAC agree
- saved file: `multiome_processed.rds`
- caveats: depth-correlation in LSI, modality-weight imbalance, GeneActivity
  is a proxy not real RNA, peak-gene linkage requires `LinkPeaks` + BSgenome

## See also

- **`seurat-wnn-multimodal`** — the general WNN pattern (RNA + ADT / CITE-seq /
  any two modalities). RNA+ATAC is one instance; switch to that recipe when
  the second modality is not ATAC.
- **`seurat-scrna-v2`** — single-sample RNA-only QC + clustering. Use it
  per-sample on a multiome's RNA modality if you want to sanity-check the RNA
  view independently first.
- **`seurat-de-testing`** — differential testing for cluster markers on the
  WNN joint clusters (FindMarkers / FindAllMarkers, with pseudobulk-DESeq2 for
  multi-sample condition effects).

## Offer an interactive view

Write a viewer-optimized store DIRECTLY from the live Seurat object with lstar (pure R,
highest fidelity — do NOT route through `.h5ad`) and **proactively offer to open it** (a
required part of delivering the result):
```r
DefaultAssay(obj) <- "RNA"        # RNA expression view; WNN clusters + modality weights ride along
d <- lstar::read_seurat(obj)
lstar::lstar_write_viewer(d, "multiome_processed.lstar.zarr")   # precomputes DE / HVGs /
                                                               # cell-major counts (optimized)
```
Then call `open_viewer(file_path="multiome_processed.lstar.zarr")` and present the returned
link so the user can explore the joint RNA+ATAC clusters on the WNN UMAP in pagoda3 — it
opens instantly (pre-optimized, no on-launch conversion). If `open_viewer` returns
`ok:false`, relay the error rather than a dead link. Format / sharing →
**`scrna-viewing-and-interchange`**.
