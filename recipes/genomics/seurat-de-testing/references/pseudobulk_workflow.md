# Pseudobulk workflow — AggregateExpression + DESeq2 / edgeR

The deep methodology behind Path C of SKILL.md. Covers what
`AggregateExpression` is actually doing, how to lay out the donor / sample
/ cell-type structure correctly, the limitation of the
`FindMarkers(test.use = "DESeq2")` wrapper for non-trivial designs, and how
to drop to native DESeq2 / edgeR for multi-factor / paired / interaction
designs.

Load this reference when Path C is in play AND any of:
- You need to verify that aggregation actually summed counts (not averaged).
- The design is not the simple two-level `~group` (i.e. you need
  `~donor + condition`, paired, or interaction terms).
- The replication is borderline (3 or fewer donors per condition) and you
  need to decide whether to proceed or fall back to description.

Methodology grounding: Squair et al. 2021 (Nat Commun, DOI
10.1038/s41467-021-25960-2) — pseudobulk methods that aggregate to the
sample level match bulk RNA-seq DE behavior on scRNA data and have
calibrated FPR; per-cell DE tools do not. Soneson & Robinson 2018 (Nat
Methods) is the earlier complementary benchmark.

## What `AggregateExpression` actually does

From the Seurat 5.5.0 source (verified via `formals(Seurat::AggregateExpression)`):

```r
formals(Seurat::AggregateExpression)
# $assays               NULL
# $features             NULL
# $return.seurat        FALSE
# $group.by             "ident"
# $add.ident            NULL
# $normalization.method "LogNormalize"
# $scale.factor         10000
# $margin               1
# $verbose              TRUE
```

For each unique combination of values in `group.by` (concatenated with `_`),
`AggregateExpression` **sums** the counts of all cells in that group, per
gene. The result is a counts matrix whose columns are the (donor × condition
× cell-type) pseudobulks. With `return.seurat = TRUE` the matrix is wrapped
in a Seurat object — column names = group keys joined by `_`, metadata
copied over from the per-cell `meta.data` (one row per group, with the
group's column value).

**Key semantics:**

- **Sum, not mean.** The pseudobulk counts are the SUM of cell counts in
  the group. DESeq2 / edgeR expect raw counts as input — they internally
  estimate size factors. If you average instead, the count semantics break
  and DESeq2 will silently misfit.
- **One assay at a time.** Pass `assays = "RNA"` to aggregate the RNA
  layer; ATAC/ADT need their own call (or a vector).
- **Counts layer is the source.** Even though `normalization.method =
  "LogNormalize"` appears in the formals, when `return.seurat = TRUE` the
  aggregated counts go into the `counts` layer of the returned object —
  that's what `FindMarkers(test.use = "DESeq2")` reads (via the `slot`
  argument routed under the hood). For native DESeq2, pull
  `GetAssayData(pseudo, assay = "RNA", layer = "counts")`.

## Donor / sample / cell-type structure — laying it out correctly

You need three things in `obj@meta.data` for Path C:

| Column | Role | Example values |
|---|---|---|
| Sample / donor identifier | The "replicate unit" — pseudobulks aggregate ALL cells from one donor × condition × cell-type | `donor_id`: `"d1"`, `"d2"`, `"d3"`, ... |
| Condition label | The variable you're testing the effect of | `stim`: `"STIM"` / `"CTRL"`; `disease`: `"PD"` / `"HC"` |
| Cell-type / cluster label | Stratifies the test (run one DE per cell-type) | `seurat_annotations`: `"CD14 Mono"`, `"CD8 T"`, ... |

The aggregation call then is:

```r
pseudo <- AggregateExpression(
  obj,
  assays      = "RNA",
  return.seurat = TRUE,
  group.by    = c("stim", "donor_id", "seurat_annotations")
)
```

The order in `group.by` controls the column-name order (`stim_donor_id_celltype`).
`paste(pseudo$seurat_annotations, pseudo$stim, sep = "_")` constructs the
`celltype.stim` ident the recipe uses for the comparison — choose the
ordering that makes the ident names readable.

### Replication requirements

A correct pseudobulk DESeq2 needs **≥3 donors per condition per cell-type**:

- Fewer than 3 → DESeq2's dispersion estimate is unstable; results are
  not robust DE — they're at best descriptive.
- Some cell-types in your object may have <3 donors per condition (rare
  populations, lost during one condition's QC). EXCLUDE those cell-types
  before running the test; print `table(pseudo$donor_id, pseudo$stim,
  pseudo$seurat_annotations)` to see the layout.

If your fixture lacks a real `donor_id` (e.g. `ifnb` ships with one
sample per condition, `orig.ident = "IMMUNE_CTRL" / "IMMUNE_STIM"`), the
Seurat v5 vignette mocks it for demonstration:

```r
set.seed(42)
obj$donor_id <- sample(c("d1","d2","d3"), size = ncol(obj), replace = TRUE)
```

This exercises the code path but **is NOT a biological claim** — the mocked
donors aren't real biological replicates. Note explicitly in any report.

## The `FindMarkers(test.use = "DESeq2")` wrapper — what it does and doesn't

Verified via `Seurat:::DESeq2DETest` (Seurat 5.5.0):

```r
function (data.use, cells.1, cells.2, verbose = TRUE, ...) {
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  dds1 <- DESeq2::DESeqDataSetFromMatrix(countData = data.use,
                                         colData = group.info,
                                         design = ~group)
  dds1 <- DESeq2::estimateSizeFactors(object = dds1)
  dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
  dds1 <- DESeq2::nbinomWaldTest(object = dds1)
  res <- DESeq2::results(dds1, contrast = c("group", "Group1", "Group2"), alpha = 0.05, ...)
  data.frame(p_val = res$pvalue, row.names = rownames(res))
}
```

What this means in practice:

- **Design is hardcoded `~group`.** Two-level factor: ident.1 vs ident.2.
  No room for additional covariates, no paired structure, no interaction.
- `fitType = "local"` for dispersion. Reasonable default; matches the
  Seurat vignette.
- The wrapper returns ONLY `p_val` — Seurat downstream adds `avg_log2FC`
  and `pct.1/pct.2` from the data matrix.
- `...` is passed to `DESeq2::results`, so you CAN tweak the contrast or
  `alpha` from the `FindMarkers` call, but not the design.

When `~group` is enough, this wrapper is the right call (one line, the
recipe's Path C). When it isn't, drop to native DESeq2 (next section).

## Multi-factor designs — drop to native DESeq2

For `~donor + condition`, paired designs, batch adjustment, interaction
terms, the Seurat wrapper can't express the model. Pull the pseudobulk
counts matrix and metadata out, build the `DESeqDataSet` yourself.

```r
# Step 1: get the pseudobulk counts and metadata.
counts <- as.matrix(GetAssayData(pseudo, assay = "RNA", layer = "counts"))
# Subset to one cell-type so the donor/condition design isn't mixed across cell-types.
keep_ct <- pseudo$seurat_annotations == "CD14 Mono"
counts  <- counts[, keep_ct]
coldata <- pseudo@meta.data[keep_ct, c("donor_id", "stim")]
stopifnot(all(rownames(coldata) == colnames(counts)))

# Step 2: drop genes with very low total counts to stabilize dispersion fit.
keep_g  <- rowSums(counts) >= 10
counts  <- counts[keep_g, ]

# Step 3: build the DESeqDataSet with the design you actually want.
# Example: condition effect adjusting for donor (paired-like for cross-condition donors).
library(DESeq2)
coldata$stim     <- factor(coldata$stim, levels = c("CTRL", "STIM"))    # CTRL = ref
coldata$donor_id <- factor(coldata$donor_id)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ donor_id + stim)
dds <- DESeq(dds)

# Step 4: extract the condition effect, log2FC shrunk for downstream interpretation.
res <- results(dds, name = "stim_STIM_vs_CTRL")
# Shrink with apeglm for visualization (DESeq2 vignette §"Log fold change shrinkage").
res_sh <- lfcShrink(dds, coef = "stim_STIM_vs_CTRL", type = "apeglm")

de_pb <- as.data.frame(res_sh)
de_pb$gene <- rownames(de_pb)
de_pb <- de_pb[order(de_pb$padj), ]
```

Design recipes for common multi-factor questions:

| Question | Design formula | Contrast |
|---|---|---|
| Stim effect adjusting for donor (cross-donor) | `~ donor_id + stim` | `stim_STIM_vs_CTRL` |
| Paired pre/post within donor | `~ donor_id + timepoint` (donor is the paired blocking factor) | `timepoint_post_vs_pre` |
| Stim × cell-type interaction | `~ stim * celltype` | `stimSTIM.celltypeCD14Mono` (the interaction term) |
| Stim adjusted for sex and batch | `~ sex + batch + stim` | `stim_STIM_vs_CTRL` |

DESeq2 vignette `bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html`
§"Multi-factor designs" is the authoritative reference; what's above is the
scRNA-pseudobulk-shaped subset.

## edgeR as an alternative — when DESeq2's dispersion fit fails

edgeR's quasi-likelihood F-test (`glmQLFit` + `glmQLFTest`) is the alternate
NB GLM workflow. Often interchangeable with DESeq2 results on the same
pseudobulk matrix; pick the one your downstream tooling expects.

```r
library(edgeR)
y <- DGEList(counts = counts, group = coldata$stim)
y <- calcNormFactors(y)                        # TMM normalization
design <- model.matrix(~ donor_id + stim, data = coldata)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = "stimSTIM")      # or whichever column
de_eR <- topTags(qlf, n = Inf)$table
de_eR$gene <- rownames(de_eR)
```

edgeR vignette `bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf`
§"Differential expression analyses" covers the full design space.

### When to pick edgeR over DESeq2 on pseudobulk

- Many cell-types × many comparisons — edgeR's `glmQLFit` is faster than
  `DESeq()` on small per-cell-type matrices.
- TMM normalization fits the data better than DESeq2's median-of-ratios
  for some compositional shifts (Squair 2021 SI uses both and reports
  similar calibration).
- Already invested in edgeR for upstream bulk analyses; consistency wins.

Both methods are pseudobulk-correct (Squair 2021). The choice is mostly
operational.

## Common pitfalls

| Pitfall | Symptom | Fix |
|---|---|---|
| Aggregated on `data` (log-normalized) instead of `counts` | DESeq2 errors `integer expected` or returns nonsense | Always `assays = "RNA"`; the counts layer is the source |
| Forgot to set the reference level on the condition factor | Contrast direction is reversed (STIM-UP appears as STIM-DOWN) | `coldata$stim <- factor(coldata$stim, levels = c("CTRL", "STIM"))` — first level is reference |
| Mixed cell-types in one `DESeqDataSet` | Dispersion fit dominated by cross-cell-type variance, not condition | Subset to one cell-type per `DESeqDataSet` (or use the interaction term) |
| <3 donors per condition × cell-type | DESeq2 errors `nrow(coldata) < nlevels(design factor)` or fits unstably | EXCLUDE that cell-type; don't fudge with mocked donors in real analyses |
| Used the per-cell counts matrix instead of pseudobulk | "Significant" hit count in thousands; FPR inflated | Aggregate first; see `anti_patterns.md` |
| Didn't shrink log2FC for visualization | Volcano plot dominated by low-count genes with huge unshrunk LFC | `lfcShrink(dds, coef = ..., type = "apeglm")` after `DESeq()` |

## Sources

- Seurat v5 DE vignette — satijalab.org/seurat/articles/de_vignette
- `formals(Seurat::AggregateExpression)`, `Seurat:::DESeq2DETest` (Seurat 5.5.0) — for function semantics
- DESeq2 vignette — bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
- edgeR users guide — bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
- Squair et al. 2021 — Confronting false discoveries in single-cell DE, Nat Commun (DOI 10.1038/s41467-021-25960-2)
- Soneson & Robinson 2018 — Nat Methods (DOI 10.1038/nmeth.4612)
- Love et al. 2014 — DESeq2 original paper, Genome Biology (DOI 10.1186/s13059-014-0550-8)
- Robinson, McCarthy & Smyth 2010 — edgeR, Bioinformatics (DOI 10.1093/bioinformatics/btp616)
