# SCT DE gotchas — `PrepSCTFindMarkers`, SCT vs RNA assay, and slot semantics

The specific gotchas around running differential expression on the SCT
assay: the mandatory `PrepSCTFindMarkers` gate, the SCT-vs-RNA choice, why
`slot = "data"` is right (and `"scale.data"` is wrong) for Wilcoxon, and
what changes when objects are merged.

Load this reference when:
- A `FindAllMarkers` / `FindMarkers` call on an SCT assay errored.
- The user asks "should I do DE on the SCT or RNA assay?"
- Multi-sample SCT and you're seeing degenerate marker results.
- A merged SCT object's DE produces nonsensical fold-changes.

## The `PrepSCTFindMarkers` gate

Source-verified signature (Seurat 5.5.0):

```r
PrepSCTFindMarkers(object, assay = "SCT", verbose = TRUE)
```

**Rule:** call `PrepSCTFindMarkers(obj)` BEFORE `FindAllMarkers` or
`FindMarkers` on an SCT assay. Always. Single-sample, multi-sample,
integrated, doesn't matter — Seurat errors otherwise.

```r
# CORRECT
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)
markers <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE,
                          min.pct = 0.25, logfc.threshold = 0.25,
                          verbose = FALSE)

# WRONG — Seurat will error with a "PrepSCTFindMarkers must be called" message
markers <- FindAllMarkers(obj, assay = "SCT", ...)
```

### What `PrepSCTFindMarkers` actually does

For a MERGED multi-sample SCT object, each sample's SCT model has its own
`median_umi` (the median UMI count across cells used for the per-model
count correction). When you merge two SCT objects, the corrected counts
between them are on DIFFERENT scales (because the medians differ).
`PrepSCTFindMarkers` re-corrects all samples to the MINIMUM `median_umi`
across models, putting them on a shared scale.

For a SINGLE-sample SCT object, there's only one model, so the
re-correction to "the minimum across one model" is a no-op — the function
runs and exits cleanly without changing data. The gate is still required
because Seurat's DE code checks for the `PrepSCTFindMarkers`-set marker;
the no-op call satisfies that check.

So:

| Object | Effect of `PrepSCTFindMarkers` |
|---|---|
| Single-sample SCT | No data change, but the call is REQUIRED (Seurat checks an internal flag) |
| Merged multi-sample SCT | Real data change — counts re-corrected to shared median |
| RNA assay only (no SCT) | Don't call this; only applies to SCT |

### What gets re-corrected

`PrepSCTFindMarkers` modifies `obj[["SCT"]]@counts` (and `@data` =
log1p(counts) downstream). The `scale.data` (residuals) is unaffected
because residuals are already on a shared scale (Pearson residuals have
~unit variance by construction).

This means:
- DE using `slot = "data"` (the default) — uses the re-corrected counts.
  This is what `FindAllMarkers` reads.
- DE using `slot = "counts"` — same; same re-correction.
- DE using `slot = "scale.data"` (residuals) — would be unchanged, but
  `FindAllMarkers` doesn't use `scale.data` for Wilcoxon (and shouldn't —
  see below).

## SCT vs RNA assay for DE — which to pick

There are TWO reasonable defaults; pick ONE per analysis and report which.

### SCT assay DE (the recipe's default)

```r
obj <- PrepSCTFindMarkers(obj, assay = "SCT", verbose = FALSE)
markers <- FindAllMarkers(obj, assay = "SCT", only.pos = TRUE, ...)
```

**Pros:**
- Uses the depth-corrected counts. Markers are not biased by per-cell
  sequencing depth.
- The same data clustering saw — internal consistency.
- Across samples (merged SCT object), the re-correction puts cells on a
  shared scale.

**Cons:**
- Fold-changes are on the corrected count scale, not raw counts. Effect
  sizes are slightly different from what a bulk-RNA paper would report.
- For very lowly expressed genes the regularization shrinks effects toward
  the mean — slight loss of sensitivity on weak markers.

### RNA assay DE (the alternative)

```r
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj, verbose = FALSE)          # RNA assay only — SCT doesn't normalize RNA
markers_rna <- FindAllMarkers(obj, assay = "RNA",
                              only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.25, verbose = FALSE)
DefaultAssay(obj) <- "SCT"                          # restore
```

**Pros:**
- Marker effect sizes match what bulk-RNA / older single-cell papers
  report — easier to cross-reference.
- No depth correction means raw biological signal is preserved.
- Required if you want to mix SCT-clustered identities with RNA-assay DE
  output (e.g. for downstream pseudobulk DE which expects raw counts).

**Cons:**
- DE results are depth-confounded. If clusters differ in mean UMI count,
  fold-changes between them inflate.
- Not internally consistent with the SCT clustering.

### What NOT to do

- ❌ **Run DE without `PrepSCTFindMarkers` on the SCT assay.** Seurat
  errors.
- ❌ **Mix SCT and RNA DE in the same figure.** Pick one. Mixing misleads
  the reader about effect sizes.
- ❌ **`FindAllMarkers(slot = "scale.data")`.** Pearson residuals are NOT
  what Wilcoxon expects — they can be negative, and the `min.pct` /
  `only.pos` semantics break. `slot = "data"` (default) is correct for
  Wilcoxon on SCT.
- ❌ **`NormalizeData` on the SCT assay.** SCT already normalized. Running
  `NormalizeData` overwrites SCT data with a log-normalized version of the
  CORRECTED counts, which is incorrect and breaks the model.

## `test.use` choices on SCT — same as RNA except `negbinom`/`poisson`

The full `test.use` enum (Wilcoxon, MAST, ROC, t, LR, …) applies on SCT
the same as on RNA, EXCEPT:

- **`negbinom`** and **`poisson`** require raw UMI counts (`slot = "counts"`).
  On SCT, the `counts` slot holds CORRECTED counts (depth-adjusted), not
  raw. Using `negbinom` on SCT counts violates the model assumptions.
  For these tests, switch to `DefaultAssay(obj) <- "RNA"` and use the
  RNA assay's `counts` layer.

For other tests on SCT, `assay = "SCT", slot = "data"` is correct.

## Multi-sample SCT — when merging matters

The single-sample recipe doesn't exercise the merge path, but the gotchas
are worth knowing:

### Sample-level SCT models survive merge, but with caveats

After `merge(obj1, obj2)`:
- `obj_merged[["SCT"]]@SCTModel.list` has length 2.
- Each cell knows which model it came from (`obj@active.ident` or
  `obj$orig.ident`).
- `obj_merged[["SCT"]]@counts` is the concatenation of the two corrected
  count matrices — but they were corrected to DIFFERENT `median_umi`s.
- `obj_merged[["SCT"]]@scale.data` is the concatenation of residuals from
  the two independent fits. Each cell's residuals are correct relative to
  its own model.

`PrepSCTFindMarkers(obj_merged)` re-corrects the counts (and data) to the
shared minimum median UMI. THEN `FindAllMarkers` works.

### When to integrate vs. merge

Merging + DE on SCT works only if the samples are biologically replicates
(same tissue, same condition). For samples that differ in condition / batch
/ donor, use proper integration:

1. SCTransform per sample.
2. `SelectIntegrationFeatures` + `PrepSCTIntegration`.
3. `FindIntegrationAnchors(normalization.method = "SCT", ...)` /
   `IntegrateData(normalization.method = "SCT", ...)` (v4 API), or
4. `IntegrateLayers(method = CCAIntegration, normalization.method = "SCT", ...)` 
   (v5 API).

Then `PrepSCTFindMarkers` + DE.

Integration is OUT OF SCOPE for this recipe — see `seurat-integration`
for the full workflow.

## Pseudobulk DE for SCT-clustered samples

When you have multi-sample SCT data and want population-level DE (case vs
control across samples), per-cell DE is the WRONG approach (false-positive
inflation from pseudo-replication). Pseudobulk DE:

1. SCTransform each sample, cluster, identify cell types.
2. Aggregate per (cell type, sample) into a pseudo-bulk count matrix
   (`AggregateExpression(obj, group.by = c("cell_type", "sample"))` or 
   `Seurat::PseudobulkExpression`).
3. Run DESeq2 / limma / edgeR on the pseudo-bulk matrix.

The `seurat-de-testing` recipe documents the pseudobulk path (Path C in
that recipe).

## DE on a non-variable gene

If you need to test a specific gene that wasn't in `VariableFeatures(SCT)`
(because `return.only.var.genes = TRUE` is the default), the gene's
residuals weren't stored in `scale.data`. Two options:

1. **Re-run SCTransform with `return.only.var.genes = FALSE`.** Doubles
   memory but stores residuals for every gene. Best for downstream
   exploration.
2. **Use the RNA assay for that gene's DE.** `FindMarkers(obj, features =
   "GENE_X", assay = "RNA", ...)`. Doesn't need re-running SCT but mixes
   assays — be explicit about which you're showing.

For `FeaturePlot` / `DotPlot` on a non-variable gene, the SCT `data` layer
(log1p of corrected counts) is computed for ALL genes regardless of
`return.only.var.genes` (only `scale.data` is restricted). So plotting a
non-variable gene on the SCT assay works without changes.

## Quick reference — when in doubt

| Question | Answer |
|---|---|
| Do I always need `PrepSCTFindMarkers`? | YES, before any FindAllMarkers / FindMarkers on the SCT assay |
| Can I skip it for single-sample SCT? | NO — Seurat checks an internal flag, errors without |
| Should I use SCT or RNA for DE? | SCT for default; RNA for raw-count effect sizes or `negbinom`/`poisson` |
| What slot for Wilcoxon? | `slot = "data"` (default). Don't override |
| What slot for `negbinom`/`poisson`? | `slot = "counts"`, AND switch to `assay = "RNA"` (SCT counts are corrected, not raw) |
| Can I mix SCT and RNA DE in one analysis? | Pick ONE per figure. Mixing misleads on effect sizes |
| What about pseudobulk DE? | Use the `seurat-de-testing` recipe — DESeq2 on AggregateExpression output, NOT FindAllMarkers |
