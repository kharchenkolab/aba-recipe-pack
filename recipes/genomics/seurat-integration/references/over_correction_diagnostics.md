# Diagnosing under- and over-integration

When to load this: the post-integration UMAP doesn't look right — samples
still form sample-specific islands (under-integration), or biologically
distinct cell types collapsed into one cluster (over-integration); the user
asks "is this enough integration?"; you need to justify a method switch.

The SKILL.md Step 6 produces the two diagnostic UMAPs and the
sample-by-cluster contingency table. This reference explains how to read
them, what to do when they look wrong, and how to confirm a fix.

## What "integrated correctly" actually looks like

Three conditions, ALL needed:

1. **By-sample UMAP is well-mixed.** Within each visible cluster, the
   points from each sample interleave — no sample-only sub-island per
   cell type.
2. **By-cluster UMAP shows distinct populations.** Clusters look like
   coherent blobs with clear boundaries, not one giant smear, not so
   many that adjacent ones overlap heavily.
3. **Sample × cluster contingency is balanced where it should be.** For
   a balanced cohort (similar cell-type compositions across samples),
   no single cluster should be >80% from one sample. For a cohort
   where the biology genuinely differs (e.g. one disease arm has a
   rare cell type), per-sample skew at THAT cluster is real biology.

```r
# The diagnostic
sample_cluster_tab <- table(cluster = obj$integrated_clusters,
                            sample  = obj$sample)
print(sample_cluster_tab)

# Per-cluster sample fractions
print(round(prop.table(sample_cluster_tab, margin = 1), 2))

# Largest single-sample share per cluster
apply(prop.table(sample_cluster_tab, margin = 1), 1, max)
```

A cluster with `max share > 0.8` in a balanced cohort either is a
sample-specific population (look up the markers — does it match a
known cell type that should be uniformly present?) or the integration
under-corrected at that cluster.

## Symptom 1 — Under-integration

**Signs:**
- By-sample UMAP shows ctrl/stim (or per-sample colors) forming
  separate sub-clusters within each cell-type island.
- Sample × cluster contingency: many clusters dominated by one sample.
- Cluster markers include cell-type markers *plus* sample-specific
  signals (e.g. ISG / IFN response genes co-driving a cluster split).

**Root-cause options:**
- Real per-sample biology — the samples really do differ (e.g. stim vs
  ctrl). If the user wants to test that, you should NOT integrate over
  the variable that carries it. Use the variable as a covariate in
  downstream DE instead.
- Method too conservative — Harmony or RPCA may leave residual batch
  effect on strongly-batched data.
- Wrong number of dims — too few dims, the integration can't span
  enough variation to align.
- Too few HVGs — `FindVariableFeatures(nfeatures = 2000)` may miss
  cross-sample-conserved features.

**Patches in order of effort:**

1. **Raise `dims` to 1:40** in both PCA and IntegrateLayers — quick
   sanity check.
2. **Switch method to CCA** (more aggressive) from RPCA/Harmony, or
   to RPCA from Harmony.
3. **Raise `k.anchor`** in the CCA/RPCA call (passes through to
   `FindIntegrationAnchors`):
   ```r
   obj <- IntegrateLayers(
     object = obj, method = CCAIntegration,
     orig.reduction = "pca", new.reduction = "integrated.cca",
     k.anchor = 20,                  # default 5 → 20 to strengthen
     verbose = FALSE
   )
   ```
4. **Add more HVGs** — `FindVariableFeatures(nfeatures = 3000)`.
5. **Confirm signal-vs-nuisance** — make sure the split variable is
   the nuisance covariate (sample / batch), not the biological one
   you wanted to test.

After patching, re-run Steps 4–6 (Integrate → cluster → assess) and
compare the new by-sample UMAP + contingency table.

## Symptom 2 — Over-integration

**Signs:**
- By-cluster UMAP collapses biologically distinct cell types into one
  cluster — e.g. T cells and NK cells share a cluster despite distinct
  marker profiles.
- Some clusters that the SINGLE-sample analysis (run via
  `seurat-scrna-v2`) had clearly separated are now merged.
- Cluster markers are weak — top genes per cluster are not specific.
- "Rare" cell types disappear — e.g. dendritic cells (DCs) were ~1%
  per sample, post-integration they fold into monocytes.

**Root-cause options:**
- Method too aggressive — CCA on small cohorts can pull rare
  populations into the larger one.
- Too few HVGs — the integration's anchor finder doesn't have enough
  features to distinguish rare types.
- The split variable is the wrong covariate — integrating over a
  variable that carries biology erases that biology.

**Patches in order of effort:**

1. **Switch CCA → RPCA** — more conservative, preserves more
   sample-specific signal.
2. **Lower `k.anchor`** — fewer anchors per query cell limits how
   much the method pulls populations together.
3. **Lower `k.weight`** in `IntegrateLayers(... k.weight = 50)` —
   tighter local anchor weighting.
4. **More HVGs** — `nfeatures = 3000`.
5. **Reduce `dims`** — too many dims can amplify a weak nonlinear
   correction.
6. **Confirm split variable** — if integrating over `donor` accidentally
   over-correlates with `condition` (e.g. donor 1 = ctrl, donor 2 = stim),
   the correction will erase the condition signal too.

## Marker preservation as a sanity check

After integration, run `FindAllMarkers` (post-`JoinLayers`) and
inspect the top markers per cluster. Compare against either:

- the per-sample marker lists from a single-sample run, or
- a known cell-type marker panel for the tissue.

Markers that survive integration are robust; markers that
**disappear** post-integration may have been signal that was
absorbed by the correction. Markers that **emerge** post-integration
that weren't present per-sample are often integration artifacts.

```r
# After integration + JoinLayers
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25)
top_per <- markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(avg_log2FC, n = 5)
print(top_per)
```

For tissue-typical markers (PBMC, etc.), see the
`annotate-celltype-scrna` recipe.

## Method comparison — same object, side by side

The cleanest "did integration help?" answer is to run two methods on the
SAME object and compare. See
`references/integration_methods.md` §"Comparing two methods on the same
object" for the exact pattern.

Visual inspection criteria for the comparison:

| Method A vs B | Read it as |
|---|---|
| A's by-sample UMAP mixes better than B's | A is the right pick for THIS dataset |
| A collapses a population B keeps separate | B is the right pick — A is over-integrating |
| Both look identical | Effect is small either way; pick the faster one (RPCA > CCA, Harmony > both) |

## What NOT to do

- **Don't iterate `k.anchor` upward until the by-sample UMAP looks
  perfect.** A fully-mixed by-sample UMAP at the cost of collapsed
  biology is worse than residual batch. Look at the by-cluster UMAP +
  markers, not just the by-sample one.
- **Don't switch methods after every iteration.** Pick one method per
  analysis; document the choice and parameters.
- **Don't integrate over a variable that carries the biology you want
  to test.** It's not just over-correction, it's erasing the signal
  before you even ask the question.
