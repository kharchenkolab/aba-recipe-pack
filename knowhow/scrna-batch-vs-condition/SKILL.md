---
name: scrna-batch-vs-condition
description: Decision guide for distinguishing batch effect from real biological variation when scRNA-seq samples separate on the UMAP. Compares per-axis diagnostics (cluster composition by sample, marker conservation, QC-confound check, designed-control test, cell-type composition). Use when the user is staring at a UMAP that splits by sample and is choosing between "integrate it away" vs "this is the signal I want to study."
when_to_use: User has multi-sample scRNA-seq, the unintegrated UMAP shows cells from different samples landing in different regions, and they're asking "is this batch (should I integrate) or biology (should I preserve)?" Or asking what diagnostics to run before reaching for an integration recipe.
avoid_when: User already knows the separation is technical (e.g. same biological sample run on two days) and just wants an integration method — send them to `scrna-integration-decision` instead. Or there's only one sample (nothing to compare).
invocation: interactive
kind: knowhow_draft
requires_tools: [WebFetch, Read]
keywords: [batch effect, biological variation, integration, condition, UMAP, sample separation, diagnostics, LISI, kBET, marker conservation, cell type composition]
domain: genomics
source: "Tran et al. 2020 (Genome Biology, batch-correction benchmark, https://doi.org/10.1186/s13059-019-1850-9); Luecken et al. 2022 (Nat Methods, scIB atlas-integration benchmark, https://doi.org/10.1038/s41592-021-01336-8); Korsunsky et al. 2019 (Nat Methods, Harmony + LISI, https://doi.org/10.1038/s41592-019-0619-0); Büttner et al. 2019 (Nat Methods, kBET, https://doi.org/10.1038/s41592-018-0254-1); Lähnemann et al. 2020 (Genome Biology, eleven grand challenges, https://doi.org/10.1186/s13059-020-1926-6)."
audience: both
produces: []
capabilities_needed: []
---

# scRNA-seq — is sample separation batch or biology?

The question this knowhow answers: **"My UMAP shows cells from
different samples landing in different regions. Is that batch effect
(which I should correct via integration) or real biology (which I
should preserve)?"**

Audience: both naive (first multi-sample dataset, hasn't built
intuition for what batch looks like) and experienced (knows the
methods but wants a defensible decision before publishing). The
decision is fundamentally **case-dependent** — there is no global rule
that says "samples mixed on UMAP = batch removed correctly". Without
running the diagnostics in §4 you cannot tell the two regimes apart,
and choosing the wrong path either (a) removes the signal you wanted
to test or (b) leaves an obvious technical confound in the figure.
This knowhow gives you the diagnostic framework, not a single answer.

The executable recipes this knowhow points at when you've made the
call: `harmony-integration-scanpy`, `seurat-integration`,
`scvi-integration` (if integrating); `seurat-scrna` and
`scrna_qc_clustering` per sample (if NOT integrating). Diagnostic
metrics (LISI, kBET) currently have **no dedicated recipe** —
inline code is provided in §6.

## Quick decision

The framework is **five diagnostics, ANDed**. Run all five before you
decide. Do NOT make this call from the UMAP alone — UMAP-by-sample is
a hypothesis generator, not a verdict.

| Signal | Reading | What it points to |
|---|---|---|
| 1. **Cluster composition by sample** (table of cluster × sample) | Each cluster ~ proportional across samples | No integration needed |
| | Each cluster sample-pure | Either batch OR sample-specific biology — diagnostics 2-5 disambiguate |
| 2. **Marker conservation** within sample-pure clusters | Top markers per cluster MATCH across samples | Batch — integrate |
| | Top markers DIFFER across samples in matched clusters | Biology — do NOT integrate |
| 3. **QC confound check** | Sample-pure clusters track depth / %MT / %ribo / doublet score | "Biology" is QC artifact — re-QC, then revisit |
| | Sample-pure clusters do NOT track QC metrics | Real signal (batch or biology) — see 2/4/5 |
| 4. **Designed-control samples** | Biological replicates of the SAME condition cluster apart | Batch by definition — integrate |
| | Replicates of the same condition co-locate; different conditions separate | Biology — do NOT integrate over condition |
| 5. **Shared cell-type composition** | Known shared cell types (e.g. T cells in all PBMC samples) appear in only ONE sample | Batch or barcode/loading artifact — fix or integrate |
| | Truly sample-specific cell types (e.g. tumor cells absent from matched-normal) | Biology — analyze samples separately or use a conservative integration that preserves rare populations |

**Default ruling.** If diagnostic 4 is available (you have biological
replicates), it is the single most decisive signal — replicates of the
same condition that cluster apart prove batch by definition (Tran
2020). Without replicates, the call is irreducibly weaker; lean on
2+3+5 jointly.

**Position this knowhow takes.** When the diagnostics conflict (e.g.
markers match across samples but QC also tracks the split), assume
batch UNTIL you can rule out the technical explanation. The cost of
mistaking batch for biology (over-interpreting a technical
confound) is higher than the cost of mistaking biology for batch
(over-correcting, losing some signal) — because the second is at
least visible at sanity-check time (§6), while the first publishes.

## When NOT to integrate

These are the cases where sample-pure clustering IS the answer and
integration is the WRONG move.

- **The biological variable you want to test IS the sample variable.**
  If the question is "does condition X differ from condition Y" and
  each sample is one condition, then "samples separate on UMAP" is
  precisely the signal you came to measure. Integrating it away
  removes the signal (Lähnemann 2020 — batch and condition are
  unidentifiable when fully confounded). Integrate over the technical
  batch (lane, day, donor within condition) if you have it; do NOT
  integrate over the condition variable.
- **Cell types are sample-specific by design.** Tumor + matched
  normal, perturbed + unperturbed where the perturbation creates a
  new state, developmental time series where stages don't overlap.
  Aggressive integration (CCA, scVI without a `categorical_covariate`
  trick) can hallucinate the missing population in the wrong sample's
  manifold. Analyze samples separately + compare cluster markers, or
  use a conservative integration (RPCA / Harmony with high theta) and
  verify with diagnostic 5.
- **The "samples" are one biological sample re-sequenced.** Then
  there's no biological variable; everything between the runs IS
  batch. Use integration, but the answer is straightforward and this
  knowhow is overkill.
- **You haven't run per-sample QC.** Integration HIDES per-sample
  quality failures — one dying sample mixes with healthy ones and you
  never know. QC each sample with `scrna_qc_clustering` first.
  Integration after.
- **You haven't checked shared cell-type composition.** Aggressive
  integration on samples with very different cell-type composition
  produces apparently-mixed UMAPs that are actually over-corrected
  (Luecken 2022). Diagnostic 5 must be run before the decision, not
  after.

## Diagnostics — the framework in detail

The five axes from the Quick-decision table, with how to run each and
how to read the result.

### Diagnostic 1 — Cluster composition by sample

**Run.** Cluster the unintegrated data (default Leiden on PCA), then
tabulate cluster × sample.

```python
import pandas as pd
ct = pd.crosstab(adata.obs['leiden'], adata.obs['sample'], normalize='index')
# Heatmap; rows = clusters, columns = samples, values = fraction
```

**Read.**
- Each row ~ uniform across columns → samples already mix at the
  cluster level. Integration is unnecessary; the UMAP separation is
  visual artifact (UMAP placement is non-linear and can exaggerate
  small PCA-space separation; Lähnemann 2020 cautions against
  UMAP-as-truth).
- Each row dominated by one column → sample-pure clusters. This is
  the ambiguous regime — could be batch, could be biology. Go to
  diagnostics 2-5.

**Limitation.** Cluster resolution matters. At very high resolution
every cluster is sample-pure trivially. Use a moderate resolution
(0.5–1.0 leiden) and report the cluster count.

### Diagnostic 2 — Marker conservation

**Run.** For each cluster, compute top markers (`sc.tl.rank_genes_groups`
with Wilcoxon, or Seurat `FindMarkers`) PER SAMPLE — i.e. within each
sample, contrast that sample's cells in the cluster vs that sample's
other cells. Compare the per-sample marker lists for the SAME cluster.

**Read.**
- Top-20 marker overlap across samples is high (Jaccard > ~0.5) AND
  the same canonical markers (e.g. CD3D/CD8A for T cells) appear in
  both → the cluster is the SAME biological cell type in both
  samples; the sample-pure clustering is BATCH on top of shared
  biology → integrate.
- Top markers differ AND no canonical shared markers → the clusters
  are DIFFERENT cell types in different samples → biology → do NOT
  integrate over the sample variable.

**Limitation.** Requires enough cells per (sample × cluster) for
markers to stabilize (rule of thumb: ≥50 cells). Sparse clusters give
noisy comparisons.

### Diagnostic 3 — QC confound check

**Run.** Color the unintegrated UMAP by `n_counts` (depth), `pct_mt`,
`pct_ribo`, and the doublet score. Or correlate cluster identity with
these per-cell QC metrics.

```python
import scanpy as sc
sc.pl.umap(adata, color=['sample', 'n_counts', 'pct_mt', 'pct_ribo', 'doublet_score'])
```

**Read.**
- Sample-pure clusters track ANY QC metric (high-MT samples cluster
  together, low-depth samples cluster together) → the "biology" you
  see is a QC artifact. Re-do per-sample QC with stricter thresholds
  (`scrna_qc_clustering`), then revisit.
- No tracking → the separation is real signal (batch OR biology;
  diagnostics 2, 4, 5 distinguish).

**Limitation.** Depth differences alone don't always mean QC failure
— some real cell types (e.g. small lymphocytes) genuinely have lower
counts. Use the per-sample QC distribution (not absolute thresholds)
to judge.

### Diagnostic 4 — Designed-control samples

**Run.** If you have biological replicates of the SAME condition
(e.g. control_1 and control_2 from different donors), tabulate
cluster × replicate within condition.

**Read.**
- Replicates of the same condition cluster APART → batch effect by
  definition (the only variable between them is technical). Integrate
  the technical covariate (`donor` / `lane` / `run`). This is the
  strongest single signal in the framework.
- Replicates co-locate; different conditions separate → the
  separation IS the biological signal. Do NOT integrate over
  condition. (Tran 2020 builds the entire benchmark logic on this:
  if shared cell types from different batches don't mix, batch
  correction is required.)
- Replicates partially co-locate (some clusters mix, some don't) →
  mixed regime — integrate the technical batch covariate (which
  Harmony, scVI handle) and then verify with diagnostic 6.

**Limitation.** Many real datasets don't have replicates within
condition. This is the irreducible weakness of single-replicate
designs (Lähnemann 2020 — batch identifiability requires replication).
Without replicates the decision relies on 2/3/5 jointly.

### Diagnostic 5 — Shared cell-type composition

**Run.** For canonical cell-type markers expected in your tissue
(e.g. T cells via CD3D, B cells via MS4A1, monocytes via CD14 in
PBMC), score each sample for the presence of each cell type.

```python
sc.tl.score_genes(adata, gene_list=['CD3D', 'CD3E', 'CD8A'], score_name='Tcell_score')
sc.pl.violin(adata, 'Tcell_score', groupby='sample')
```

**Read.**
- Canonical shared cell types appear in ALL samples (the violin
  distributions look similar) but the cells nonetheless cluster
  apart on UMAP → batch (the same cells, technically separated).
  Integrate.
- A canonical shared cell type is missing in one sample entirely
  (the violin is flat for that sample) → either (a) loading /
  barcode / capture problem (technical, worth investigating before
  any analysis), or (b) genuinely sample-specific biology (e.g. a
  cell type ablated by the perturbation). The diagnosis depends on
  whether the absence is expected.
- A genuinely sample-specific cell type that's expected (tumor
  cells only in the tumor sample) → biology — analyze separately
  or use conservative integration that explicitly preserves rare /
  unique populations.

**Limitation.** Requires knowing the expected cell-type composition
of the tissue. For unknown / exploratory tissues, fall back on
diagnostics 1-4. # REVIEW(scrna-biology): is there a tissue-agnostic
proxy for "expected shared cell-type composition" beyond known
marker lists? Co-expression module preservation across samples?

## Anti-patterns

- **"Samples mix on the UMAP, so integration worked."** UMAP can
  hide misalignment (Luecken 2022 — both bio-conservation and
  batch-removal metrics are needed; visual inspection alone is
  insufficient). Integration that mixes samples WHILE collapsing
  distinct cell types (over-correction) ALSO produces a mixed UMAP.
  Verify in PCA space + with LISI/kBET (§6), not just UMAP-by-sample.
- **Integrating over the biological variable.** If the question is
  "stim vs ctrl" and each sample is one condition, integrating with
  `sample` (or `condition`) as the batch key removes the signal you
  came to measure. Use the technical covariate (donor, lane, day)
  instead. If technical and biological are fully confounded
  (one-donor-per-condition), the experiment is unidentifiable
  (Lähnemann 2020) — say so, don't integrate-and-pretend.
- **Reaching for integration as the first move.** Sample-pure
  clustering is a symptom; integration is one of several treatments.
  Diagnose first (the five axes), THEN pick a treatment (integrate,
  re-QC, analyze separately, or accept the biology).
- **Skipping diagnostic 5 (shared composition) before integration.**
  If samples have radically different cell-type composition,
  integration over-corrects and hallucinates shared populations.
  Verify the populations you expect to share actually exist in each
  sample first.
- **Choosing integration method by paper recency.** "scANVI is
  newer" is not a benchmark. Use Luecken 2022 / Tran 2020 to pick;
  default to Harmony for first-pass and verify (Korsunsky 2019,
  Tran 2020 — Harmony top-performing across non-identical-cell-type
  benchmarks).
- **Reporting integrated counts.** Integration corrects the
  embedding (PCA / UMAP / latent), not the count matrix. DE on the
  raw counts with sample as a covariate (per-cell) or pseudobulk
  (per-sample). See `scrna-de-methodology` knowhow when it lands.
- **Treating "sample-pure clustering = biology" without checking
  QC.** A dying sample makes its own cluster (high MT, low depth).
  That cluster is QC artifact, not biology. Run diagnostic 3
  first.

## Sanity checks — verifying the call

After you've decided (integrate vs not), run these checks.

### If you DECIDED TO INTEGRATE

- **UMAP by sample after integration.** Samples should be mixed
  within clusters that DO share biology. Some sample-purity may
  remain if some cell types are genuinely sample-specific — that's
  fine, IF it matches what diagnostic 5 predicted.
- **UMAP by cell type after integration.** Same cell type from
  different samples should co-locate. If T cells from sample A and
  T cells from sample B are in different clusters post-integration,
  the integration failed.
- **Marker conservation post-integration.** Re-run marker analysis
  on the integrated clusters; top markers should match per-sample
  marker analysis (diagnostic 2). Heavy marker drift indicates the
  integration is altering biology — switch to a more conservative
  method (RPCA, Harmony with higher `theta`).
- **iLISI / cLISI** (Korsunsky 2019 — Local Inverse Simpson Index).
  iLISI on the sample/batch label should approach the number of
  batches (good mixing); cLISI on the cell-type label should stay
  near 1 (preserved biology). Both, not either — Luecken 2022
  emphasizes the trade-off.
- **kBET** (Büttner 2019 — k-nearest-neighbour batch effect test).
  Reports a rejection rate per neighborhood; low post-integration
  = well-mixed. Use on the integrated PCA / latent embedding, not
  the UMAP.

### If you DECIDED NOT TO INTEGRATE

- **Per-sample cluster-marker comparison.** Run
  `seurat-scrna` (or equivalent) per sample, then manually match
  clusters by canonical markers across samples. The clusters that
  match are shared cell types; the clusters that don't are
  sample-specific.
- **Don't compare DE results across separately-analyzed samples
  cell by cell.** Per-cell DE across unaligned objects is
  ill-defined. Aggregate to pseudobulk (per sample × cluster) and
  compare; or do per-sample DE and compare the LISTS of hits, not
  the per-cell statistics.
- **Document the decision.** If the paper reviewers ask why you
  didn't integrate, the answer is the diagnostics — write down
  which axes told you not to.

### Inline LISI / kBET code (no recipe yet)

LISI via `harmonypy`:

```python
import harmonypy as hm
import numpy as np
# adata.obsm['X_pca'] from sc.tl.pca, adata.obs['sample'] the batch label
lisi = hm.compute_lisi(adata.obsm['X_pca'], adata.obs[['sample', 'celltype']], ['sample', 'celltype'])
print('iLISI (sample, want high, max = n_samples):', np.median(lisi[:, 0]))
print('cLISI (celltype, want low, min = 1):       ', np.median(lisi[:, 1]))
```

kBET via the R package (`theislab/kBET`):

```r
library(kBET)
# pca <- prcomp(...)$x  (or output of any embedding)
# batch <- factor(metadata$sample)
res <- kBET(pca, batch)
# res$summary$kBET.observed — rejection rate; low = well-mixed
```

# TODO(recipe): standalone `scrna-batch-diagnostics` recipe that
runs all five diagnostics + LISI + kBET in one shot. Currently the
user has to copy these snippets.

## See also

**Benchmarks + reviews (cited):**
- Tran et al. 2020, Genome Biology — batch-correction benchmark
  across 14 methods; the canonical "should I integrate" empirical
  reference (https://doi.org/10.1186/s13059-019-1850-9)
- Luecken et al. 2022, Nat Methods — scIB atlas-scale integration
  benchmark; bio-conservation vs batch-removal trade-off framework
  (https://doi.org/10.1038/s41592-021-01336-8)
- Lähnemann et al. 2020, Genome Biology — eleven grand challenges
  in single-cell data science; identifiability of batch vs
  biology under confounding
  (https://doi.org/10.1186/s13059-020-1926-6)

**Method + metric papers:**
- Korsunsky et al. 2019, Nat Methods — Harmony method, introduces
  LISI for batch mixing + cell-type preservation
  (https://doi.org/10.1038/s41592-019-0619-0)
- Büttner et al. 2019, Nat Methods — kBET k-nearest-neighbour
  batch-effect test
  (https://doi.org/10.1038/s41592-018-0254-1)

**Recipes that execute paths from §4-§6:**
- `harmony-integration-scanpy` — Harmony via scanpy (first-pass
  integration default)
- `seurat-integration` — Seurat v5 IntegrateLayers (CCA / RPCA /
  Harmony / scVI options)
- `scvi-integration` — scVI standalone (atlas-scale)
- `seurat-scrna` — per-sample Seurat workflow for the
  "do NOT integrate, analyze separately" path. # REVIEW(catalogue):
  the assignment referenced `seurat-scrna-v2` but the catalogue
  has `seurat-scrna` — confirm which recipe should be the canonical
  per-sample pointer.
- `scrna_qc_clustering` — per-sample QC + clustering (runs before
  any of the above)

**Adjacent knowhow:**
- `scrna_pipeline.md` — upstream framing: multiple samples → keep
  separate by default, do NOT naively concatenate
- `scrna-integration-decision` — sibling knowhow: WHICH integration
  method (Harmony vs CCA vs scVI vs conos) once you've decided
  integration is the answer. # TODO(knowhow): this sibling does not
  yet exist in the catalogue at draft time — placeholder cross-ref.
- `scrna-de-methodology` — downstream: DE on the integrated /
  per-sample-analyzed objects (use pseudobulk for condition
  contrasts; per-cell Wilcoxon only for cluster markers within
  sample). Currently being drafted in parallel.

**Coverage gaps (see REVIEW_LOG.md):**
- No standalone diagnostic recipe — LISI / kBET / cross-tabulation
  given as inline snippets only.
- The sibling `scrna-integration-decision` knowhow is referenced but
  not yet authored.
- Reference to `seurat-scrna-v2` in the original framing is unresolved
  against the actual catalogue entry `seurat-scrna`.
