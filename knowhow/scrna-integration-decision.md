---
name: scrna-integration-decision
description: Decision guide for choosing a scRNA-seq sample/batch integration method. Compares Harmony, Seurat CCA, Seurat RPCA, scVI, conos, and the BBKNN/Scanorama/FastMNN tier. Use when the user has multiple scRNA-seq samples whose cells separate by sample in the unintegrated UMAP and is choosing between methods.
when_to_use: User has 2+ scRNA-seq samples, the unintegrated PCA/UMAP shows cluster-by-sample rather than cluster-by-cell-type, and they are asking "which integration method?" or "Harmony vs CCA vs scVI for my data?". Also use when an agent is binding an integration recipe and needs a defensible reason to prefer one over the others.
avoid_when: User has one sample (no integration to do — see scrna-qc-clustering), or the biological question IS the sample/condition contrast (integrating over it removes the signal — see anti-patterns §5), or the user wants a runnable workflow (route to the executable recipe directly: harmony-integration / harmony-integration-scanpy / seurat-integration / scvi-integration / conos-integration / bp-data-integration).
invocation: interactive
kind: knowhow_draft
requires_tools: [WebFetch, Read]
keywords: [integration, batch correction, Harmony, CCA, RPCA, scVI, scANVI, conos, BBKNN, Scanorama, FastMNN, LIGER, scRNA-seq, sample effect, donor effect, batch key, decision guide]
domain: genomics
source: "Luecken et al. 2022 Nat Methods benchmark; Korsunsky 2019 (Harmony); Stuart 2019 (Seurat v3 CCA); Lopez 2018 (scVI); Polanski 2020 (BBKNN); Hie 2019 (Scanorama); Heumos et al. 2023 single-cell best practices; ABA executable recipes in aba-recipes/recipes/genomics/."
audience: both
---

# scRNA-seq integration — which method for which dataset?

The question this knowhow answers: **"I have multiple scRNA-seq samples and the
cells from different samples cluster apart in the unintegrated UMAP — which
integration method should I use?"**

Audience: both naive (first multi-sample analysis) and experienced
(disambiguating between methods for a specific dataset shape). The decision is
multi-axis (cell count × language stack × downstream goal × cell-type overlap
across samples × whether cell labels exist). The executable recipes this
knowhow connects to are `harmony-integration` (R/Seurat),
`harmony-integration-scanpy` (Python), `seurat-integration` (R/Seurat v5
IntegrateLayers — CCA/RPCA/Harmony), `scvi-integration` (Python), and
`conos-integration` (R), plus `bp-data-integration` (the scIB best-practice
multi-method Python flow that covers BBKNN/Scanorama/scANVI alongside).

## Quick decision

For the typical multi-sample dataset (2–10 samples, ≤200k cells, shared
cell-type composition, no cross-species/cross-modality complication):

| Your situation | Use | Recipe |
|---|---|---|
| Default — fast, conservative, well-validated, **you just want it to work** | **Harmony** | `harmony-integration-scanpy` (Python) / `harmony-integration` (R) |
| You're already in R/Seurat v5 | **Harmony via IntegrateLayers** OR **RPCA** (more conservative than CCA) | `seurat-integration` |
| Atlas-scale (>500k cells), complex/nested batch effects, OR you'll do downstream scANVI label transfer / scArches mapping | **scVI** (or **scANVI** if you have labels) | `scvi-integration` |
| You have trusted cell-type labels and want the integration to preserve them | **scANVI** (label-aware variant of scVI) | `scvi-integration` then scANVI step (see `bp-data-integration`) |
| You want a multi-method scIB-ranked comparison before committing | **Run several; rank with scIB** | `bp-data-integration` |
| You want a joint kNN graph without forcing a shared embedding (e.g. trajectory work where the embedding is fragile) | **conos** | `conos-integration` |

If none of these fit (cross-species, very low depth, rare-population focus,
condition-imbalanced cell types), read §4 before picking. Harmony first +
sanity-check the result (§6) is rarely the wrong default for a first pass —
the literature consistently lists it among top performers on "simple"
integration tasks (Luecken 2022 [1]; sc-best-practices Heumos et al. 2023 [8]).

## When NOT to integrate

- **One sample.** Integration is correction for sample-effect; one sample has
  no sample-effect to correct. Run `scrna-qc-clustering` (Python) or
  `seurat-scrna` (R) on the sample directly. (See the existing knowhow
  `scrna_pipeline.md` for the same guardrail.)
- **The biological question IS the sample/condition variable.** If you are
  asking "does stim differ from ctrl", integrating *over* `condition` removes
  the very signal you want to test. Integrate over the technical nuisance
  variable (donor / lane / 10x run / day), not the biological one. Per-cell
  DE on the integrated embedding is also wrong here — use pseudobulk DESeq2
  on the raw counts (Squair 2021; see ABA memory item `scrna_de_method`).
- **Cell types are sample-specific by design** (e.g. tumor-only vs
  matched-normal-only populations, or a cell line absent from the donor
  sample). Aggressive integration over-pulls these into the wrong manifold.
  Use a conservative method (Harmony with higher `theta`, or RPCA — see
  §4) — or analyze each sample separately and *compare* clusters.
  # REVIEW(integration-methods): the "Harmony with high theta" guidance is
  # community lore for conservative integration; is there a specific
  # benchmark or Korsunsky-lab recommendation we should cite for the
  # theta-as-conservatism-knob claim?
- **Per-sample QC hasn't been done.** Integration hides per-sample quality
  failures (a low-quality sample looks like a "batch effect" the method
  faithfully erases). Run `scrna-qc-clustering` per sample first; integrate
  the survivors.
- **The unintegrated UMAP already mixes samples within clusters.** If
  cells from different samples co-locate by cell type without correction,
  there is no batch effect worth correcting — adding integration only
  risks over-correction and bookkeeping cost (Heumos et al. 2023 [8] —
  "visualize first; you may not need it"). Skip integration; proceed to
  clustering / annotation on the raw PCA.

## Alternatives — full matrix

| Method | Characterization | When it wins | When it fails | Cost | Recipe |
|---|---|---|---|---|---|
| **Harmony** (Korsunsky 2019 [3]) | Iterative soft-clustering in PCA space; learns dataset correction factors per soft cluster, leaves counts untouched, writes a corrected embedding | Most "simple-batch" datasets; fast; conservative; well-validated by independent benchmarks (Luecken 2022 [1]; sc-best-practices [8]) | Highly imbalanced cell-type composition across samples; when correction has to be very large (Luecken 2022 [1] ranks scANVI/scVI/Scanorama above Harmony on "complex" tasks) | Linear in cell count; ~10⁶ cells on a laptop per [3] | `harmony-integration-scanpy`, `harmony-integration`, or `seurat-integration` with `HarmonyIntegration` |
| **CCA / anchors (Seurat v3 default)** (Stuart 2019 [4]) | Canonical-correlation across pairs of samples; identifies mutual-nearest-neighbor "anchors"; builds correction vectors | Datasets where most cell types are shared across samples, small-to-medium scale, cross-modality (CITE / ATAC) integration via the same anchor framework | Aggressive — over-corrects rare or sample-specific populations; scales poorly above ~100k cells (CCA is quadratic-ish in pairwise-anchor work). Most modern best-practice guidance prefers Harmony / scVI as defaults [8]. # REVIEW(integration-methods): is CCA still the right *first* recommendation for any specific case in 2025-2026, or has it been fully superseded by Harmony+RPCA+scVI as the default tier? | Worst of the matrix at large N | `seurat-integration` with `CCAIntegration` (Seurat v5 default) |
| **RPCA (Seurat v5)** | Reciprocal PCA; faster + more conservative cousin of CCA — each pair-of-samples is projected into the other's PCA before anchor finding | Larger datasets where CCA is too slow; partly-shared cell-type composition; rare-population preservation matters | Less aggressive than CCA — under-corrects when batch effects are strong | Substantially faster than CCA at moderate N | `seurat-integration` with `RPCAIntegration` |
| **scVI** (Lopez 2018 [5]) | Variational autoencoder on raw counts; conditions on batch covariate; outputs a per-cell latent representation `X_scVI` that you treat like a PCA embedding | Atlas-scale; complex/nested batch effects (multiple technologies); downstream generative tasks (scANVI label transfer, scArches reference mapping, scVI DE); when an opaque-but-powerful embedding is OK | Heavy compute (effectively GPU-bound for any non-tiny dataset); failure mode is opaque ("did it work?" requires sanity checks); needs a `batch_key` decision the user must own | GPU minutes-to-hours; CPU days-to-impractical | `scvi-integration` |
| **scANVI** (label-aware scVI) | scVI + cell-type labels as a semi-supervised signal so integration is constrained to preserve known biology | When you have trusted labels (atlas references, reliable manual annotation) and want integration to be label-aware — Luecken 2022 [1] ranks it as a top overall performer on complex tasks; sc-best-practices [8] explicitly recommends it when labels are available | No labels available; labels are unreliable (then it just confidently propagates the bad labels) | Same as scVI | `scvi-integration` then scANVI fit (see `bp-data-integration`) |
| **conos** | Builds per-sample expression graphs, joins them into a cross-sample mNN graph; no shared embedding by default — joint clustering / label propagation / DE happen ON the joint graph | Trajectory / continuous-state work where forcing a shared embedding distorts geometry; label propagation across samples; you want to preserve per-sample structure | Many downstream tools (CellxGene, scIB metrics, Seurat-only workflows) expect an integrated *embedding*, not a joint graph. Less widely used → less community tooling | Linear-ish; pre-process per sample then graph-join | `conos-integration` |
| **FastMNN** (Haghverdi 2018 / batchelor) | Mutual-nearest-neighbor batch correction in PCA space; outputs corrected coordinates | Lightweight R/Bioconductor session; pairwise batch correction; small-medium scale | Older method; community has largely shifted to Harmony / scVI / Seurat v5; rarely the top pick in recent benchmarks. # REVIEW(integration-methods): worth keeping in §4 or move to §7? | Linear-ish | No dedicated recipe — `bp-data-integration` mentions it; route via batchelor if needed (recipe gap — see REVIEW_LOG) |
| **BBKNN** (Polanski 2020 [6]) | Modifies the kNN graph construction to balance neighbors across batches; corrects the *graph* only, not an embedding or expression | Massive datasets where you only need the graph (clustering, UMAP) — orders-of-magnitude speedup over MNN-style methods [6]; "at least some cells of the same type exist across batches" assumption holds | You need a corrected expression matrix or a continuous embedding for downstream tasks (BBKNN doesn't produce one [8]); deep imbalance | Fastest in the matrix; linear in cell count | `bp-data-integration` (multi-method Python flow includes BBKNN) — no dedicated single-method recipe (gap — see REVIEW_LOG) |
| **Scanorama** (Hie 2019 [7]) | Pairwise nearest-neighbor matching across all dataset pairs ("panorama stitching"); produces both an integrated embedding and corrected expression | Heterogeneous collections (datasets that don't all share the same cell types — [7]'s headline claim); atlas-scale (1M cells in ~9 h per [7]); high-ranking on Luecken 2022 "complex" tasks [1] | When all samples share the same composition, Harmony/scVI are simpler and usually as good or better; opaque parameter surface | Hours at atlas-scale; faster than CCA, slower than BBKNN | `bp-data-integration` (multi-method) — no dedicated single-method recipe (gap — see REVIEW_LOG) |

Notes on the matrix:
- **Luecken 2022 [1] top tier on complex integration tasks: scANVI, Scanorama,
  scVI, scGen** — verbatim from the benchmark's RNA-task ranking. On *simple*
  batch tasks Harmony and Seurat were also "consistently top performers"
  per [8]. This is why Harmony stays the default in the Quick-decision
  table even though it does not top the complex-task ranking.
- **scGen** is omitted from the matrix because it requires labels AND a
  perturbation-modeling framing; for most users the relevant label-aware
  pick is **scANVI**, which we recommend instead. # REVIEW(integration-methods):
  is this scGen-omission defensible, or should the matrix include it for
  completeness?
- **LIGER (rliger)** is omitted: long-standing method (Welch 2019), but
  Luecken 2022 [1] does not place it in the top tier for either simple or
  complex RNA tasks, and ABA has no dedicated LIGER recipe. Mention in §7
  for completeness.

For per-method internals, the executable recipes themselves carry the
deep-dive (parameter knobs, output schema) — see the recipe pointers above.

## Anti-patterns

- **`sc.concat()` then cluster, calling it "integration".** Concatenating count
  matrices without a batch-aware correction is the most common mistake. The
  joint PCA then separates cells by sample (batch), not by cell type. Use one
  of the methods in §4 — never raw concat-then-cluster.
  (This guardrail is already documented in `scrna_pipeline.md` — it stays
  here too because §5 is where users look when picking a method.)
- **Integrating over the biological variable.** If `condition` (stim vs
  ctrl, disease vs healthy) IS the question, integrating over it erases
  the signal. Integrate over the *technical* nuisance variable
  (donor/lane/run/day) and keep the biological variable as a separate
  covariate or contrast. Cited as a top failure mode in [8].
- **Reporting "integrated counts" for differential expression.** Every
  method except Seurat v4 `IntegrateData` corrects an *embedding* (Harmony
  → `X_pca_harmony`; scVI → `X_scVI`; BBKNN → kNN graph), not the count
  matrix. Cluster + visualize on the integrated space; do DE on the *raw*
  counts with sample as a covariate (per-cell) — and prefer pseudobulk
  DESeq2/edgeR for cross-condition contrasts (Squair 2021; ABA memory item
  `scrna_de_method`).
- **Picking a method by paper recency.** "scANVI is newer so it must be
  better" is not a benchmark. scANVI is in fact top-tier on Luecken 2022's
  *labelled complex* tasks [1], but that recommendation is conditional on
  *having reliable labels* — without them, scVI/Harmony/Scanorama are the
  defensible picks.
- **Skipping per-sample QC then "integrating away the difference".**
  Integration hides per-sample quality failures rather than fixing them.
  A high-MT or low-UMI sample looks like a "batch effect" the method
  obediently corrects. Per-sample QC first, integration second.
- **Choosing CCA because "Seurat tutorials use it".** Seurat v5's tutorials
  default to `CCAIntegration` for backward compatibility, but the same
  `IntegrateLayers` call accepts `RPCAIntegration` (faster, more
  conservative) and `HarmonyIntegration` (Luecken/sc-best-practices
  consensus default). The drop-in cost is one keyword argument; don't pick
  CCA just because it's the default.

## Sanity checks — how to know your choice was right

After running the integration:

- **UMAP coloured by sample.** Samples should be *mixed* across clusters,
  not segregated. A cluster that's 95% one sample is either real
  sample-specific biology (confirm it; e.g. matched tumor-only cluster) or
  an integration failure (try a more aggressive method, OR check sample
  QC).
- **UMAP coloured by known cell type / canonical marker.** Same cell type
  from different samples should co-locate; canonical markers should be
  contiguous on the UMAP rather than fragmented by sample.
- **Cluster-marker stability per sample.** The top markers per integrated
  cluster should match the top markers when each sample is analyzed
  separately. If markers drift dramatically post-integration, the method
  is altering biology — switch to a more conservative method (Harmony
  with higher `theta`; or RPCA instead of CCA).
- **Quantitative mixing metrics — scIB.** The `scib-metrics` package
  (Luecken 2022 [1]) reports both batch-removal scores (graph iLISI,
  kBET, ASW-batch) and bio-conservation scores (ARI/NMI vs known labels,
  cLISI, ASW-cell-type). The same biological dataset can rank methods
  differently — that's expected. Pick the method whose *bio-conservation*
  is high enough, then prefer the one with better batch removal among
  acceptable bio-conservation candidates. `bp-data-integration` automates
  this comparison.
- **Don't trust just the UMAP.** UMAP visualizes; it does not validate.
  Check the integrated embedding directly (`X_pca_harmony`, `X_scVI`)
  via scIB metrics or a per-sample marker-overlap audit.

## See also

**Benchmarks + reviews:**
1. Luecken et al. 2022, Nat Methods — *Benchmarking atlas-level data
   integration in single-cell genomics*; 16 methods × 13 tasks × 1.2M
   cells. (https://doi.org/10.1038/s41592-021-01336-8)
2. Tran et al. 2020, Genome Biology — *A benchmark of batch-effect
   correction methods for single-cell RNA sequencing data*; older but
   still cited. (https://doi.org/10.1186/s13059-019-1850-9)
   # REVIEW(citations): I could not retrieve the Tran 2020 abstract via
   # WebFetch (Springer auth wall). The references to it in §4 are
   # general/community-knowledge; an expert should verify they hold or
   # remove the inline citation.

**Method papers:**
3. Korsunsky et al. 2019, Nat Methods — Harmony.
   (https://doi.org/10.1038/s41592-019-0619-0)
4. Stuart et al. 2019, Cell — Seurat v3 / CCA / anchors.
   (https://doi.org/10.1016/j.cell.2019.05.031)
5. Lopez et al. 2018, Nat Methods — scVI.
   (https://doi.org/10.1038/s41592-018-0229-2)
6. Polanski et al. 2020, Bioinformatics — BBKNN.
   (https://doi.org/10.1093/bioinformatics/btz625)
7. Hie et al. 2019, Nat Biotechnol — Scanorama.
   (https://doi.org/10.1038/s41587-019-0113-3)

**Community consensus:**
8. Heumos et al. 2023, Nat Rev Genet / *Single-Cell Best Practices*
   (sc-best-practices.org) — integration chapter
   (https://www.sc-best-practices.org/cellular_structure/integration.html);
   tiered recommendation (Harmony/Seurat for simple, scVI/scANVI/Scanorama
   for complex, scANVI when labels exist).

**Recipes that execute paths from §2 / §4:**
- `harmony-integration` (R/Seurat + harmony) — `aba-recipes/recipes/genomics/harmony-integration.md`
- `harmony-integration-scanpy` (Python/scanpy + harmonypy) — `aba-recipes/recipes/genomics/harmony-integration-scanpy.md`
- `seurat-integration` (R/Seurat v5 `IntegrateLayers` — CCA / RPCA / Harmony) — `aba-recipes/recipes/genomics/seurat-integration.md`
- `scvi-integration` (Python/scvi-tools) — `aba-recipes/recipes/genomics/scvi-integration.md`
- `conos-integration` (R/conos) — `aba-recipes/recipes/genomics/conos-integration.md`
- `bp-data-integration` (Python multi-method, scIB-ranked, includes
  BBKNN/Scanorama/scANVI alongside scVI/Harmony) —
  `aba-recipes/recipes/genomics/bp-data-integration.md`
- `create-harmony-embeddings-scrna` (embedding-only step when PCA already
  exists) — `aba-recipes/recipes/genomics/create-harmony-embeddings-scrna.md`

**Adjacent knowhow:**
- `scrna_pipeline.md` — upstream framing: multiple samples → keep
  separate by default; never naive concat.
- (gap) `scrna_de_methodology.md` — downstream framing: DE across
  conditions uses pseudobulk on raw counts, NOT per-cell on the
  integrated embedding. Memory item `scrna_de_method` covers this; a
  dedicated knowhow doc is the natural next addition.
  # TODO(knowhow): author scrna-de-methodology knowhow.

**Recipe gaps surfaced during authoring** (see `REVIEW_LOG.md`):
- No dedicated single-method **BBKNN** recipe (only the multi-method
  `bp-data-integration`). Users who want a BBKNN-only pipeline are
  currently routed via the best-practices recipe.
- No dedicated **Scanorama** recipe (same — only inside
  `bp-data-integration`).
- No dedicated **FastMNN / batchelor** recipe (mentioned in
  `bp-data-integration`).
- No dedicated **LIGER** recipe — and no current literature pressure to
  add one (it does not top recent benchmarks).
- No dedicated **scANVI** recipe — currently reached via `scvi-integration`
  then the scANVI step from `bp-data-integration`. Worth promoting to its
  own recipe given Luecken 2022 [1] top-tier ranking when labels exist.
