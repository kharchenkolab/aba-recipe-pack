---
name: scrna-analysis
description: Top-level navigation hub for single-cell RNA-seq (scRNA-seq) analysis. Walks through the pipeline stages from raw data to biological interpretation; at each stage names the decision the user faces, points at the specific knowhow doc covering that decision in depth, and names the canonical executable recipe(s). Read THIS first when starting any scRNA-seq analysis or when unsure where to find the right method-choice guidance.
when_to_use: User has scRNA-seq data (single sample, multi-sample, integrated, multimodal) and is figuring out how to analyze it. Or user is mid-analysis and needs to decide what to do next (cluster? integrate? annotate?). Or agent is helping a biologist orient across the many possible paths through an scRNA-seq pipeline. This hub points outward to specific knowhow + recipes — it does NOT replace them.
avoid_when: User has already picked a path and wants step-by-step execution — go straight to the named recipe. User is asking a specific methodology question (e.g. "Wilcoxon vs MAST?") — go straight to the relevant per-class knowhow (e.g. scrna-de-methodology).
invocation: interactive
kind: knowhow_draft
requires_tools: [Skill, Read]
keywords: [scRNA-seq, single cell, pipeline overview, navigation, where to start, analysis hub, stage selection, scrna analysis steps, what next]
domain: genomics
source: "Top-level pipeline framing drawn from Heumos et al. 2023 (Nat Rev Genet) sc-best-practices, Luecken & Theis 2019 (Mol Syst Biol) Current best practices, and the structure established by aba-recipes' existing per-stage knowhow docs (scrna_pipeline.md, scrna-de-methodology, scrna-integration-decision, etc.)."
audience: both
---

# scRNA-seq analysis — pipeline navigation hub

This is the **entry point** for single-cell RNA-seq analysis. The pipeline has
many stages, and at most stages the user has a real choice between methods.
This doc maps the stages, names the decisions, and routes you to the right
deeper material:

- For **principles + opinion** on a decision → the relevant **knowhow** doc
- For **executable workflow** at a stage → the relevant **recipe**

This is not a tutorial. It's a navigator. If you're at stage X with question
Y, this doc tells you where the answer lives.

## Quick orientation — where are you in the pipeline?

```
                              ┌──── doublet detection ─────┐
[load 10x]──[QC]──[filter]──┬─┤                             ├──[normalize]──[HVG]
                            │ └─── (optional, post-filter) ─┘
                            ↓
              [single sample? multi-sample?]
                            ↓
             ┌─── single: continue per-sample ──→ [PCA]──[cluster]──[UMAP]──[markers]──[annotate]
             ↓
             └─── multi: [integrate? or not?]
                            ↓
                            ├── NO → per-sample analysis, compare clusters across
                            └── YES → [pick method]──[integrate]──[joint cluster]──[markers]──[annotate]
                                                                         ↓
                                                                  [DE: per-cell? pseudobulk?]
                                                                         ↓
                                                                  [trajectory? cell-cell? etc.]
```

Each stage decision below has its own knowhow doc with the depth + citations.

## Pipeline stages — decision table

For each stage: the decision you face, the knowhow doc with the answer, and
the recipe that executes it.

| Stage | Decision | Knowhow (depth) | Recipe (execute) |
|---|---|---|---|
| **0. Experimental design** | How many samples? How many cells per sample? Loading concentration? | (out of scope here; consult statistician + Heumos 2023 §"Study design") | — |
| **1. Load + initial QC** | What thresholds for nFeature / nCount / percent.mt? Per-tissue defaults? | `scrna_pipeline.md` (general principles); each recipe's QC step is data-driven from quantile tables | `seurat-scrna-v2`, `scrna-qc-clustering` (scanpy) |
| **2. Sample-structure** | Single-sample workflow, OR multi-sample-aware from the start? | `scrna_pipeline.md` — **CRITICAL rule**: multiple samples are NOT one matrix; don't `sc.concat()` raw then cluster | Per-sample first via `seurat-scrna-v2` / `scrna-qc-clustering`; then revisit at stage 5 |
| **3. Doublet detection** | Do I need one? Which caller? Pre- or post-empty-droplet filter? | `scrna-doublet-detection` | TBD recipe — currently inline in QC steps |
| **4. Normalization** | LogNormalize, SCTransform, or Pearson residuals? | `scrna-normalization-decision` | `seurat-scrna-v2` (LogNormalize), `seurat-sctransform` (SCT) |
| **5. Integration — yes or no?** | Is the sample-separation in my UMAP batch or biology? | `scrna-batch-vs-condition` — diagnostic framework | (decision, not execution) |
| **6. Integration method** | Harmony, CCA, RPCA, scVI, conos, FastMNN — which? | `scrna-integration-decision` | `seurat-integration`, `harmony-integration-scanpy`, `scvi-integration` |
| **7. Clustering** | Louvain, Leiden, SLM? What resolution? | (mostly recipe-level; method comparison is mature, resolution is empirical) | Recipe steps cover this — `seurat-scrna-v2` Step 6 documents Louvain default + Leiden alternative |
| **8. Cell-type annotation** | Manual markers, SingleR, Azimuth, scArches, CellTypist? | `scrna-celltype-annotation` | `seurat-scrna-v2` Step 7 (manual); `seurat-multimodal-reference-mapping` (Azimuth); standalone SingleR recipe TBD |
| **9. Marker discovery** | Which test? Which thresholds? | `scrna-de-methodology` (covers markers as one of three paths) | `seurat-scrna-v2` Step 7, `seurat-de-testing` Path A |
| **10. Differential expression** | Per-cell or pseudobulk? Which test? | `scrna-de-methodology` — **CRITICAL rule**: pseudobulk for condition effects across samples; per-cell DESeq2 is wrong (Squair 2021) | `seurat-de-testing` (all three paths) |
| **11. Multimodal extensions** | CITE-seq? Multiome (RNA+ATAC)? Spatial? | Stage-specific knowhow (TBD; ask in chat) | `seurat-cite-seq`, `seurat-wnn-multimodal`, `seurat-rna-atac-integration` |
| **12. Trajectory / pseudotime** | Monocle3, Slingshot, PAGA? | (knowhow TBD) | (recipe TBD) |
| **13. Cell-cell communication** | CellChat, NicheNet, LIANA? | (knowhow TBD) | (recipe TBD) |
| **14. View / share the result** | How do I explore it interactively? Which format to save / share / archive? | `scrna-viewing-and-interchange` | pagoda3 (built-in viewer); `get_viewer_url` |

For stages marked TBD, the depth doesn't yet exist in the catalogue. The
agent should call this out to the user rather than improvising.

## Cross-cutting decisions (apply across stages)

### Software ecosystem — Seurat (R) or scanpy (Python)?

Both ecosystems can do essentially the same analysis end-to-end. The choice
is mostly:

- **R/Seurat** — picked when downstream tools are Bioconductor (DESeq2,
  edgeR, ComplexHeatmap, …), when the lab is R-native, or when the user
  names Seurat. Strong default for clinical / biostatistics-heavy projects.
- **Python/scanpy** — picked when downstream tools are
  Python/PyTorch-native (scVI, scArches, CellTypist, cell2location, …),
  when the lab is Python-native, or when the user names scanpy. Strong
  default for atlas-scale / deep-learning-adjacent projects.

The recipes are organized so that the same analysis has parallel paths in
both — e.g. `seurat-scrna-v2` (R) and `scrna-qc-clustering` (Python). Pick
your ecosystem at the start and stick with it within a project.

### Multi-sample by default

This is in `scrna_pipeline.md` already but worth emphasizing here:
**multiple scRNA-seq samples are multiple datasets, not one matrix.** Never
`sc.concat()` raw count matrices into one AnnData and cluster — the joint
PCA then separates by sample, not by cell type. Either analyze per-sample
+ compare, or use a batch-aware integration recipe (stage 5–6 above).

### Modalities (RNA only, CITE-seq, multiome, spatial)

If your data is multimodal, several stages above branch:
- **CITE-seq (RNA + ADT protein)**: see `seurat-cite-seq` recipe + (when
  joint clustering is needed) `seurat-wnn-multimodal` for WNN
- **10x Multiome (RNA + ATAC)**: see `seurat-rna-atac-integration` recipe;
  the ATAC side uses Signac with TF-IDF + SVD/LSI rather than PCA
- **Reference mapping** for either modality: see
  `seurat-multimodal-reference-mapping`

Spatial transcriptomics is currently out of scope (was in the catalogue
briefly; dropped pending re-scoping).

## Anti-patterns at the pipeline level

These are mistakes the navigator's job is to prevent — they recur across
stages because they're about HOW the analysis is approached, not about a
specific method choice.

- **Concatenating samples and calling it integration.** Raw `sc.concat()`
  or `merge()` followed by joint PCA + clustering is NOT integration;
  it's the failure mode integration was invented to fix. If you have
  multiple samples and want joint analysis, use stage 5–6.
- **Picking a method by recency.** "scANVI is newer than Harmony so it
  must be better" is not a benchmark. Cite a benchmark
  (Luecken 2022, etc.) or default to a well-validated conservative
  method (Harmony for integration, Wilcoxon for markers, manual + a
  curated panel for annotation).
- **Integrating over the biological variable.** If condition X vs Y IS
  the biological question, integrating over `condition` removes the
  signal. Integrate over the technical variable (donor, lane, day),
  not the biological one. See `scrna-batch-vs-condition`.
- **Trusting the UMAP for biology.** UMAP is a visualization; clusters
  come from the SNN/wsnn graph. Read cluster identity from
  `obj$seurat_clusters` / `adata.obs['leiden']`, never by drawing
  polygons on the UMAP.
- **DESeq2 on per-cell counts for condition effects.** Squair 2021
  demonstrated 10–100× inflated FPR. The hard rule: bulk DE tools
  (DESeq2/edgeR/limma/pydeseq2) on per-cell counts is wrong. See
  `scrna-de-methodology` — the single most-load-bearing knowhow doc
  in this catalogue.
- **Skipping per-sample QC before integration.** Integration hides
  per-sample quality failures. Run the single-sample QC + clustering
  per sample first (`seurat-scrna-v2` or `scrna-qc-clustering` per
  sample); integrate after, not before.
- **Cell-type names from a single marker.** "FCGR3A+ monocytes" needs
  FCGR3A *plus* the rest of the non-classical-monocyte signature.
  See `scrna-celltype-annotation`.
- **Reporting integrated counts.** Integration corrects the EMBEDDING
  (PCA, UMAP), not the counts. Cluster + visualize on the integrated
  space; do DE on the raw counts via the appropriate per-cell or
  pseudobulk path.

## When to leave this hub

This is the scRNA hub. If your data is:

- **Bulk RNA-seq** → go to `bulk_rnaseq_de.md` knowhow + `deseq2-r` recipe
- **Spatial transcriptomics** → out of catalogue scope right now
- **ATAC-only (no RNA)** → see Signac single-modality vignettes; no
  catalogue recipe yet
- **Genomics other than RNA** → out of scope here

## See also

### Sub-knowhow (per-decision depth)

- `scrna_pipeline.md` — sample-structure rule (multiple samples → keep
  separate by default); the focused position that predates this hub
- `scrna-integration-decision` — which integration method for which
  dataset
- `scrna-de-methodology` — per-cell vs pseudobulk DE (the Squair 2021
  hard rule)
- `scrna-normalization-decision` — LogNormalize vs SCTransform vs
  Pearson residuals
- `scrna-doublet-detection` — which caller for which protocol
- `scrna-celltype-annotation` — manual vs SingleR vs Azimuth vs
  scArches vs CellTypist
- `scrna-batch-vs-condition` — diagnostic framework for "is this
  sample-separation batch or biology"

### Canonical executable recipes

- `seurat-scrna-v2` — R/Seurat single-sample QC + clustering + markers
- `seurat-sctransform` — R/Seurat with SCTransform normalization
- `seurat-integration` — R/Seurat v5 layer-based multi-sample
- `seurat-reference-mapping` — R/Seurat reference mapping (RNA only)
- `seurat-cite-seq` — R/Seurat CITE-seq (RNA + ADT)
- `seurat-wnn-multimodal` — R/Seurat WNN multimodal clustering
- `seurat-multimodal-reference-mapping` — R/Seurat WNN-reference mapping
- `seurat-rna-atac-integration` — R/Seurat + Signac 10x Multiome WNN
- `seurat-de-testing` — R/Seurat DE methodology (3 paths)
- `scrna-qc-clustering` — Python/scanpy single-sample equivalent of
  seurat-scrna-v2
- `harmony-integration-scanpy` — Python/scanpy + Harmony
- `scvi-integration` — Python/scvi-tools deep-generative integration
- `conos-integration` — R/conos sample-graph integration
- `deseq2-r` — bulk RNA-seq DE (NOT scRNA per-cell; see anti-patterns)

### Reading list (literature anchors)

- **Heumos et al. 2023, Nat Rev Genet** — "Best practices for single-cell
  analysis across modalities" — comprehensive review covering most
  stages above (https://www.nature.com/articles/s41576-023-00586-w).
  The companion sc-best-practices book is at
  https://www.sc-best-practices.org/
- **Luecken & Theis 2019, Mol Syst Biol** — "Current best practices in
  single-cell RNA-seq analysis: a tutorial"
  (https://doi.org/10.15252/msb.20188746). Older but still cited as
  the pipeline-shape paper.
- **Luecken et al. 2022, Nat Methods** — integration benchmark cited
  by `scrna-integration-decision`
- **Squair et al. 2021, Nat Commun** — DE-methodology benchmark cited
  by `scrna-de-methodology`

## REVIEW notes

# REVIEW(scrna-overall-framing): does the stage breakdown above match the
expected biologist mental model? An expert reviewer should sanity-check
that the stage ordering + decisions align with how they'd present the
pipeline in a lab onboarding doc. If a stage is missing (e.g. cell-cycle
regression, ambient RNA correction with SoupX/CellBender) or out of order,
flag it.

# REVIEW(catalogue-gaps): the table flags TBD knowhow + TBD recipes for
several stages (doublet detection, trajectory, cell-cell, ambient RNA).
These gaps are the next-priority knowhow + recipe to author. Confirm the
priority before scheduling.

# REVIEW(spatial-scope): spatial transcriptomics was dropped from the
catalogue during this session pending re-scoping. The hub currently says
"out of catalogue scope right now" — confirm that framing matches
intended product direction (it implies spatial is coming back later;
adjust if not).
