---
name: scrna-celltype-annotation
description: Decision guide for assigning cell-type labels to scRNA-seq clusters. Compares manual marker-based annotation (PanglaoDB/CellMarker), automated reference-correlation (SingleR), supervised reference mapping (Azimuth), pre-trained linear classifiers (CellTypist), and deep-learning reference mapping (scArches / scANVI). Use when the user has a clustered scRNA-seq object and is choosing how to label the clusters.
when_to_use: User asks "what cell types are these clusters", "how do I annotate my scRNA-seq", "should I use SingleR / Azimuth / CellTypist / scANVI", or "manual markers vs automated annotation". The deliverable is a recommendation, not a how-to walkthrough.
avoid_when: User wants the actual execution (point them at the recipe named in §4); user is doing spatial or non-scRNA-seq annotation (different references and assumptions); user is in QC/clustering territory (point at `scrna_qc_clustering` / `seurat-scrna-v2`).
invocation: interactive
kind: knowhow_draft
requires_tools: [WebFetch, Read]
keywords: [cell type annotation, scRNA-seq, marker genes, PanglaoDB, CellMarker, SingleR, Azimuth, CellTypist, scArches, scANVI, panhumanpy, reference mapping, label transfer, decision guide]
produces: []
capabilities_needed: []
domain: genomics
source: "Aran 2019 (SingleR), Hao 2021 (Azimuth/WNN), Lotfollahi 2022 (scArches), Domínguez Conde 2022 (CellTypist), Pasquini 2021 + Abdelaal 2019 benchmarks, PanglaoDB (Franzén 2019), CellMarker 2.0 (Hu 2023), plus ABA's seurat-scrna-v2 Step 7 and bp-annotation for the executable paths."
audience: both
---

# scRNA-seq cell-type annotation — manual markers, automated references, or deep learning?

The question this knowhow answers: **"My scRNA-seq is clustered. How should I
assign cell-type labels — eyeball the markers against PanglaoDB, run SingleR
or Azimuth on a reference, train CellTypist, or do deep-learning mapping with
scArches / scANVI?"**

Audience: both — the naive biologist who has clusters and wants names on
them, and the bioinformatician disambiguating between annotation tools for
a specific dataset shape (rare tissue, atlas-scale query, cross-species,
disease state). The decision is multi-axis: reference availability × tissue
specificity × pipeline (R vs Python) × confidence requirement
(exploratory vs paper-grade) × granularity (broad vs fine subtypes). The
executable recipes this knowhow points at are `seurat-scrna-v2`,
`annotate-celltype-with-panhumanpy`, `bp-annotation`,
`scvi-label-transfer-scanvi`, `scvi-reference-mapping`,
`annotate-celltype-scrna`, and `unsupervised-celltype-transfer-between-scrna-datasets`.

## Quick decision

For the typical case (human or mouse tissue, clustered scRNA-seq, you want
both broad and fine labels you can defend in a paper):

| Your situation | Use |
|---|---|
| **Default first pass — always do this** | **Manual marker check** against PanglaoDB / CellMarker 2.0 / a published tissue atlas (`seurat-scrna-v2` Step 7 or `bp-annotation` §1). Even if you automate, you still validate manually. |
| Human PBMC / immune, Python pipeline, want a fast pre-trained call | **CellTypist** with `Immune_All_Low.pkl` (`bp-annotation` §2). Linear classifier, multi-resolution, majority voting. |
| Human PBMC / multimodal reference exists, R pipeline | **Azimuth** (PBMC / lung / motor cortex / heart references). Currently invoked in ABA via `annotate-celltype-with-panhumanpy` (the pan-human Azimuth neural net) or in Seurat directly. |
| Any tissue with a well-curated bulk/sorted reference, R / Bioconductor pipeline | **SingleR** against a Bioconductor reference (HumanPrimaryCellAtlas, MonacoImmune, BlueprintEncode, ImmGen for mouse). # REVIEW(annotation-methods): no dedicated `singler-celltype-annotation` recipe exists in ABA yet — this row currently lacks an executable pointer. |
| You have your OWN annotated reference (or a saved scVI / scANVI hub model) and a new query you want mapped onto it | **scANVI label transfer** (`scvi-label-transfer-scanvi`) if you can retrain together; **scArches reference mapping** (`scvi-reference-mapping`) if the reference is frozen. |
| Atlas-scale query, cross-batch / cross-tissue, want uncertainty propagated | **scArches + scANVI**, with explicit uncertainty thresholding (`scvi-reference-mapping`). |
| You want an ensemble that hedges between methods | **popV** (`unsupervised-celltype-transfer-between-scrna-datasets`) — runs CellTypist + KNN-scVI + KNN-Harmony + scANVI + others and votes. |
| Rare tissue with no published reference; novel or developmental cell types | **Manual only.** Don't trust any automated tool here — see §3. |

If your situation isn't in the table OR your dataset has unusual features
(non-model organism, heavy disease perturbation, very sparse depth), read
§4 (Alternatives matrix) before picking.

**The non-negotiable rule:** every automated label is a hypothesis until
you have eyeballed the canonical markers for that population on the same
clusters (Hao 2021; Heumos 2023). Section 5 spells out the anti-patterns
that bite users who skip this.

## When NOT to do automated cell-type annotation

- **Rare tissue with no published reference.** If nobody has produced an
  annotated single-cell atlas of your tissue/organ/developmental stage,
  there is nothing for SingleR / Azimuth / CellTypist / scArches to map
  onto. Running them anyway will return the closest match in whatever
  reference you forced — usually wrong, sometimes confidently wrong.
  Annotate manually from canonical markers + a published bulk / sorted
  cell-type reference for the constituent lineages.
- **You suspect novel cell types or cell states not in any reference.**
  Tumor cells, developmental intermediates, perturbation-induced states,
  rare disease populations. Reference-based methods will silently map
  these onto the nearest "normal" cell type. Manual marker analysis is
  the only path; treat reference-based labels (if you run them) as a
  one-line sanity-check, never as ground truth.
- **Very high inter-sample heterogeneity (e.g. heavy disease, perturbation,
  patient-specific subclones).** Reference mapping assumes the query's
  underlying biology matches the reference's; large condition shifts
  break that assumption (Lotfollahi 2022 motivates scArches precisely
  because vanilla mapping fails on out-of-distribution queries). Even
  with scArches, validate per-condition that mapped labels make
  biological sense.
- **You haven't clustered (or clustered badly).** Annotation operates on
  clusters or on per-cell calls projected onto clusters; bad clustering
  → bad annotation. Fix QC + clustering first (`scrna_qc_clustering` /
  `seurat-scrna-v2`).
- **You're asking the wrong question.** "What are these cells?" sometimes
  reframes to "what STATE are these cells in?" — interferon response,
  cell cycle phase, stress signature. Those aren't cell types; they're
  states. Use module scoring (`sc.tl.score_genes`, `AddModuleScore`),
  not cell-type annotation.

## Alternatives — full matrix

| Method | Characterization | When it wins | When it fails | Cost | Recipe |
|---|---|---|---|---|---|
| **Manual marker panels** (PanglaoDB, CellMarker 2.0, atlas paper marker lists) | Score curated marker sets per cluster (dotplot / module score), assign by inspection. Database-only or atlas-derived gene lists. | Always the floor — required regardless of automated method chosen. Wins outright when no reference exists, or when the user needs paper-grade annotation. | Slow; subjective; misses fine subtypes (the dotplot resolution caps at "I can see a difference"). Sensitive to marker-set quality — PanglaoDB has known false-positive markers. # REVIEW(annotation-methods): is the Franzén 2019 PanglaoDB marker set still considered the canonical reference for broad lineages, or has CellMarker 2.0 (Hu 2023) overtaken it for breadth? | Minutes for the dotplot + human time | `seurat-scrna-v2` Step 7 (R); `bp-annotation` §1 (Python/scanpy) |
| **SingleR** (Bioconductor) | Per-cell Spearman correlation between each cell's expression and bulk/sorted reference profiles; assigns highest-scoring label. Reference-correlation, not classifier. | Bioconductor pipeline; broad lineage calls (T/B/NK/Mono/etc.); when good bulk-sorted refs exist (HumanPrimaryCellAtlas, MonacoImmune, BlueprintEncode, ImmGen) (Aran 2019). | Fine subtypes the reference can't resolve; rare tissues without a curated SingleR reference; the per-cell pruned-vs-unpruned distinction is easy to misinterpret. | Linear in cells × ref labels; minutes to ~hour on 100k cells | # TODO(recipe): `singler-celltype-annotation` — no dedicated ABA recipe yet; users can call SingleR ad-hoc from `run_r` against the Bioconductor `celldex` references |
| **Azimuth** (Seurat lab; supervised, supports WNN multimodal refs) | Anchor-based label projection onto a Seurat-lab curated reference; returns labels at multiple granularities + a mapping-score + a prediction-score per cell (Hao 2021). | The reference's exact tissue (PBMC, lung, motor cortex, heart, bone marrow, pancreas, kidney, tonsil); user wants hierarchical labels + per-cell confidence; CITE-seq queries that match a WNN reference. | Off-reference tissue / species; query batches that drift from reference; reliance on Azimuth's confidence score WITHOUT marker validation (see §5). | Minutes for the standard refs; longer for the pan-human ANN | `annotate-celltype-with-panhumanpy` (the pan-human Azimuth neural net via panhumanpy); ad-hoc Seurat `MapQuery` for the original tissue-specific refs |
| **CellTypist** (consortium-trained linear logistic-regression classifiers + majority voting) (Domínguez Conde 2022) | Pre-trained per-tissue models (Immune_All_Low/High, Pan_Fetal_Human, Adult_Mouse_Gut, etc.); CP10k + log1p input; majority voting smooths per-cell calls within local neighborhoods. | Python pipeline; immune tissues (the original publication target); fast scaling; when you want a calibrated confidence column out of the box. | Tissues outside the published model catalogue; severe batch effects relative to model training data; when fine subtypes you need aren't in the model's output classes. | Seconds-to-minutes on 100k cells | `bp-annotation` §2; ensemble path in `unsupervised-celltype-transfer-between-scrna-datasets` |
| **scANVI** (scvi-tools; semi-supervised VAE classifier head on top of scVI) | Trains a classifier inside the scVI latent space using cells that ARE labeled; predicts the unlabeled ones in the SAME integrated embedding. | You have a partially-labeled object (some cells annotated, others Unknown) and want labels propagated within the same dataset; you control both reference and query end-to-end. | You need to map onto a SAVED reference you can't retrain — use scArches instead; weak when labeled set is too small or too imbalanced. | GPU recommended; tens of minutes to hours | `scvi-label-transfer-scanvi` |
| **scArches** (architectural surgery on a frozen scVI/scANVI reference) (Lotfollahi 2022) | Freezes the reference encoder, adds per-query-batch weights, fine-tunes only those. Maps a new query into the reference's latent space without retraining the reference; with a scANVI reference, transfers labels + uncertainty. | Atlas-scale references (HLCA, Tabula Sapiens, etc.); you want to add a new sample to a published atlas; cross-condition queries where re-training the full reference is impractical. | Reference not trained scArches-ready (`encode_covariates=True`); query too far out of distribution; user trusts the mapped labels without uncertainty thresholding. | GPU recommended; query-fine-tune is faster than full scVI train | `scvi-reference-mapping` |
| **popV ensemble** | Runs CellTypist + KNN-scVI + KNN-Harmony + KNN-Scanorama + scANVI + Random-Forest + XGBoost + Onclass against a reference, then votes per cell. | When you want a hedged label with a built-in disagreement signal; reproducibility studies. | Heavy compute; ensembles hide the failure modes of individual methods if you don't inspect the per-method outputs. | High — runs multiple methods | `unsupervised-celltype-transfer-between-scrna-datasets` |
| **LLM-from-markers** (ABA-specific) | Per-cluster top markers fed to an LLM constrained to the CZI census cell-type vocabulary. Hypothesis-generation. | Quick first-pass naming on small numbers of clusters when no obvious reference fits; useful as a marker-narrative generator. | NOT a substitute for any of the above. Treat output as a marker-summarization aid; ALWAYS validate against the dotplot and the literature. | LLM tokens + the dotplot | `annotate-celltype-scrna` |
| **scGate** (narrow alternative) | Hierarchical marker-gate purification — designed to *gate* populations (e.g. "purify T cells") more than to annotate the whole object. | Cleaning a population for downstream analysis; sanity-check on a single lineage. | Not a whole-object annotation tool; mention rather than recommend for the general question. | Minutes | # TODO(recipe): no scGate recipe yet |
| **GeneCT** (narrow alternative) | Gene-centric classifier trained on cell-type-specific gene panels. | Niche use; less widely adopted than the methods above. | Limited tissue coverage relative to CellTypist / Azimuth. # REVIEW(annotation-methods): is GeneCT actively maintained as of 2025? If not, drop from §4. | Modest | # TODO(recipe): no GeneCT recipe yet |

## Anti-patterns

- **Naming a cell type from a SINGLE marker.** "FCGR3A+ monocytes" needs
  FCGR3A *plus* the rest of the non-classical monocyte signature (LST1,
  CDKN1C, MS4A7, FCN1-low). One marker is a hypothesis; the full
  signature is the label. The `seurat-scrna-v2` Step 7 callout names
  this anti-pattern explicitly.
- **Over-trusting Azimuth (or CellTypist, or scANVI) confidence scores.**
  The score is calibrated against the reference's distribution; if your
  query is out of distribution, the score is over-confident. Validate
  every "high-confidence" call against the marker dotplot for that
  population. Heumos 2023 (single-cell best practices) says automated
  calls are a *starting point*, not the answer.
- **Running a brain reference on intestinal data (or any tissue
  mismatch).** Reference-based methods will return SOMETHING — usually
  the nearest available label, which can be biologically nonsensical.
  Tissue must match reference tissue, or you're doing pan-tissue
  annotation (Tabula Sapiens / pan-human Azimuth) which is the
  exception, not the default. # REVIEW(annotation-methods): is there a
  named tissue-mismatch failure paper to cite here, or is this
  community consensus only?
- **Assigning labels per-cell when clusters are coherent.** Cluster-level
  annotation is more robust than per-cell, given marker sparsity and
  dropout (`bp-annotation` §1). Per-cell labels are tempting because
  CellTypist/Azimuth give them; aggregate to cluster majority before
  reporting, except in the cases where per-cell label is the point
  (e.g. detecting label heterogeneity within a cluster).
- **Forgetting to inspect mapping quality before transferring labels.**
  scArches and Azimuth both expose per-cell mapping-quality / prediction
  scores. Cells below threshold should be "Unknown", not the
  best-effort label. The book threshold for scArches is uncertainty
  > 0.2 → "Unknown" (Heumos 2023).
- **Ensembling without inspecting per-method disagreement.** popV will
  return a vote even when its members disagree wildly. Treat per-cell
  method disagreement as a flag for manual review, not noise to be
  averaged out.
- **Using the LLM-marker-naming recipe as ground truth.** ABA's
  `annotate-celltype-scrna` produces an LLM-narrated name from the top
  markers — it is a hypothesis generator that should feed §1 (manual
  validation), never replace it. The recipe itself warns that LLM
  labels are hypotheses.
- **Re-labeling on the integrated counts.** Annotation, like DE, should
  use the corrected EMBEDDING (clusters from the integrated space) but
  the original counts for marker expression. Don't run
  `rank_genes_groups` / `FindAllMarkers` on the batch-corrected matrix.
- **Skipping QC and re-clustering after relabeling.** If annotation
  reveals that "cluster 7" is two populations smushed (e.g. CD4 and
  CD8 T cells co-clustered), the answer is sub-cluster, not "annotate
  it ambiguously".

## Sanity checks

After running your chosen annotation:

- **Per-cluster dotplot of canonical markers for the assigned labels.**
  Every label should be defensible from this single figure. If T cells
  don't show CD3D/CD3E/CD8A/CD4 or B cells don't show CD79A/MS4A1, the
  label is wrong.
- **UMAP coloured by assigned label.** Labels should form contiguous
  patches on the embedding. Salt-and-pepper labels (one cell-type call
  scattered across clusters) usually mean per-cell calls without
  majority-voting / smoothing.
- **Confidence column inspected.** CellTypist `conf_score`, Azimuth
  `mapping.score` + `predicted.celltype.l2.score`, scArches/scANVI
  prediction probability. Plot the histogram; if everything is at the
  ceiling, the model is over-confident and you should down-weight it.
- **Cluster-level majority + uncertainty.** Per cluster, what fraction
  of cells got the dominant label? <70% means the cluster is mixed
  (sub-cluster) or the call is uncertain (mark "Unknown").
- **Cross-check between methods.** Run CellTypist and a manual marker
  dotplot, or SingleR and Azimuth, on the same object. Agreement at
  the broad-label level (lymphoid / myeloid / stromal) is the floor;
  disagreement at the fine level is informative — investigate, don't
  pick the prettier answer.
- **Held-out marker sanity.** Pick a marker NOT used in any reference
  the method consumed (e.g. a recent paper's marker for a subtype the
  reference may pre-date) and check that it lights up the expected
  cluster.
- **Map back to biological expectation.** For PBMC: ~60% T, ~15% B,
  ~15% mono, ~5–10% NK is the textbook range. Wildly off proportions
  (e.g. "90% NK") usually means a labeling bug, not biology.

## See also

**Method papers (primary):**
- Aran et al. 2019, *Nat Immunol* — SingleR
  (https://doi.org/10.1038/s41590-018-0276-y)
- Hao et al. 2021, *Cell* — Azimuth + WNN reference paper (Seurat v4)
  (https://doi.org/10.1016/j.cell.2021.04.048)
- Lotfollahi et al. 2022, *Nat Biotechnol* — scArches architectural surgery
  (https://doi.org/10.1038/s41587-021-01001-7)
- Domínguez Conde et al. 2022, *Science* — CellTypist + cross-tissue immune
  cell atlas (https://doi.org/10.1126/science.abl5197)
- Xu et al. 2021, *Mol Syst Biol* — scANVI semi-supervised model
  (https://doi.org/10.15252/msb.20209620)
  # REVIEW(citations): confirm Xu 2021 is the correct primary citation for scANVI vs the scvi-tools framework paper (Gayoso 2022).

**Benchmarks + reviews:**
- Pasquini et al. 2021, *Comput Struct Biotechnol J* — automated cell-type
  annotation benchmark across methods
  (https://doi.org/10.1016/j.csbj.2021.01.015)
- Abdelaal et al. 2019, *Genome Biology* — earlier comparison of
  automated annotation tools
  (https://doi.org/10.1186/s13059-019-1795-z)
- Heumos et al. 2023, *Nat Rev Genet* / Single-cell Best Practices —
  annotation chapter (https://sc-best-practices.org/cellular_structure/annotation.html)
  # REVIEW(annotation-methods): is there a 2023–2025 benchmark that has
  superseded Pasquini 2021 / Abdelaal 2019 for annotation method choice?

**Canonical marker catalogues:**
- Franzén et al. 2019, *Database* — PanglaoDB
  (https://doi.org/10.1093/database/baz046)
- Hu et al. 2023, *Nucleic Acids Res* — CellMarker 2.0
  (https://doi.org/10.1093/nar/gkac947)

**Recipes that execute paths from §2 / §4:**
- `seurat-scrna-v2` (Step 7) — R/Seurat marker discovery + manual annotation
- `bp-annotation` — Python best-practice annotation: manual + CellTypist + scArches/scANVI
- `annotate-celltype-with-panhumanpy` — pan-human Azimuth neural net
- `annotate-celltype-scrna` — LLM-from-markers hypothesis (always validate manually)
- `scvi-label-transfer-scanvi` — scANVI semi-supervised label transfer
- `scvi-reference-mapping` — scArches reference mapping onto a frozen scVI/scANVI atlas
- `unsupervised-celltype-transfer-between-scrna-datasets` — popV ensemble (CellTypist + KNN-scVI + KNN-Harmony + scANVI + …)

**Coverage gaps (no recipe yet — flagged for follow-up):**
- `# TODO(recipe): singler-celltype-annotation` — dedicated SingleR + celldex path is a known gap; users currently call SingleR ad-hoc from `run_r`.
- `# TODO(recipe): azimuth-seurat-direct` — Seurat-native Azimuth (`MapQuery` against the tissue-specific Seurat-lab references) currently only available through ad-hoc `run_r`; pan-human is covered by panhumanpy.
- `# TODO(recipe): scgate-purify` — narrow alternative, useful for population purification rather than whole-object annotation.

**Adjacent knowhow:**
- `scrna_pipeline.md` — upstream framing (do not naively concatenate; per-sample QC first; integrate before annotating multi-sample objects).
- `bulk_rnaseq_de.md` — out of scope; bulk DE has nothing to do with cell-type annotation but is the closest sibling currently in `knowhow/`.
- (Future) `scrna-integration.md` — annotation usually follows integration; the integration choice affects which cells co-cluster and therefore what gets annotated together.
- (Future) `scrna-de-methodology.md` — once labels are assigned, DE between conditions within a cell type lives there (pseudobulk, not per-cell, for cross-condition tests; Squair 2021).
