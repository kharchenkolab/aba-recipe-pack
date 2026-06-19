# Review queue — scrna-batch-vs-condition

**Draft status:** AUTHORED (awaiting expert review)
**Author:** Phase δ Worker 6 (Claude / aba-knowhow-authoring)
**Drafted on:** 2026-06-05
**Confidence:** MEDIUM — the analytical question is genuinely
case-dependent; the deliverable is a diagnostic FRAMEWORK rather than
a single decision rule. Where the framework takes a position
("default to assume batch until ruled out"), it is the author's
opinion grounded in the cited literature, not a settled benchmark
finding.
**Estimated review time:** 60-90 min (5 high-priority + 4 medium + 3
low markers, plus 3 coverage gaps).

## High priority — claims that change the recommendation if wrong

1. **`§Quick-decision default-ruling` ("assume batch until you can
   rule out the technical explanation")** — this is the author's
   framework position, not a benchmark finding. The justification
   given (publish-cost asymmetry: misread batch-as-biology survives,
   misread biology-as-batch is visible at sanity-check time) is a
   reasoning argument, not a citation. # REVIEW(integration-methods,
   scrna-biology): is "default to batch" the right asymmetry for the
   typical biologist workflow, or does it bias against detecting
   sample-specific biology? Affects: the entire framing of the doc.

2. **`§Diagnostic-2 marker-conservation reading`** — the doc claims
   "Top-20 marker overlap with Jaccard > ~0.5 indicates the same
   cell type." # REVIEW(scrna-de): is 0.5 the right threshold, or
   should it be 0.3 / 0.7? The number is the author's heuristic, not
   sourced. Affects: how often the diagnostic returns "batch" vs
   "biology." Replace with a cited threshold or remove the specific
   number.

3. **`§Diagnostic-4 designed-controls "strongest single signal"`** —
   the doc elevates biological replicates of the same condition as
   the decisive diagnostic. # REVIEW(scrna-biology): is this
   defensible when the replicate-vs-replicate signal exists in some
   cell types but not others (a common real outcome)? The framework
   currently routes mixed-regime to "integrate the technical
   covariate"; verify that's the right call.

4. **`§Anti-patterns "treat samples-mix-on-UMAP as proof"`** — the
   doc says UMAP-by-sample mixing can be over-correction in
   disguise. This is cited to Luecken 2022. # REVIEW(integration-
   methods): does Luecken 2022 specifically demonstrate over-
   correction producing apparently-mixed UMAPs, or is this an
   inference? Spot-check the paper's Fig 2/3 + discussion.

5. **`§Sanity-checks integrated-counts anti-pattern`** — the doc
   prohibits "reporting integrated counts" and routes DE to
   pseudobulk. # REVIEW(scrna-de): this position is consistent with
   Squair 2021 (cited in `scrna-de-methodology`) but the cross-ref
   to that knowhow is forward-looking. Confirm the prohibition is
   universal (some methods, e.g. scVI, do produce a corrected count-
   like output `scvi.model.SCVI.get_normalized_expression` — what's
   the guidance on those?).

## Medium priority — claims that change the explanation if wrong

6. **`§Diagnostic-1 cluster-resolution caveat`** — the doc suggests
   "0.5–1.0 leiden" as a moderate resolution. # REVIEW(scrna-
   methods): are these the right defaults for the typical 5k-50k-
   cell PBMC / tissue dataset? Cite a source or remove the specific
   numbers.

7. **`§Diagnostic-3 QC tracking interpretation`** — claims sample-
   pure clusters that track depth/MT are "QC artifact." # REVIEW(
   scrna-biology): is this always true? Some biological cell types
   genuinely have low depth (small lymphocytes, late erythroid).
   The doc acknowledges this but the headline reading may be too
   strong. Soften?

8. **`§Diagnostic-5 inline REVIEW marker on tissue-agnostic proxy`** —
   the doc asks whether co-expression module preservation across
   samples could serve as a tissue-agnostic proxy for "expected
   shared composition." # REVIEW(scrna-biology): is there a
   published method (BRGMM, hdWGCNA module preservation, ...) that
   already does this? If yes, name it.

9. **`§See-also entry on integrated-counts`** — wording "Integration
   corrects the embedding (PCA / UMAP / latent), not the count
   matrix." # REVIEW(integration-methods): scVI / scANVI DO produce
   a corrected count-like output. The statement is true for Harmony
   / CCA / RPCA but not strictly for VAE-based methods. Refine.

## Low priority — minor wording / cosmetic

10. **`§Frontmatter "audience: both"`** — the doc is dual-audience
    by design. Confirm the dual-track structure (Quick-decision for
    naive, full diagnostics for experienced) reads cleanly to both
    audiences in a real review session.

11. **`§Citations`** — Tran 2020 is now ~6 years old. # REVIEW(
    citations): is there a more recent multi-method integration
    benchmark that supersedes Tran for the "should I integrate"
    question specifically? scIB-metrics (2024+) extends Luecken;
    others?

12. **`§7 See-also cross-refs`** — the doc forward-refs to
    `scrna-de-methodology` and `scrna-integration-decision`. Verify
    these neighbors materialize so the cross-refs resolve.

## Coverage gaps surfaced during authoring

- `# TODO(recipe): scrna-batch-diagnostics` — the five diagnostics
  in §4 + LISI + kBET in §6 are currently inline snippets. A
  dedicated recipe that produces a one-shot diagnostic report
  (cross-tab heatmap, per-cluster marker comparison, QC overlay,
  replicate-vs-replicate, cell-type-score-per-sample, LISI, kBET)
  would close the gap. Until it exists, the user copy-pastes.

- `# TODO(knowhow): scrna-integration-decision` — the sibling
  knowhow ("WHICH integration method") is referenced but not yet
  authored at draft time. Confirm whether another Phase δ worker
  is producing it; if not, add to the queue. Without it, this
  knowhow's "if you decided to integrate" arm ends in the
  alternatives without a method-selection guide.

- `# TODO(catalogue): seurat-scrna naming`** — the original
  assignment referenced `seurat-scrna-v2` as the per-sample
  recipe, but the catalogue actually has `seurat-scrna`. Confirm
  whether v2 was renamed, deprecated, or never existed; the
  cross-ref currently points at `seurat-scrna` with a REVIEW
  marker.

## Approval criteria

A reviewer approves when:
- Every High priority marker (1-5) is answered or removed —
  these change the framework's recommendations.
- Every Medium priority marker (6-9) is answered, accepted as
  is, or softened to a hedge.
- The TODO(recipe) gap is either filled (a recipe is authored)
  or annotated as "user runs inline snippets" with an issue
  tracker reference.
- The TODO(knowhow) sibling is confirmed in flight or this
  knowhow's See-also is rewritten to not forward-ref a
  non-existent knowhow.
- The `seurat-scrna` vs `seurat-scrna-v2` discrepancy is
  resolved.
- The reviewer's name + date are added to the frontmatter
  (`reviewed_by`, `reviewed_on`).
- The frontmatter `kind` is changed from `knowhow_draft` to
  `knowhow`.
