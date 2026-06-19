# Review queue — scrna-analysis (top-level hub)

**Draft status:** AUTHORED (awaiting expert review)
**Author:** Claude (orchestrator, this session)
**Drafted on:** 2026-06-05
**Estimated review time:** ~30 minutes (small — most claims are structural, not method-comparative)

## High priority — claims that change the hub's recommendation if wrong

1. **Stage breakdown completeness** — `# REVIEW(scrna-overall-framing)` —
   the stage table covers: load+QC, sample structure, doublet detection,
   normalization, integration, clustering, annotation, marker discovery,
   DE, multimodal, trajectory, cell-cell. Verify nothing critical is
   missing (specifically: ambient RNA correction via SoupX/CellBender,
   cell-cycle regression, embedded denoising methods like MAGIC). If
   these belong, add as stages 1.5 / 4.5 / etc.

2. **Software ecosystem framing** — "Seurat for R-native and Bioconductor-
   adjacent projects; scanpy for Python/atlas-scale" is a working
   heuristic, not citation-backed. An expert may want to soften this
   or split by sub-axis (project scale, downstream tooling).

## Medium priority — claims that change the wording if wrong

3. **Multi-modal scope** — currently routes spatial to "out of catalogue
   scope right now". Confirm: is spatial coming back, in which case the
   wording should hint at "later", or permanently out of scope?

4. **Cross-cutting anti-patterns list** — 8 anti-patterns. Verify the
   list is comprehensive at the pipeline level (vs stage-specific
   anti-patterns which belong in per-stage knowhow). Specific candidates
   to add or remove based on expert judgment.

## Low priority — cosmetic / structural

5. **ASCII pipeline diagram** — included as a quick visual. Verify it
   conveys the right mental model; replace with a different shape if
   not (e.g. a decision tree, a sequential flow with branch points).

6. **Sub-knowhow naming** — the names referenced (`scrna-de-methodology`,
   `scrna-integration-decision`, etc.) are the names being drafted in
   parallel by the δ workers. If any final name differs from what's
   written here, update the cross-references.

## Coverage gaps surfaced during authoring

- `# TODO(knowhow): scrna-trajectory-analysis` — stage 12 has no knowhow doc
- `# TODO(knowhow): scrna-cell-cell-communication` — stage 13 has no knowhow doc
- `# TODO(knowhow): scrna-ambient-rna` (if SoupX/CellBender is added as a stage)
- `# TODO(recipe): standalone doublet-detection recipe` (currently inline)
- `# TODO(recipe): SingleR-annotation recipe` (currently only Azimuth via multimodal-ref-mapping)

## Approval criteria

A reviewer approves when:
- Stage breakdown verified or amended
- Software-ecosystem framing accepted or refined
- Spatial scope clarified
- Frontmatter `kind: knowhow_draft` changed to `knowhow`, `reviewed_by` + `reviewed_on` added
