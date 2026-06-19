# Review queue — scrna-celltype-annotation

**Draft status:** AUTHORED (awaiting expert review)
**Author:** Phase δ Worker 5 (Claude Opus 4.7, ABA knowhow-authoring skill)
**Drafted on:** 2026-06-05
**Estimated review time:** 60–90 min (six REVIEW markers + three TODO(recipe) gaps + cross-check against current 2024–2025 annotation literature)

## High priority — claims that change the recommendation if wrong

1. **`§4 SingleR row — "no dedicated singler-celltype-annotation recipe exists"`** —
   `# REVIEW(annotation-methods): confirm there is no SingleR-specific recipe in ABA's catalogue. seurat-scrna-v2 mentions SingleR but does not execute it; bp-annotation skips SingleR in favor of CellTypist/scArches.`
   Affects: §2 Quick-decision row pointing R/Bioconductor users at SingleR — currently has no executable recipe. Either author a `singler-celltype-annotation` recipe (recommended; SingleR is widely used in the R single-cell community and the gap is real) or down-weight the SingleR recommendation in §2.

2. **`§2 Quick-decision row "ensemble that hedges between methods" → popV`** —
   `# REVIEW(annotation-methods): is popV the right ensemble recommendation in 2025, or has the field moved on (e.g. CZ CELLxGENE's reference annotation, or scvi-tools' integrated annotation framework)? unsupervised-celltype-transfer-between-scrna-datasets is biomni-distilled and may be stale.`
   Affects: whether the ensemble row should stay in the quick-decision table at all.

3. **`§5 Anti-pattern "over-trusting Azimuth confidence scores"`** —
   `# REVIEW(annotation-methods): is Azimuth's mapping.score + predicted.celltype.l2.score considered systematically miscalibrated on out-of-distribution queries, or is this folklore? Looking for a specific paper or Satija-lab note documenting the calibration failure mode.`
   Affects: anti-pattern strength. If miscalibration is documented, cite the paper; if folklore, soften to "validate against markers regardless of score" without the miscalibration claim.

## Medium priority — claims that change the explanation if wrong

4. **`§4 Manual marker panels row — PanglaoDB vs CellMarker 2.0 status`** —
   `# REVIEW(annotation-methods): is the Franzén 2019 PanglaoDB marker set still considered the canonical reference for broad lineages, or has CellMarker 2.0 (Hu 2023) overtaken it for breadth? Both are cited in §7 but the §4 row treats PanglaoDB as primary.`
   Affects: which database the manual-annotation row recommends as the default.

5. **`§5 Anti-pattern "brain reference on intestinal data"`** —
   `# REVIEW(annotation-methods): is there a named tissue-mismatch failure paper to cite here, or is this community consensus only? If consensus, mark as such; if there's a paper documenting this failure mode systematically, cite it.`
   Affects: citation discipline (every claim cited or hedged).

6. **`§4 GeneCT row — maintenance status`** —
   `# REVIEW(annotation-methods): is GeneCT actively maintained as of 2025? If not, drop from §4 — the user-task framing only mentioned it as a narrower alternative to mention, but listing an unmaintained tool as an option is misleading.`
   Affects: whether GeneCT (and the corresponding `# TODO(recipe): scgate / GeneCT`) stays in scope.

## Low priority — minor wording / cosmetic

7. **`§7 scANVI citation`** —
   `# REVIEW(citations): confirm Xu 2021 (Mol Syst Biol) is the correct primary citation for scANVI vs the scvi-tools framework paper (Gayoso 2022, Nat Biotechnol). The scvi-tools docs cite both; pick the one the methods chapter of a paper would normally cite.`
   Affects: which paper is named as the scANVI primary source.

8. **`§7 Benchmarks — Pasquini 2021 / Abdelaal 2019 recency`** —
   `# REVIEW(annotation-methods): is there a 2023–2025 benchmark that has superseded Pasquini 2021 / Abdelaal 2019 for annotation method choice? Older benchmarks may pre-date CellTypist and scArches updates.`
   Affects: which benchmark anchors the alternatives matrix's claims.

## Coverage gaps surfaced during authoring

These are method recommendations made in §2 / §4 that DO NOT currently have an executable recipe. Each is flagged in §7 of the draft with `# TODO(recipe):`. Decision options for each: (a) author the recipe and remove the TODO, (b) remove the recommendation from §2 / §4 until a recipe exists, (c) keep with explicit "no recipe yet; use upstream docs" annotation.

- **`# TODO(recipe): singler-celltype-annotation`** — SingleR is recommended for R/Bioconductor users with bulk-sorted references (HumanPrimaryCellAtlas, MonacoImmune, ImmGen, BlueprintEncode via the `celldex` package). The §2 quick-decision row points here. Current state: ad-hoc only via `run_r`. **Recommendation: author this recipe — it is the largest gap.**

- **`# TODO(recipe): azimuth-seurat-direct`** — The Seurat-native Azimuth path (`Seurat::MapQuery` against Satija-lab tissue-specific references such as PBMC / lung / motor cortex / heart / bone marrow / pancreas / kidney / tonsil) currently has no ABA recipe. `annotate-celltype-with-panhumanpy` covers the pan-human ANN but not the original tissue-specific Azimuth references. **Recommendation: author as a second R/Seurat annotation recipe alongside the SingleR one above.**

- **`# TODO(recipe): scgate-purify`** — Narrow alternative for population purification (e.g. "purify T cells from a mixed object"). Mentioned in §4 as a narrower alternative; not a quick-decision row. **Recommendation: low priority; only author if user feedback drives it.**

- **`# TODO(recipe): genect-annotation`** — Narrower alternative; flagged for §4 GeneCT row removal pending REVIEW marker #6 above.

## Approval criteria

A reviewer approves when:
- Every High priority marker (1–3) is answered or removed. In particular, marker #1 either spawns the `singler-celltype-annotation` recipe or downgrades the SingleR recommendation in §2.
- Every Medium priority marker (4–6) is answered, accepted as is, or softened to a hedge with named citation or "community consensus".
- TODO(recipe) gaps are filled or annotated with explicit "no recipe yet — use <upstream-docs-URL>" notes in the draft, so the agent isn't stranded mid-plan when a user picks one of those rows.
- The reviewer's name + date are added to the frontmatter (`reviewed_by`, `reviewed_on`).
- The frontmatter `kind` is changed from `knowhow_draft` to `knowhow`.

## Notes for the reviewer

- The Quick-decision table (§2) is the most-read section. If you change a §2 row, propagate the change to §4 and §5 — they were authored to be internally consistent.
- The "always do manual marker validation" framing is intentionally repeated across §2, §3, §4 (LLM-from-markers row), and §5 (anti-patterns). This is the load-bearing position the draft takes — if you want a different position, that's a framing change, not a wording fix.
- §3 (When NOT to do automated annotation) is intentionally short and high-leverage; resist the temptation to add more edge cases. The four bullets cover the cases users actually attempt; everything else belongs in §5 (anti-patterns) or in the future `scrna-best-practices` knowhow once one exists.
- Citations were drawn from the user-task brief + cross-checked against ABA's existing `bp-annotation` recipe (which itself cites the Single-cell Best Practices chapter). No paper was read end-to-end during authoring — every cited claim is a candidate for spot-check, not just the explicitly REVIEW-marked ones.
