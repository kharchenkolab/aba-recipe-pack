# Review queue — scrna-de-methodology

**Draft status:** AUTHORED (awaiting expert review)
**Author:** Phase δ Worker 2 (Claude / aba-knowhow-authoring)
**Drafted on:** 2026-06-05
**Estimated review time:** ~45 min (8 substantive REVIEW markers + 3 recipe-gap TODOs)

The central claim — "for cross-sample condition DE in scRNA-seq, use
pseudobulk + DESeq2/edgeR; never bulk methods on per-cell counts" — is
settled by Squair 2021 and re-confirmed by Crowell 2020 and the Heumos
2023 best-practices book. This draft takes that position firmly in §2
and §5. Most REVIEW markers are about secondary magnitudes (the "10-100×"
inflation figure) or about whether catalogue recipes match the
recommended paths.

## High priority — claims that change the recommendation if wrong

1. **`§Quick-decision hard-rule paragraph` + `§Sanity-checks "thousands of
   significant genes"`** —
   `# REVIEW(scrna-de): is "10-100×" the right magnitude to quote? Squair's
   figures show inflation depends on the cells-per-subject ratio — give the
   range with the regime that produces each end.`
   Affects: the central magnitude claim repeated in §2, §5, §6. If the
   actual range in Squair 2021 is narrower (e.g. "3-30×") or wider, the
   sanity-check rule of thumb ("5000+ hits = inflation") needs adjustment.
   Reviewer action: read Squair Figure 2 / S2 and report the range +
   conditions; update both §2 and §6.

2. **`§3 "One sample, no replicates, but you want a condition effect"`** —
   `# REVIEW(scrna-de): is this the right framing for the n=1-per-arm case?
   Some labs use MAST with a subject random effect even at n=1; Squair
   argues this still inflates FPR.`
   Affects: §3 escape-valve recommendation. If MAST+random-effect at n=1
   per arm is community-accepted as a "best you can do" option, §3 should
   say so explicitly rather than rejecting the comparison outright.

3. **`§4 row "MAST"` + the broader Q2 vs Q3 split** —
   `# REVIEW(scrna-de): MAST with a subject random effect is listed as
   "acceptable for Q3". Heumos 2023 says so; Squair 2021 is more
   skeptical. Reconcile — does the knowhow recommend pseudobulk DESeq2
   as the ONLY good Q3 path, or pseudobulk-first + MAST-with-RE as a
   secondary?`
   Affects: §4 Q3 recommendations + §5 anti-patterns. Current draft makes
   pseudobulk the strict pick; if MAST+RE is also acceptable, the
   "Q3-can-only-use-pseudobulk" framing softens.

## Medium priority — claims that change the explanation if wrong

4. **`§5 "Pseudobulk with too few cells"` — the 30-cell threshold** —
   `# REVIEW(scrna-de): is 30 cells the current consensus threshold,
   or has it tightened?`
   Affects: §5 anti-pattern threshold + a parameter recommendation
   downstream users will copy. Heumos 2023 uses 30; some recent papers
   recommend higher (50-100). Reviewer action: pick a single number and
   cite it.

5. **`§6 "Stability across resampling" — 70%/20% overlap rule`** —
   `# REVIEW(scrna-de): do you have a citation for the 70%/20% rule of
   thumb, or is this community lore?`
   Affects: a specific quantitative sanity-check claim. If uncited,
   either find a citation, soften to "the pseudobulk result should be
   notably more stable", or remove.

6. **`§7 Crowell 2020 citation DOI`** —
   `# REVIEW(scrna-de): verify the DOI; reading from memory.`
   Affects: a load-bearing supporting citation. Quick to verify.

7. **`§7 Soneson & Robinson 2018 citation DOI + framing`** —
   `# REVIEW(scrna-de): confirm DOI; per-cell methods looked better here
   than under Squair 2021, but the ground truth set was different.`
   Affects: explanation of why the field reversed position between
   Soneson 2018 and Squair 2021. If the reviewer can articulate it in
   one sentence, add it to §7 — improves the credibility of the central
   position.

## Low priority — minor wording / coverage

8. **`§7 Adjacent knowhow — scrna-integration-knowhow`** —
   `# REVIEW(coverage): is the integration knowhow drafted yet?`
   Affects: cross-link existence. If the integration knowhow isn't
   drafted, change the bullet to a `# TODO(knowhow):` flag.

## Coverage gaps surfaced during authoring

These are recipe gaps named in §7 and inline in §2 / §4. Each leaves an
agent following the knowhow without an executable path for some branch.

- **`# TODO(recipe): seurat-de-testing`** — the prompt assumed this
  recipe exists ("Path A / Path B / Path C in `seurat-de-testing`"); it
  does not. Currently `seurat-scrna` Step 7 covers Path A (Wilcoxon)
  and exposes `test.use` for Paths B (MAST/LR), but does not factor
  out Path C (Seurat v5 `AggregateExpression` → DESeq2 pseudobulk).
  **Decision needed:** either (a) author `seurat-de-testing` as a
  standalone recipe consolidating Paths A/B/C, or (b) update §2 / §4
  pointers to clarify "Q3 path needs the cross-language handoff:
  pseudobulk in Seurat, then `deseq2-r`". Option (a) is cleaner for
  users; option (b) is cheaper now.
- **`# TODO(recipe): edgeR-r`** — §4 lists edgeR but the catalogue has
  DESeq2 (`deseq2-r`) + limma (`limma_voom`) + pydeseq2 (`bulk-rnaseq-de`)
  only. edgeR is a one-line call inside `run_r`, but a dedicated recipe
  would close the §4 row.
- **`# TODO(recipe): scrna-pseudobulk-de-r`** — the R equivalent of
  `bp-differential-expression` (the Python pseudobulk pipeline). A
  natural companion for Seurat users; would make Q3 entirely R-resident.

## Approval criteria

A reviewer approves when:

- Markers 1, 2, 3 (High priority) are answered or removed — these
  change the load-bearing recommendation.
- Markers 4, 5, 6, 7 (Medium) are answered, accepted, or softened to a
  hedge.
- The `# TODO(recipe): seurat-de-testing` gap has a decision: either
  the recipe is authored, or §2 / §4 are updated to use the cross-recipe
  handoff explicitly without naming a recipe that doesn't exist.
- The `# REVIEW(...)` markers in the body are removed (the reviewer
  edits the prose to the verified claim).
- Frontmatter `kind: knowhow_draft` → `knowhow`; `reviewed_by` +
  `reviewed_on` added.

## Notes for the reviewer

- **The position is intentionally strong.** Squair 2021 is the textbook
  result for this question — the knowhow doesn't hedge on it. If the
  reviewer disagrees, the right move is to start over (per
  `review_process.md` — when framing is wrong, redraft rather than
  patch). Most editing should be magnitude precision, recipe-pointer
  cleanup, and MAST-with-RE nuance for Q3.
- **Q1 vs Q2 vs Q3 framing is load-bearing.** Almost every common
  user mistake reduces to "applied the answer for question X to
  question Y". If the reviewer wants to change this framing, expect
  cascade edits through §2, §4, §5.
- **The anti-pattern list is intentionally redundant on the core rule.**
  Four of five anti-patterns are restatements of "bulk methods on
  per-cell counts is wrong, in different languages and frameworks".
  This is deliberate — the rule is the most-violated pattern in the
  literature and the user prompt explicitly asked for emphasis in §2
  AND §5.
