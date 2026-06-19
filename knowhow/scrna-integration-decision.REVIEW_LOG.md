# Review queue — scrna-integration-decision

**Draft status:** AUTHORED (awaiting expert review)
**Author:** Claude (Phase δ Worker 1)
**Drafted on:** 2026-06-05
**Estimated review time:** ~45 minutes (5 inline REVIEW markers + 1 TODO + 5 recipe-coverage gaps; method-comparative judgment required for the High-priority items)

## High priority — claims that change the recommendation if wrong

1. **`§Quick-decision — Harmony as default for the typical case`** —
   `# REVIEW(integration-methods)`: the table makes Harmony the default
   for "typical" datasets and routes scVI/scANVI to atlas-scale + complex
   batch + label-aware tasks. This mirrors the sc-best-practices [8] /
   Luecken 2022 [1] consensus (Harmony "simple", scVI/scANVI "complex"),
   but is THIS still the right tradeoff for an ABA biologist user in
   2026, given that scVI is now broadly available and ABA already has a
   working `scvi-integration` recipe with GPU pinning in place? An
   expert should confirm: keep Harmony as default, or promote scVI to
   default when a GPU is available. Affects: §2 entire table + §4
   "default" framing.

2. **`§4 — CCA still in the matrix?`** —
   `# REVIEW(integration-methods)`: CCA is included with a "rarely the
   right *first* recommendation" hedge but kept in the matrix because
   Seurat v5 defaults to it and many published papers use it. Should
   the row be MORE strongly demoted ("present for awareness; do not
   pick"), kept as written (informative but discouraged), or remain
   neutral? An expert should set the temperature of the CCA row.
   Affects: §4 CCA row + the §5 anti-pattern "Choosing CCA because
   Seurat tutorials use it".

3. **`§3 — Harmony with high theta for conservative integration`** —
   `# REVIEW(integration-methods)`: the "use Harmony with high theta to
   be conservative when cell types are sample-specific" claim is
   community lore (the harmonypy parameter is well-known). Is there a
   citation we can pin this to (Korsunsky lab follow-up, an issue
   thread, a published worked example), OR should it be softened to
   "Harmony parameters can be tuned for conservatism — consult the
   harmonypy docs"? Affects: §3 third bullet.

## Medium priority — claims that change the explanation if wrong

4. **`§4 — scGen omission`** —
   `# REVIEW(integration-methods)`: scGen is in the Luecken 2022 [1]
   top tier on complex tasks alongside scANVI/scVI/Scanorama, but is
   omitted from the matrix on the rationale that "for most users the
   relevant label-aware pick is scANVI, which we recommend instead."
   Is the omission defensible (scGen requires a perturbation framing
   that most users don't have), or should it be included for
   completeness with a "rarely the right pick for vanilla integration"
   note? Affects: §4 matrix completeness.

5. **`§4 — FastMNN row, is it worth keeping?`** —
   `# REVIEW(integration-methods)`: FastMNN/batchelor is in the matrix
   for completeness but flagged as superseded. Either keep it where it
   is (matrix completeness), demote it to §7 (mention-only), or remove
   it entirely. The "no dedicated ABA recipe" gap should drive the
   decision: if we won't add a FastMNN recipe, the row should not
   stay in the alternatives matrix. Affects: §4 FastMNN row + recipe
   coverage gap below.

## Low priority — minor wording / cosmetic

6. **`§7 — Tran 2020 [2] citation`** —
   `# REVIEW(citations)`: I could not retrieve the Tran 2020 abstract
   via WebFetch (Springer auth wall). Inline citations to it in §4 are
   generic / community-knowledge ("Harmony performance, scaling,
   limitations"). An expert should either confirm those claims hold
   per the Tran 2020 text, soften the inline citation, or remove it
   (the document is well-supported by Luecken 2022 [1] alone).

## Coverage gaps surfaced during authoring

These are recipe-catalogue gaps where the knowhow recommends or mentions
a method that does not have a dedicated executable recipe in
`aba-recipes/recipes/genomics/`. The user is currently routed via
`bp-data-integration` (multi-method) or has no recipe path at all.

- `# TODO(recipe): scanvi-integration` — scANVI is recommended in §2
  and §4 (label-aware integration; Luecken 2022 [1] top tier). It is
  currently reached via `scvi-integration` then a scANVI step inside
  `bp-data-integration`. Worth promoting to its own recipe given how
  often it's the *correct* method when labels exist.
- `# TODO(recipe): bbknn-integration` — BBKNN (Polanski 2020 [6]) is
  in §4 as the fastest atlas-scale option. Only reachable today via
  `bp-data-integration`. A dedicated single-method recipe would clean
  up the "you don't need an embedding, just a graph" path.
- `# TODO(recipe): scanorama-integration` — Scanorama (Hie 2019 [7]) is
  in §4 as a top-tier complex-task method. Same situation as BBKNN.
- `# TODO(recipe): fastmnn-integration` — FastMNN / batchelor is in
  §4 but no recipe. Tied to High-priority item #5 (whether to keep the
  row at all).
- `# TODO(knowhow): scrna-de-methodology` — referenced in §3 (don't
  integrate over the biological variable) and §5 (per-cell vs
  pseudobulk DE) as the natural downstream knowhow. The memory item
  `scrna_de_method` covers this in ABA's project memory; a dedicated
  knowhow doc would close the loop.

LIGER is intentionally omitted from this gap list — the literature
does not push us toward adding a LIGER recipe given Luecken 2022 [1]
rankings.

## Approval criteria

A reviewer approves when:
- Every High priority marker (1–3) is answered or removed — the
  default-method tradeoff (Harmony vs scVI), the CCA temperature, and
  the high-theta-Harmony claim are the three load-bearing decisions in
  this doc.
- Every Medium priority marker (4–5) is answered, accepted, or
  softened.
- The Low priority Tran 2020 citation is either confirmed or removed.
- TODO(recipe) gaps are tracked — the knowhow can ship without these
  recipes existing (they are flagged in §7 with "no dedicated recipe
  yet" notes), but the author of the next sweep should pick at least
  scanvi-integration up.
- The reviewer's name + date are added to the frontmatter (`reviewed_by`,
  `reviewed_on`).
- The frontmatter `kind` is changed from `knowhow_draft` to `knowhow`.

## Notes for the reviewer

- The Quick-decision table (§2) is the most-read section. If you change
  one thing in this doc, change it there first — most users will stop
  reading at the bottom of §2.
- §4's matrix is intentionally one-sentence-per-cell. If a cell needs
  expanding, the meta-skill (`aba-knowhow-authoring`) says factor into
  `references/<method>.md`; do not let the matrix grow into prose.
- The doc takes a clear position (Harmony first, scVI/scANVI when
  atlas-scale or labels-present, CCA discouraged) — please don't soften
  it to "it depends on your data" without naming the dataset features
  that drive the decision. Section §2's table is the codification of
  "what it depends on, specifically".
