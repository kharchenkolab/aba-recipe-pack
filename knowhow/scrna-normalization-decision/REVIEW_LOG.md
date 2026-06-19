# Review queue — scrna-normalization-decision

**Draft status:** AUTHORED (awaiting expert review)
**Author:** Phase δ Worker 3 (ABA agent)
**Drafted on:** 2026-06-05
**Estimated review time:** ~60–90 min (7 REVIEW markers; ~3 citation spot-checks; 2 recipe-gap TODOs)

This knowhow is explicitly MEDIUM-CONFIDENCE. The literature on scRNA-seq
normalization is less settled than on DE-method choice — the field's most
defensible empirical claim (Booeshaghi 2022) is that depth normalization
matters MORE than the choice among shifted-log / SCT / Pearson residuals,
so several positions in this doc are framed as "regime-specific" rather
than "X dominates Y". The high-priority markers below ask the expert to
either harden those hedges or confirm them.

## High priority — claims that change the recommendation if wrong

1. **`§2 Quick-decision row "Low UMI + heterogeneous → SCTransform v2"`** —
   `# REVIEW(scrna-normalization): is the "median UMI ≲2k AND heterogeneous
   per-cell-type depth" precondition the right cutoff?` The threshold and
   heterogeneity definition come from a paraphrase of Choudhary & Satija
   2022, not their literal Figure 2 / Table 1. If the precondition is
   actually broader (e.g. SCT v2 wins at median UMI ≲5k), the matrix's
   default-to-shifted-log recommendation for the typical PBMC dataset is
   wrong. Pull the paper's quantitative claims.
   *Affects:* the TL;DR table — the highest-leverage content in the doc.

2. **`§2 framing "shifted-log is the safe default"`** —
   `# REVIEW(scrna-normalization): is shifted-log still the field
   consensus default in 2025, or has analytic Pearson residuals overtaken
   it for first-pass clustering on modern droplet data?` Lause 2021 +
   Heumos 2023 are favourable to Pearson residuals but in specific
   contexts (HVG selection, rare cells). If recent benchmarks have shown
   Pearson residuals as the broader default, the matrix's top row needs
   to flip. Cite the most recent benchmark.
   *Affects:* the default recommendation. If wrong, every "use shifted-log
   unless..." line in this doc inverts.

3. **`§4 atlas-scale row (absent)`** —
   `# REVIEW(scrna-normalization): for atlas-scale data (>1M cells), is
   there a community consensus on normalization?` The doc currently does
   not name a recommendation for >1M cell datasets. Some atlas papers
   feed raw counts to scVI and skip normalization for the embedding.
   *Affects:* the doc is silent on a regime where users have real
   choices to make. Either add a row to §4 or harden the §2 row
   "Downstream is scVI" into "for atlas-scale, prefer scVI on raw counts;
   normalization happens only for plotting/marker DE".

## Medium priority — claims that change the explanation if wrong

4. **`§5 anti-pattern "Pearson residuals as the expression layer"`** —
   `# REVIEW(scrna-normalization): is this a methods claim or just an
   idiom mismatch?` Currently framed as an anti-pattern, but it may be
   purely a UX/idiom issue (tools expect log-norm). If no benchmark shows
   interpretive harm from plotting on residuals, soften to "convention,
   not error".
   *Affects:* anti-pattern severity wording.

5. **`§6 sanity check "PC1 vs depth |r| > 0.4"`** —
   `# REVIEW(scrna-normalization): is the 0.4 threshold sourced?` It is
   a practitioner rule-of-thumb without a clear citation; either cite it
   or soften to "any visibly strong correlation".
   *Affects:* a single sanity-check threshold; reviewer's call to keep,
   replace, or soften.

## Low priority — citation hygiene

6. **`§7 citation list — Choudhary & Satija 2022`** —
   `# REVIEW(citations): confirm this IS the SCT v2 paper.` The title
   was paraphrased from memory. Cross-check the title + DOI; the
   `s13059-021-02584-9` DOI may instead be the Satija-lab error-model
   paper, not the SCT v2 paper itself. If they're separate, list both.

7. **`§7 citation list — Booeshaghi et al. 2022`** —
   `# REVIEW(citations): venue + final DOI.` Listed without DOI; need to
   confirm whether this is a published paper or a still-cited preprint,
   and that the headline claim ("depth normalization is the most
   important step") appears in the abstract / conclusion.

## Coverage gaps surfaced during authoring

- `# TODO(recipe): seurat-sctransform` — the assignment named a
  standalone `seurat-sctransform` recipe. The catalogue currently carries
  SCT only as an alternative block inside `seurat-scrna-v2`. Two
  defensible resolutions:
  1. Split the SCT path into its own recipe and re-point §2 / §4 at it.
  2. Leave SCT inlined and keep the current pointer.
  Reviewer to decide; the matrix points at the current state.

- `# TODO(recipe): scanpy-pearson-residuals-clustering` — Analytic
  Pearson residuals are accessible through `bp-normalization` (a stage,
  not a full pipeline). A user wanting "scrna-qc-clustering, but with
  Pearson residuals for HVG/PCA" must compose two recipes manually. A
  combined recipe would close this gap and let §2 row 3 point at a
  single executable target instead of "compose these two".

## Approval criteria

A reviewer approves when:
- Both High-priority markers 1 and 2 are answered (these set the TL;DR
  table's defaults; cannot ship a knowhow on a wrong default).
- High-priority marker 3 (atlas-scale) is either answered or explicitly
  scoped out ("this knowhow does not cover >1M cells; see X").
- Medium-priority markers 4–5 are answered, softened, or accepted.
- Low-priority citation checks 6–7 are resolved.
- TODO(recipe) gaps are either filled or annotated as "deferred,
  knowhow points at the current state".
- Frontmatter `reviewed_by` + `reviewed_on` are filled in.
- Frontmatter `kind: knowhow_draft` is changed to `kind: knowhow`.
