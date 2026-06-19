# Review queue — scrna-doublet-detection

**Draft status:** AUTHORED (awaiting expert review)
**Author:** Claude (Phase δ Worker 4, ABA knowhow auto-gen)
**Drafted on:** 2026-06-05
**Estimated review time:** ~45 minutes (5 inline `# REVIEW(...)` markers + 4 `# TODO(recipe):` flags + the central Xi-vs-Germain tension)

## High priority — claims that change the recommendation if wrong

1. **`§4 closing paragraph "Xi & Li vs Germain tension"`** —
   `# REVIEW(doublet-detection-methods): the Xi & Li vs Germain tension is the
   single most important judgment in this doc. Confirm the current (2026)
   consensus is still "scDblFinder > DoubletFinder", or surface a more recent
   benchmark.`
   Affects: §1 Quick Decision (the default recommendation), §4 alternatives
   matrix (the row ordering), and the entire framing of the doc. If a 2024-2026
   benchmark has flipped the ordering — or if a newer caller (e.g. vaeda) has
   overtaken scDblFinder — the top-of-doc recommendation must change.

2. **`§3 "The doublets ARE the biology"`** —
   `# REVIEW(doublet-detection-methods): cite a representative PIC-seq /
   designed-doublet paper here so the reader has a concrete anchor (Giladi
   2020?).`
   Affects: §3 escape-valve credibility. Without a real anchor citation, the
   "deliberate doublet recovery" carve-out is hand-wavy. A reviewer with
   single-cell-interaction experience should name 1-2 canonical papers.

3. **`§3 "Plate-based protocols (Smart-seq2/3, MARS-seq)"`** —
   `# REVIEW(doublet-detection-methods): confirm — does the field actually skip
   doublet calling for Smart-seq2? Or is it just less common?`
   Affects: §3 escape-valve. If plate-based protocols still benefit from
   computational doublet calling (some labs do run scDblFinder on Smart-seq3
   data), this carve-out becomes too aggressive and should be softened to "less
   important" rather than "skip".

## Medium priority — claims that change the explanation if wrong

4. **`§5 Anti-patterns "Ensembling 5+ callers naively"`** —
   `# REVIEW(doublet-detection-methods): is there a published consensus
   ensemble strategy beyond Germain's "average two methods"?`
   Affects: §5 anti-pattern severity. If a more sophisticated published
   ensemble exists (e.g. weighted by per-sample dbr, or rank-based across
   3+ methods), the anti-pattern should be softened.

5. **`§4 "10x Cell Ranger built-in"` row** —
   The row claims standard 3'/5' v3.1 Cell Ranger does NOT emit per-droplet
   doublet calls — only the expected multiplet rate from loading. Verify
   against current Cell Ranger releases (2024-2026). Specifically check
   whether Cell Ranger 8+ or GEM-X added a per-cell doublet column.
   `# REVIEW(doublet-detection-methods, 10x-genomics): does current Cell
   Ranger emit per-cell doublet flags for standard 3'/5' chemistry?`

## Low priority — minor wording / cosmetic

6. **`§2 10x multiplet table`** — figures sourced from CG000315 Rev E user guide
   via web search summary, not direct PDF parse (the PDF binary content
   couldn't be extracted by WebFetch). Spot-check 2-3 rows against the actual
   user-guide PDF table before promotion.

7. **`§7 citation list — Heumos 2023 venue`** — the Single-cell Best Practices
   book has been variously cited as Nat Rev Genet (the companion review) and
   the sc-best-practices.org online book. The two are not the same document;
   our citation references the online book (correct for the doublet section)
   but the Nat Rev Genet tag could mislead.

## Coverage gaps surfaced during authoring

- `# TODO(recipe): scdblfinder-r` — recommended as the default in §1 for R
  users; currently only mentioned inline in `bp-quality-control`. Authoring a
  dedicated recipe would let the agent bind §1's recommendation to a real
  executable plan step. **Highest-priority recipe gap.**
- `# TODO(recipe): doubletfinder-r` — recommended in §4 for ensemble use and
  for users in legacy Seurat pipelines. No standalone recipe in the catalogue.
- `# TODO(recipe): solo-doublet-scvi` — recommended in §1 / §4 for users who
  already trained an scVI model. No standalone recipe.
- `# TODO(recipe): htodemux-r` (or a future demultiplex recipe) — for
  hashing-based confirmed-doublet calling. The §4 row points at a nonexistent
  recipe.
- `# TODO(recipe): scds-doublet-r` — covered for completeness in §4 but no
  recipe; low priority since scds is rarely the default in modern pipelines.

## Approval criteria

A reviewer approves when:
- All three High-priority markers (#1, #2, #3) are answered or removed.
- Medium-priority markers (#4, #5) are answered, accepted, or softened.
- At least the **scdblfinder-r** recipe TODO is either filled or annotated
  with a "no recipe yet; use the R snippet in `bp-quality-control`" pointer.
- 10x multiplet table figures (§2) are spot-checked against the actual
  CG000315 user-guide PDF.
- `reviewed_by` + `reviewed_on` added to frontmatter.
- Frontmatter `kind` changed from `knowhow_draft` to `knowhow`.

## Author notes for the reviewer

- The doc takes a strong position: **scDblFinder is the default for R,
  scrublet for Python.** This contradicts the older Xi & Li 2021 ranking
  (which placed DoubletFinder #1) but matches the Germain 2021 +
  Best Practices book consensus. If the reviewer's expertise diverges,
  flipping the default is the single largest edit that will be needed.
- The §3 ("when NOT to do this") section is intentionally aggressive —
  doublet calling is over-applied in practice, and the assignment briefed
  this section as especially important. If the reviewer thinks §3 is too
  permissive about skipping the analysis, soften the bullets but keep the
  section's structural prominence.
- §5 anti-patterns are ordered by expected user-impact, not by literature
  novelty. The "before strict QC, after empty droplets" ordering rule is
  the most-cited pitfall in the scDblFinder vignette and gets the longest
  treatment.
- All citations have been WebFetch'd or WebSearch'd; no abstract-only
  citations. Spot-check ~10% per `review_process.md`.
