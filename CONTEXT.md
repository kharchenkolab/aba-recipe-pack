# The ABA standard bioinformatics cookbook

This content layer carries the broad bioinformatics + analytical
recipes that the ABA platform itself does not ship. The platform
provides the runtime and a small set of operating skills; this
layer provides the biology: scRNA-seq workflows (scanpy, Seurat,
scVI, harmony, pagoda2), bulk RNA-seq DE (DESeq2, limma-voom),
variant calling, alignment, chromatin (ATAC, ChIP-seq), database
fetches (GEO, ENCODE, Ensembl, STRING), pharmacology, microbiology,
imaging, and more.

When an institution overlay sits on top of this layer, it inherits
all of these as defaults and adds/overrides only the institution-
specific pieces.

## Where to look for things

- `recipes/<domain>/` — hand-curated recipes (genomics, database,
  genetics, molecular_biology, support_tools).
- `recipes/biomni-derived/<domain>/` — recipes distilled from upstream
  [Biomni](https://github.com/snap-stanford/Biomni) (pharmacology,
  microbiology, pathology, immunology, biochemistry, …). Implementations
  use ABA's own libraries; Biomni is referenced via `source:` provenance
  only, never imported at runtime.
- `knowhow/` — reference text that recipes cross-link to.
- `catalog/python_bio.yaml`, `catalog/r_bioconductor.yaml` — hand-curated
  capability seed YAMLs.
- `catalog/biomni-derived/` — reference catalogue mined from Biomni's
  tool descriptions (discovery metadata; not runnable here).

## Conventions

- Recipes use **`Skill(skill=<name>)`** to load procedures.
- Plans **`present_plan(steps=[{title, description, skill?, …}])`**
  with one step per major operation. Bind `skill:` to a recipe name
  when one applies.
- Cross-references between recipes use the canonical `name:` of the
  target recipe. The `aliases:` mechanism lets a renamed recipe
  keep responding to its historical names without breaking every
  recipe that mentions it.
- Capability seeds use `ensure_capability(name)` semantics — pip /
  conda / Bioconductor provisioning.
