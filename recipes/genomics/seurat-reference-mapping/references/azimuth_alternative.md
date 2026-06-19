# Azimuth — pre-built tissue references

When to load this: the user names Azimuth, asks for a pre-built human
tissue reference (PBMC, lung, kidney, heart, ...), or you want a
faster, lower-knob mapping than the manual `FindTransferAnchors` +
`MapQuery` flow.

Azimuth wraps Seurat's reference-mapping flow into a one-call
`RunAzimuth(query, reference = "<refname>")` against a curated set of
human tissue references hosted by the Satija lab. The references are
pre-validated, multi-level annotated (commonly `predicted.celltype.l1`
/ `.l2` / `.l3`), and known to produce sane mappings without parameter
tuning. Use Azimuth when a matching tissue reference exists; fall back
to the manual SKILL.md flow when none does, or when you have a
custom-built reference.

## Pre-built references (as of 2026)

The reference catalogue lives at <https://azimuth.hubmapconsortium.org>
and the `Azimuth` R package. Names callable via `RunAzimuth(reference
= "<name>")`:

| Reference name | Tissue / cohort | Hierarchy depth | Notes |
|---|---|---|---|
| `pbmcref` | Peripheral blood mononuclear cells | l1 (8) / l2 (~30) / l3 (~58) | The canonical immune reference |
| `lungref` | Human lung (LungMAP / HLCA) | l1 / l2 / l3 | Healthy human lung |
| `kidneyref` | Kidney | l1 / l2 / l3 | |
| `heartref` | Heart | l1 / l2 | |
| `humancortexref` | Human cortex / brain | l1 / l2 / l3 | |
| `tonsilref` | Tonsil | l1 / l2 | |
| `bonemarrowref` | Bone marrow | l1 / l2 | |
| `fetusref` | Human fetus | l1 / l2 | Multi-organ fetal atlas |
| `adiposeref` | Adipose tissue | l1 / l2 | |
| `mousecortexref` | Mouse cortex | l1 / l2 | Mouse, not human |
| `pancreasref` | Pancreas | l1 / l2 | |

Verify the current list with `Azimuth::AvailableData()` once the
package is loaded — Satija's catalogue updates periodically.

## Install

```r
# Azimuth is on GitHub, not CRAN.
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("Seurat", quietly = TRUE))  install.packages("Seurat")
if (!requireNamespace("BPCells", quietly = TRUE)) {
  remotes::install_github("bnprks/BPCells/r")
}
if (!requireNamespace("Azimuth", quietly = TRUE)) {
  remotes::install_github("satijalab/azimuth", ref = "master")
}
```

`BPCells` is a transitive dep; install it first if it's not already
present. The first `RunAzimuth(reference = "pbmcref")` call downloads
the reference (a few hundred MB) and caches it under
`~/.cache/R/Azimuth/`.

## The one-call flow

```r
library(Azimuth)

# Single call. Replaces this recipe's Steps 1–4 entirely.
query <- RunAzimuth(query, reference = "pbmcref")
```

After `RunAzimuth`:

- `query$predicted.celltype.l1` / `.l2` / `.l3` — the multi-level
  predicted labels (depending on the reference's hierarchy depth).
- `query$predicted.celltype.l1.score` / `.l2.score` / `.l3.score` —
  per-level confidence.
- `query[["ref.umap"]]` — query cells projected onto the reference's
  UMAP. Reduce dimensions for plotting use the same `DimPlot(query,
  reduction = "ref.umap")` pattern.
- `query[["ref.spca"]]` — projection into the reference's supervised PCA.
- `query[["impADT"]]` (PBMC ref only) — imputed antibody-derived tag
  expression from the reference's CITE-seq data.

From there, jump to the SKILL.md's Step 5 visualization — the same
`DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype.l2")`
pattern works.

## Mapping Azimuth's outputs to this recipe's outputs

| Azimuth | This recipe's manual flow |
|---|---|
| `query$predicted.celltype.l2` | `query$predicted.celltype` (after `MapQuery(refdata = list(celltype = ...))`) |
| `query$predicted.celltype.l2.score` | `query$predicted.celltype.score` |
| `query[["ref.umap"]]` | `query[["ref.umap"]]` (same name, same shape) |
| `query[["ref.spca"]]` | `query[["ref.pca"]]` (analogous; Azimuth uses supervised PCA, manual uses standard) |
| Multi-level labels (l1, l2, l3) | Single label level per `MapQuery` call (or multi-list in one call) |

The visualization step (SKILL.md Step 5) and the persistence step
(Step 6) work identically — Azimuth's outputs slot into the same
column names with the `.l2` suffix.

## When Azimuth is the right pick

Use Azimuth when ALL of:

- Your query is **human** (with the one mouse exception above) and the
  tissue matches one of the catalogue entries.
- You want **multi-level annotations** (l1/l2/l3) out of the box.
- You don't have a custom reference, and there's no reason to insist
  on one.
- A first-pass mapping is acceptable; you don't need to tune
  `k.anchor` / `k.weight` / etc.

Use the **manual** SKILL.md flow when ANY of:

- You have a **custom labeled reference** (e.g. an integrated atlas
  you built yourself, a tissue-specific reference not in Azimuth's
  catalogue).
- Your query is from a **non-human, non-mouse-cortex** organism.
- You need to map onto **specific labels** that Azimuth's hierarchy
  doesn't carry (e.g. donor-specific subtypes).
- The Azimuth reference's tissue is close-but-not-right (e.g. PBMC ref
  for a query of sorted T cells — the manual flow lets you build a
  tighter T-cell-only reference).
- You need to tune anchor or transfer parameters (k.anchor, k.weight,
  k.filter).

## Caveats

- **Score distributions differ.** Azimuth's `predicted.celltype.l2.score`
  is comparable to the manual flow's `predicted.celltype.score` but the
  exact thresholds for "high confidence" depend on the reference. The
  general guidance — ≥0.8 confident, 0.4–0.8 ambiguous, <0.4 alien — still
  applies but tune per reference.
- **The reference download is large.** First call per reference is slow.
  Subsequent calls reuse the cached reference.
- **Azimuth's CCA-based anchor finding is more aggressive.** For some
  queries (e.g. heavily perturbed by drug treatment), Azimuth may
  confidently assign labels that the manual `reduction = "pcaproject"`
  flow would flag as low-confidence. Compare both if the result seems
  too clean.
- **Network dependency.** Azimuth pulls the reference at call time;
  offline runs require pre-cached references.
- **Mismatched organism / tissue silently produces nonsense.** A heart
  query against `pbmcref` will run end-to-end and produce confidently-wrong
  labels. Always verify the reference matches the query's tissue.

## Inspecting an Azimuth reference's metadata

To see what labels a pre-built reference carries before mapping:

```r
library(Azimuth)
# Lists the multi-level labels available, e.g. "celltype.l1", "celltype.l2", "celltype.l3"
ref_meta <- Azimuth::AvailableData()
print(ref_meta)
```

Or after the first `RunAzimuth` call, inspect the cached reference
directly:

```r
ref <- LoadReference(path = file.path("~/.cache/R/Azimuth", "pbmcref"))
table(ref$map$plot.layer)
```

For deep diagnostics + per-reference label vocabularies, the upstream
docs at <https://azimuth.hubmapconsortium.org> are more current than
this reference file.
