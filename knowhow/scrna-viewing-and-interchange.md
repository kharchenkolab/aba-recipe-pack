---
name: scrna-viewing-and-interchange
description: How to interactively VIEW a single-cell (scRNA-seq) result and how to choose an interchange FORMAT for sharing, archiving, and portability. Covers pagoda3 (ABA's built-in interactive viewer for .h5ad / .lstar.zarr), the .h5ad vs .lstar.zarr vs RDS/Seurat vs loom trade-offs, how to save/convert, and how to hand an .h5ad to external platforms (cellxgene, UCSC Cell Browser) that ABA does not launch itself.
when_to_use: User has produced a single-cell result (a clustered / annotated / integrated object — an .h5ad or native store) and wants to look at it interactively (UMAP, per-gene expression, cluster labels), share it with a collaborator, archive it, or move it to another tool. Also use when an agent has just finished a single-cell analysis and should offer the user a way to explore the result, or when deciding which format to save a single-cell object in.
avoid_when: The result is a plain figure / table / PDF (those open inside ABA directly — no viewer needed). The user is still mid-processing (QC / clustering not done) — finish the analysis recipe first. The data is not single-cell (this doc is scoped to scRNA-seq objects).
invocation: interactive
kind: knowhow_draft
requires_tools: [Read]
keywords: [view single cell, scRNA-seq viewer, interactive UMAP, explore clusters, expression browser, visualize cells, pagoda3, cellxgene, UCSC cell browser, h5ad, AnnData, lstar, lstar.zarr, zarr, RDS, Seurat, loom, interchange format, export single cell, share results, portability, save anndata, open viewer]
domain: genomics
source: "ABA viewer framework + pagoda3 integration; AnnData / .h5ad spec; lstar (L* Zarr interchange, lstar-sc on PyPI); cellxgene + CELLxGENE Census docs; pagoda2 export_and_interop reference (recipes/genomics/pagoda2-scrna-v2)."
audience: both
---

# Viewing & interchange for single-cell results

Two related questions this answers: **"how do I look at this single-cell result
interactively?"** and **"what format should I save it in to share / move / keep?"**

## Viewing interactively — pagoda3 (built into ABA)

ABA ships **pagoda3** as its interactive single-cell viewer. It opens a result
as a live UMAP / embedding you can color by cluster, metadata, or per-gene
expression — in a new browser tab, over ABA's own login. It handles single-cell
results saved as **`.h5ad`** or **`.lstar.zarr`** (see formats below).

There is nothing to install or export — three ways to open one:

- **Files tab** — click the `.h5ad` / `.lstar.zarr` file; its viewer is pagoda3.
- **Dataset card** — a dataset whose file is a single-cell object shows an
  **"↗ Explore in pagoda3"** button on its focus card.
- **Ask Guide** — Guide can hand you a launch link directly (the `get_viewer_url`
  tool). This is the fastest path when you're already in chat.

On first open ABA prepares the data store (a brief "preparing…" screen), then
loads the viewer; subsequent opens are instant (the prepared store is cached).

**For the agent:** producing a clustered / annotated / integrated single-cell
object? Handing over the view link is a **closing step of the analysis itself,
not an optional follow-on.** In the **same turn** you write the result — before you
summarize or suggest a next sample — call `get_viewer_url(path="processed.h5ad")`
(a bare filename is fine — it's resolved against the project's files) or
`get_viewer_url(entity_id=…)` for a registered dataset, then present the returned
link. Do **not** wait to be asked, and do **not** let it drop because your plan's
last step was named "markers" or "save": writing the result and offering the view
are one gesture, and this offer is part of the recipe's contract (not scope creep to
skip at the plan's end). Offer this for single-cell result objects only;
a figure or table already renders inside ABA and needs no viewer. What to point
`get_viewer_url` at depends on the object:

- **`.h5ad` / `.h5mu` / `.lstar.zarr`** — hand these to the viewer as-is. If the
  analysis was in R, prefer producing one of these from the live session (highest
  fidelity) and view that.
- **Seurat / SingleCellExperiment `.rds`** — you *can* hand the `.rds` directly (ABA
  converts on launch), but that's a lower-fidelity fallback for installs without the
  R stack; if you just made the object in R, export `.h5ad`/`.lstar.zarr` in-session
  and point the viewer at that instead.
- **pagoda2 / conos** — do **not** hand the raw pagoda2/conos `.rds` to the viewer
  (lstar's converter reads only Seurat/SCE `.rds`). Point `get_viewer_url` at the
  `.h5ad` or `.lstar.zarr` the analysis writes from the live object (the pagoda2
  recipe already produces `pagoda2_processed.h5ad`; see "How to save / convert").

Don't nag — offer once, when the result is ready. If `get_viewer_url` returns
`ok:false`, relay the error; never hand out a link that didn't resolve.

## Interchange formats — what to save

| Format | What it is | Reach for it when |
|---|---|---|
| **`.h5ad`** (AnnData) | The single-cell lingua franca; scanpy-native, read by nearly every Python/R tool and by cellxgene | Default for sharing, archiving, and cross-tool moves. When in doubt, save `.h5ad`. |
| **`.lstar.zarr`** (lstar, **directory**) | Chunked, disk-backed, multimodal-friendly Zarr **directory** store | You want **portability + instant viewing** (pagoda3's zero-conversion fast path — opens with no per-launch conversion), or the object is large and shouldn't be fully loaded into RAM. **Prefer the directory form**: it loads much faster in the viewer (parallel chunk fetches) and stays updatable. Round-trips with AnnData / Seurat / pagoda2 via lstar. |
| **`.lstar.zarr.zip`** (single-file) | The same store packed into one STORED (range-readable) file | Only when you need a **single file** to move/host. Loads slower than the directory and is immutable — unpack it back to the directory to work with or view it. |
| **RDS / Seurat** (`.rds`, `.h5Seurat`) | R-native serialization | Staying inside an R/Seurat workflow. Convert to `.h5ad` before handing to a viewer or a Python tool. |
| **loom** | Older HDF5-based format | Only for a tool that specifically requires it — otherwise prefer `.h5ad`. |

Recommendation, scoped to single-cell: **save results as `.h5ad` for
interoperability, and additionally as `.lstar.zarr` when portability or instant
viewing matters.** `.lstar.zarr` is a good interchange target on its own merits
(chunked, range-readable, multimodal) — not only because pagoda3 reads it.

## How to save / convert

**Preferred pattern — emit the viewer store IN-SESSION and point `get_viewer_url` at
it.** A recipe that just produced a result should write a *viewer-optimized*
`.lstar.zarr` from the in-memory object at the end of the session, and offer
`get_viewer_url(path="…​.lstar.zarr")`. Then the link opens instantly (already
optimized — no banner, no on-launch conversion) and it's highest-fidelity (built
from the live object, using its raw counts). ABA *can* also convert an `.h5ad`/`.rds`
on the fly when the user opens one directly — but that's the **fallback** for
ad-hoc files, not the path recipes should rely on.

From scanpy / AnnData (Python) — the analog of R's `lstar_write_viewer`:

```python
adata.write_h5ad("result.h5ad")                 # universal interchange
import lstar
lstar.write(lstar.read_anndata(adata), "result.lstar.zarr", viewer=True)  # optimized store
```

`viewer=True` precomputes DE / variable genes / cell-major counts. Keep raw counts
in `adata` (`.layers['counts']` or `.raw`) so those stats use real counts; lstar-sc
≥0.1.7 falls back to log-normalized only when no raw counts are present.

From R objects, **write the `.lstar.zarr` store directly from the live object with
lstar** — it's pure R, highest-fidelity (all reductions/layers/metadata carried
across), and the uniform path for every R format:

```r
d <- lstar::read_seurat(obj);   lstar::lstar_write_viewer(d, "result.lstar.zarr")  # Seurat
d <- lstar::read_sce(sce);      lstar::lstar_write_viewer(d, "result.lstar.zarr")  # SingleCellExperiment
d <- lstar::read_pagoda2(p2);   lstar::lstar_write_viewer(d, "result.lstar.zarr")  # pagoda2
d <- lstar::write_conos(con);   lstar::lstar_write_viewer(d, "joint.lstar.zarr")   # conos (joint)
```

Use `lstar_write_viewer` (not plain `lstar_write`) for the viewer path: it
precomputes DE / variable genes / cell-major counts in R so pagoda3 opens the store
**optimized** (no "Not viewer-optimized" banner) — and does it in pure R, avoiding
the WASM prep step. It needs a clustering/grouping present on the object.

**Do NOT detour through `.h5ad` just to view an R object.** Writing `.h5ad` from R
(`zellkonverter`/`sceasy`) spins up a **basilisk/reticulate Python environment** —
slow and pointless when lstar goes straight to the store in R. Export `.h5ad` only
when the target is a *different* tool (scanpy, cellxgene). (`lstar::read_*` are
pure Rcpp — no Python. Note `library(Seurat)` itself loads the `reticulate`
namespace, but that alone spins up no Python env; the store path never calls into
Python.)

**On-the-fly `.rds` is a fallback, not the preferred path.** ABA *can* convert a
Seurat/SCE `.rds` at viewer-launch (hand the `.rds` straight to `get_viewer_url`),
but that path is a compromise for installs **without** the heavy R stack — it
reads the serialized `.rds` without a full live session and is lower-fidelity
than an in-session export. When you're already in R, export `.h5ad`/`.lstar.zarr`
and view that. For **pagoda2/conos** there is no `.rds` fallback at all — lstar's
`.rds` converter reads only Seurat/SCE, so their ingest *must* go through the R
functions above on the live object. See the export reference bundled with the
pagoda2 recipe — `recipes/genomics/pagoda2-scrna-v2` →
`references/export_and_interop.md`.

Register the saved file as a dataset so it appears as an entity (and picks up
the "Explore in pagoda3" affordance) rather than living only on disk.

## Other viewers (export, not integrated)

ABA launches pagoda3 directly; it does **not** launch these — but `.h5ad` is the
handoff format, so exporting once serves all of them:

- **cellxgene** — open the `.h5ad` in cellxgene Desktop or a hosted Explorer, or
  submit to the CELLxGENE Census for public sharing.
- **UCSC Cell Browser / Cirrocumulus** — build/point them at the exported `.h5ad`.

If a lab standardizes on one of these, the workflow is the same: produce the
result, `write_h5ad`, hand off the file.

## See also

- `scrna-analysis` — the pipeline navigation hub (this doc is the "view / share
  the result" endpoint of that pipeline).
- Recipes that produce a viewable object: `scrna-qc-clustering`, `bp-annotation`,
  `bp-clustering`, `scvi-integration`, `seurat-scrna-v2`.
- ABA's pagoda3 integration — how the viewer + launcher are wired.
