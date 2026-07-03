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
source: "ABA viewer framework + pagoda3 integration (misc/pagoda3_integration.md); AnnData / .h5ad spec; lstar (L* Zarr interchange, lstar-sc on PyPI); cellxgene + CELLxGENE Census docs; pagoda2 export_and_interop reference (recipes/genomics/pagoda2-scrna-v2)."
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
- **Ask Guide** — Guide can hand you a launch link directly (the `open_viewer`
  tool). This is the fastest path when you're already in chat.

On first open ABA prepares the data store (a brief "preparing…" screen), then
loads the viewer; subsequent opens are instant (the prepared store is cached).

**For the agent:** when you finish a single-cell analysis that yields a
clustered or annotated object, *offer* to open it in pagoda3 — call
`open_viewer(entity_id=…)` and present the returned link. Only offer this for
single-cell result objects (`.h5ad` / `.lstar.zarr`); a figure or table already
renders inside ABA and needs no viewer. Don't nag — offer once, when the result
is ready.

## Interchange formats — what to save

| Format | What it is | Reach for it when |
|---|---|---|
| **`.h5ad`** (AnnData) | The single-cell lingua franca; scanpy-native, read by nearly every Python/R tool and by cellxgene | Default for sharing, archiving, and cross-tool moves. When in doubt, save `.h5ad`. |
| **`.lstar.zarr`** (lstar) | Chunked, HTTP-range-readable, disk-backed, multimodal-friendly Zarr store | You want **portability + instant viewing** (it's pagoda3's zero-conversion fast path — opens with no per-launch conversion), or the object is large and shouldn't be fully loaded into RAM. Round-trips with AnnData / Seurat / pagoda2 via lstar. |
| **RDS / Seurat** (`.rds`, `.h5Seurat`) | R-native serialization | Staying inside an R/Seurat workflow. Convert to `.h5ad` before handing to a viewer or a Python tool. |
| **loom** | Older HDF5-based format | Only for a tool that specifically requires it — otherwise prefer `.h5ad`. |

Recommendation, scoped to single-cell: **save results as `.h5ad` for
interoperability, and additionally as `.lstar.zarr` when portability or instant
viewing matters.** `.lstar.zarr` is a good interchange target on its own merits
(chunked, range-readable, multimodal) — not only because pagoda3 reads it.

## How to save / convert

From scanpy / AnnData (Python):

```python
adata.write_h5ad("result.h5ad")                 # universal
import lstar
lstar.convert_anndata("result.h5ad", "result.lstar.zarr")   # lstar-sc (in ABA's base env)
```

From Seurat or pagoda2 (R): see the export reference bundled with the pagoda2
recipe — `recipes/genomics/pagoda2-scrna-v2` → `references/export_and_interop.md`
(covers `.h5ad` export and the lstar Zarr path from R objects).

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
- `misc/pagoda3_integration.md` (ABA repo) — how the viewer + launcher are wired.
