# Modality-weight diagnostics — reading the `RNA.weight` column

Deep detail behind Step 5 of `seurat-wnn-multimodal`. Loaded when the agent
needs to interpret the per-cell `RNA.weight` distribution, defend a
"WNN was worth it" claim against a critical user, or diagnose a
modality-weight pattern that argues for re-tuning `dims.list`.

The `RNA.weight` metadata column is **the** diagnostic for WNN. It
answers "for each cell, how much did the RNA reduction contribute to
its joint neighborhood — vs the ADT (or ATAC) reduction?" Values:

| `RNA.weight` value | Meaning |
|---|---|
| **1.0** | Pure RNA — neighborhood is fully determined by RNA. The ADT/ATAC modality added nothing for this cell |
| **0.5** | Balanced — both modalities contributed equally |
| **0.0** | Pure ADT/ATAC — neighborhood is fully determined by the second modality. RNA added nothing |

For each cell, `RNA.weight + ADT.weight = 1` (or `+ ATAC.weight = 1`).

## The headline diagnostic: sorted VlnPlot

```r
VlnPlot(obj, features = "RNA.weight", group.by = "seurat_clusters",
        sort = TRUE, pt.size = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed",
             colour = "red", linewidth = 0.4)
```

`sort = TRUE` orders clusters by median weight — RNA-driven on one
end, ADT/ATAC-driven on the other, balanced in the middle. The dashed
red line at 0.5 is the reference; clusters above are RNA-driven,
below are ADT/ATAC-driven.

What to look for in the resulting figure:

| Pattern | Reading |
|---|---|
| **Most clusters near 0.5, narrow spread (0.4–0.6)** | The two modalities are roughly equivalent informants. WNN didn't change much vs averaging. If the joint UMAP also looks like the RNA UMAP, WNN added little value here — consider whether the panel/dims justify the joint workup |
| **Clear bimodal — some clusters at ~0.2, others at ~0.8** | WNN doing real work. RNA-driven and ADT-driven populations are biologically distinct (e.g. progenitors vs T cells in bone marrow) |
| **Long tail of clusters with `RNA.weight < 0.2`** | ADT/ATAC carries lineage signal RNA misses. These are the wins — fine-grained T-cell / B-cell / monocyte subsets that the protein panel resolves but RNA cannot |
| **Long tail of `RNA.weight > 0.8`** | RNA carries signal ADT misses. Typical for proliferating cells, stress states, ribosomal-program clusters — surface proteins are bottoming out on these |
| **Cluster with weight straddling 0.5 (wide IQR)** | Heterogeneous cluster — some cells are RNA-driven, others ADT-driven. Sub-cluster candidate; consider raising `resolution` |
| **All cells at exactly 0.5** | Bug — `modality.weight.name` collision (single-string fallback) OR `dims.list` mismatch (too many dims on one modality drowning the other) |
| **All cells at exactly 1.0** | One modality empty — verify `Reductions(obj)` includes the second-modality reduction AND `dims.list[[2]]` indexes existing dimensions |

## Expected `RNA.weight` ranges by tissue

These are empirical ranges from validated WNN runs in the published
Seurat / collaborator literature. Use as ballpark "is my number
sensible":

| Tissue | Typical median `RNA.weight` | Notes |
|---|---|---|
| Bone marrow (bmcite, healthy) | 0.40–0.50 | Slight ADT favor; lineage markers (CD34, CD38, CD3/4/8, CD19, CD14, CD16) carry strong signal for hematopoietic populations |
| Peripheral blood (PBMC, healthy) | 0.45–0.55 | Roughly balanced; T-cell subsets benefit from CD4/CD8/CCR7, monocyte subsets benefit from CD14/CD16 |
| Tonsil / lymph node | 0.30–0.45 | ADT-favored; B-cell developmental stages have sharp surface-marker progressions (IgD, CD27, CD38) |
| Tumor / TME | 0.50–0.65 | RNA-favored; activation/exhaustion programs (transcript) dominate over surface markers, especially for T cells |
| Healthy dissociated organ (liver, lung, kidney) | 0.55–0.70 | RNA-favored; surface markers less established outside hematopoietic compartments |
| RNA + ATAC multiome (PBMC) | 0.55–0.70 | RNA-favored; ATAC signal is sparser per cell, but contributes at lineage-defining loci |
| RNA + ATAC multiome (developing brain) | 0.45–0.60 | More balanced; ATAC captures cell-type identity through accessibility of regulatory regions |

If your median is **far outside** the expected range:

- Too high (close to 1.0) → second modality reductions are noise.
  Check `Reductions(obj)` and `dims.list` — the ADT or ATAC range
  may be empty or all noise. Check the per-modality elbow.
- Too low (close to 0.0) → RNA reduction is noise. Re-check
  `nFeature_RNA` filtering, batch effects, or the HVG count.

## The "is WNN worth it" decision

A user reasonably asks: "did joint clustering give us anything the RNA
clustering didn't?" Three signals to gather:

1. **Cluster count.** WNN clusters vs RNA-only clusters at the same
   resolution. Higher count (typical: +25% to +50%) is the headline
   number. If WNN gives the SAME count, the second modality didn't
   resolve anything new.

2. **UMAP topology comparison.** The Step-5 `umap_rna_only.png`
   shows the WNN cluster labels overlaid on the RNA-only UMAP. If
   the cluster polygons are intermingled on the RNA UMAP but separate
   on the WNN UMAP, those are the WNN wins. Quantify:

   ```r
   # Pairs of clusters whose ADT-driven separation is invisible on RNA UMAP:
   # ratio of intra-cluster distance vs inter-cluster distance in each space.
   library(dplyr)
   wnn_emb <- Embeddings(obj, "wnn.umap")
   rna_emb <- Embeddings(obj, "rna.umap")
   ids <- as.character(Idents(obj))
   centroids_wnn <- aggregate(wnn_emb, list(id = ids), median)
   centroids_rna <- aggregate(rna_emb, list(id = ids), median)
   # Clusters that are CLOSE on RNA but FAR on WNN are the wins.
   ```

3. **Per-cluster `RNA.weight` distribution.** Clusters with median
   `RNA.weight` < 0.4 are the ones the ADT/ATAC modality is driving;
   if those clusters wouldn't exist (or would be merged with others)
   under RNA-only clustering, WNN earned its keep.

The headline summary line for the user is something like:

> "WNN produced 23 clusters vs 15 on RNA-only (+53%). Median
> `RNA.weight` is 0.42 (slight ADT favor, typical for bone marrow).
> Clusters 13/16/17/21 are intermingled on the RNA-only UMAP and
> separate clearly on the WNN UMAP — those are the populations the
> protein panel resolved."

## When the weights argue for re-tuning

| Observation | Likely cause | Remedy |
|---|---|---|
| All cells at `RNA.weight ≈ 0.5`, no variation across clusters | `sd.scale` too high — bandwidth is washing out per-cell signal | `sd.scale = 0.5` in `FindMultiModalNeighbors` |
| RNA.weight bimodal (cells at 0.0 or 1.0, nothing in between) | `sd.scale` too low — kernel is over-sharpened | `sd.scale = 2` |
| Median weight near 1.0 with a narrow spread | Second modality has too few informative dims | Increase `dims.list[[2]]` upper bound (check elbow) OR re-run per-modality preprocessing |
| Median weight near 0.0 with a narrow spread | RNA reduction is noisy | Re-check `nFeature_RNA` filter, increase HVG count, re-run PCA |
| `RNA.weight` missing entirely | `modality.weight.name` single-string + non-RNA DefaultAssay (SCT.weight landed instead) | Pass 2-element vector: `c("RNA.weight", "ADT.weight")` |
| Weights look correct but joint UMAP looks identical to RNA UMAP | `dims.list[[2]]` underspecified — too few dims on the second modality | Bump second-modality dims |

## A note on `ADT.weight` (or `ATAC.weight`)

The partner column (`ADT.weight` for CITE-seq, `ATAC.weight` for
multiome) is mechanically `1 - RNA.weight` per cell when there are
exactly two modalities. It carries the same information as
`RNA.weight`; the recipe uses `RNA.weight` for plots by convention.
For three-modality WNN (RNA + ADT + ATAC), the three columns sum to
1 per cell and each carries independent information.

## Cross-checking with biological annotation

If you have cell-type labels (from a reference mapping run or curated
annotation), grouping `RNA.weight` by cell type — not cluster — is
sometimes more interpretable:

```r
# Assuming obj$celltype.l2 is the curated annotation.
VlnPlot(obj, features = "RNA.weight", group.by = "celltype.l2",
        sort = TRUE, pt.size = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red")
```

In bmcite, the Hao 2021 paper shows: progenitors (HSPC, GMP) cluster
at high `RNA.weight` (~0.7) — transcripts are the signal; T cells
(CD4 naive, CD4 memory, CD8 naive) cluster at low `RNA.weight`
(~0.3) — surface proteins are the signal. The recipe reproduces this
pattern on the validation fixture.

## References

- Hao Y., et al. (2021). *Integrated analysis of multimodal single-cell
  data.* Cell 184, 3573–3587. doi:10.1016/j.cell.2021.04.048 — Figure 1
  has the canonical interpretation of progenitor (RNA-driven) vs
  T-cell (ADT-driven) WNN weights in bmcite.
- Seurat WNN vignette — the "Visualize modality weights" section
  is the worked-example reference for the sorted VlnPlot pattern.
