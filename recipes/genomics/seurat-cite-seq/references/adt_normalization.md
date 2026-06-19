# ADT normalization — CLR semantics, margin choice, DSB alternative

Deep detail behind Step 3 of `seurat-cite-seq`. Loaded when the agent needs
to defend the normalization choice, switch to DSB on an isotype-equipped
panel, or explain why a `margin=1` legacy script disagrees with the current
recipe.

ADT counts behave nothing like RNA counts. Each protein has its own
bimodal distribution: a low-background peak (non-specific binding,
ambient antibody in the droplet) and a specific-binding peak (cells that
genuinely express the antigen). Log-normalization of total counts — the
RNA default — flattens that bimodality. The CITE-seq community
(Stoeckius et al. 2017, Nat. Methods) settled early on the
**centered log-ratio (CLR)** transform from compositional-data theory,
which preserves the per-protein bimodality and produces ~0-centered
values where the bandwidth tracks the signal-to-background ratio.

## CLR `margin` semantics — `margin=2` is the v5 default

`NormalizeData(obj, normalization.method = "CLR", margin = ?, assay = "ADT")`
takes one of two values. The math is the same — `log(x_i / geometric_mean(x))`
— but the geometric mean changes:

| `margin` | Geometric mean across | Interpretation | When |
|---|---|---|---|
| `2` (default) | **cells**, per feature | Stabilizes per-protein bimodality; each protein becomes a ~0-centered distribution whose spread is the signal-to-background ratio. | **Seurat v3.2+ default, recommended for almost every panel.** Used in the multimodal vignette. Makes inter-protein comparisons in a FeaturePlot honest. |
| `1` (legacy) | **features**, per cell | Per-cell normalization across the panel — older convention. | Only when you're replicating an older pipeline or a paper that explicitly used `margin=1`. Otherwise legacy. Reading scripts that pre-date v3.2 — they often default to `margin=1` silently. |

`margin = 2` is what the Seurat multimodal vignette explicitly calls:

```r
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2)
```

Pick `margin=2` unless you have a hard reason. If the user supplies a script
that used `margin=1`, the comparison plots will look different
(per-cell-normalized) but they're not wrong — they're a different question.
Don't silently switch.

### How `margin=1` and `margin=2` differ in practice

| Question | `margin=2` (per-feature) | `margin=1` (per-cell) |
|---|---|---|
| "Is CD3 high on this cell relative to other cells?" | YES, directly | indirectly — relative to other proteins on the same cell |
| "Is CD3 high relative to CD4 on this cell?" | indirectly | YES, directly |
| FeaturePlot of `adt_CD3` | gradient is "high-vs-low CD3 across cells" | gradient is "CD3 relative to this cell's panel mean" — confounded by panel size |
| Biaxial CD4 vs CD8 | classical flow-style quadrants | shifted; quadrants depend on per-cell sums |

The biaxial flow-style framing (Step 5 of the SKILL.md) implicitly assumes
`margin=2` — that's what makes the quadrants stable across cells.

## Panel-size effects on normalization

CLR variance scales with the geometric mean of the panel. Three regimes:

- **Tiny panel (<10 markers).** The geometric mean is dominated by 2–3
  highly-expressed lineage proteins; CLR can over-correct. Use
  `margin=2` AND skip ADT-only clustering — the panel doesn't carry
  enough independent signal. Treat ADT as overlay only.
- **Mid panel (10–30 markers).** Sweet spot for CLR `margin=2`.
  Per-protein bimodality stabilizes; downstream ADT-PCA is sensible
  with `npcs = min(50, n - 1)`.
- **Large panel (>30 markers).** CLR is fine but DSB starts to win when
  isotypes + empty droplets are available — it explicitly models the
  background instead of just centering it.

If you do not know the panel size yet, the recipe's default — CLR
`margin=2`, `ScaleData(assay="ADT")`, optional `RunPCA` capped at
`npcs = n_adt - 1` with `approx=FALSE` — works for all three regimes;
the only difference is downstream clustering aggressiveness.

## DSB — the alternative when isotypes are in the panel

**Mulè, Martins, Tsang (2022, Nat. Commun.)** describe DSB
("Denoised and Scaled by Background"), the protein-aware alternative
to CLR. Two-step denoise:

1. Subtract the empty-droplet background (from the **raw**, unfiltered
   matrix — the droplets the cell-caller rejected).
2. Use isotype controls (e.g. Mouse-IgG1, IgG2a) to estimate residual
   technical variation per cell, then regress it out.

DSB ships in the `dsb` CRAN package; it is NOT in base Seurat.

```r
# DSB drop-in alternative to NormalizeData(..., method="CLR", margin=2).
# Requires the RAW (unfiltered) 10x matrix + isotype-control names.
if (!requireNamespace("dsb", quietly = TRUE)) install.packages("dsb")
library(dsb)

# raw_adt   <- full unfiltered ADT matrix (every droplet, not just cells)
# cells_adt <- the cell-passing ADT matrix (== adt_counts[, colnames(obj)])
# isotype_controls <- e.g. c("Mouse-IgG1", "Mouse-IgG2a", "Rat-IgG2b")
dsb_norm <- DSBNormalizeProtein(
  cell_protein_matrix         = cells_adt,
  empty_drop_matrix           = raw_adt[, !colnames(raw_adt) %in% colnames(obj)],
  denoise.counts              = TRUE,
  use.isotype.control         = TRUE,
  isotype.control.name.vec    = isotype_controls
)

# Write into the ADT data layer; downstream steps unchanged.
obj <- SetAssayData(obj, assay = "ADT", layer = "data", new.data = dsb_norm)
```

DSB requirements (all three):

1. **Raw unfiltered matrix on disk.** `raw_feature_bc_matrix` from
   CellRanger — the one with all droplets, not just cells.
2. **Isotype controls in the panel.** Antibodies against a target the
   cells don't express (e.g. mouse IgG isotypes on a human sample).
   Without isotypes DSB can run (`use.isotype.control = FALSE`) but
   loses its main denoising lever.
3. **A clearly bimodal QC plot showing empty vs cell droplets.** If
   the empty droplets carry as much signal as the cells, the panel has
   a contamination problem CLR can't fix either.

If you only have #1 (no isotypes), CLR `margin=2` is still the right
default — DSB without isotypes is worse than CLR with `margin=2`.

### When to NOT use DSB

- Panel has no isotype controls (most published panels under 25
  markers).
- Only filtered matrix on disk (no `raw_feature_bc_matrix`).
- Cross-sample comparison where empty-droplet backgrounds differ across
  runs — DSB's per-run subtraction can introduce batch effects. CLR
  `margin=2` is more comparable across samples for that case.

## Quick decision tree

```
panel size?
├── < 10 markers ─→ CLR margin=2; overlay only; no ADT clustering
├── 10–30 markers ─→ CLR margin=2 (default)
└── > 30 markers
    └── isotype controls + raw matrix available?
        ├── yes ─→ DSB
        └── no  ─→ CLR margin=2
```

## Sanity-check after normalization

Regardless of normalization choice, the per-protein `summary()` should
land in a sensible band:

```r
adt_norm <- GetAssayData(obj, assay = "ADT", layer = "data")
# CLR margin=2 — expect median near 0, IQR ~0.5-2, range ~ -2 to 6.
# DSB — expect median near 0, IQR ~0.5-1.5, range ~ -3 to 8.
# margin=1 (legacy) — expect median near 0 per CELL, panel-wide spread
#   depends on panel size — not directly comparable to margin=2 numbers.
print(t(apply(adt_norm[seq_len(min(6, nrow(adt_norm))), , drop = FALSE],
              1, summary)))
```

If the values are an order of magnitude off, the most common cause is a
forgotten `assay = "ADT"` argument — the call ran against the RNA assay
silently. Repeat with the explicit assay routing.

## References

- Stoeckius M., et al. (2017). *Simultaneous epitope and transcriptome
  measurement in single cells.* Nat. Methods 14, 865–868.
  doi:10.1038/nmeth.4380. The CITE-seq paper itself; introduces CLR for
  ADT.
- Mulè M.P., Martins A.J., Tsang J.S. (2022). *Normalizing and denoising
  protein expression data from droplet-based single cell profiling.*
  Nat. Commun. 13, 2099. doi:10.1038/s41467-022-29356-8. DSB paper.
- Seurat multimodal vignette
  <https://satijalab.org/seurat/articles/multimodal_vignette> —
  authoritative source for `margin=2` as the v5 default.
