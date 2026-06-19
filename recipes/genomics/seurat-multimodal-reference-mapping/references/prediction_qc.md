# Prediction-score QC and ambiguity detection

How to interpret `predicted.<label>.score` values, when to flag
low-confidence cells, what bimodal score distributions mean biologically,
and the failure mode where a tissue-mismatched reference yields confident
wrong labels. Load this when the agent needs to defend a confidence
cutoff or explain why some cells look mislabeled despite high scores.

Sources consulted: Stuart et al. *Cell* 177 (2019) "Comprehensive
integration of single-cell data" (the `TransferData` machinery the
score comes from); Hao et al. *Cell* 184 (2021) (multimodal mapping);
Seurat source `R/integration.R` (`TransferData`'s prediction.score
construction).

## What the prediction.score is, precisely

After `TransferData` (or `MapQuery` calling it internally), each query
cell gets:

- `predicted.<key>` — the predicted label (the argmax of the prediction
  probability distribution).
- `predicted.<key>.score` — the maximum prediction probability across
  all labels in the reference's label set.

So `predicted.celltype.l1.score = 0.85` means "the most likely label
for this cell has weight 0.85" — i.e. the cell is confidently assigned
to that label and the next-best alternative has weight ≤ 0.15.

Mechanism (simplified): the score is computed from anchor-pair scores,
weighted by anchor-to-query-cell distances. A cell with many anchors
all pointing to one reference label gets a high score. A cell with
anchors split between two labels gets a score near 0.5.

The score is NOT a posterior probability under a calibrated model — it
is a heuristic, and its absolute value depends on the reference's label
set size and anchor density. So **cutoffs should be calibrated against
the actual score distribution of the query**, not a universal number.

## Calibrating a cutoff from the distribution

The recipe's default cutoff is 0.5. Refine by looking at the histogram
(produced as `predicted_score_hist.png` in Step 5):

```r
# Per-label distribution
df <- data.frame(
  celltype = query$predicted.celltype.l1,
  score    = query$predicted.celltype.l1.score
)

# Median + interquartile range per label
aggregate(score ~ celltype, df, function(x) c(med = median(x), q1 = quantile(x, 0.25), q3 = quantile(x, 0.75)))
```

Reading the histogram, look for:

| Pattern | What it means |
|---|---|
| Tight peak at >0.95 | Robust mapping — the reference covers this population well |
| Broad distribution centered at ~0.7 | Marginal — the reference is on-topic but not a perfect match |
| Bimodal (peaks at ~0.9 and ~0.4) | Ambiguous — query cells are split between two reference labels |
| Mass concentrated near 0.5 | Reference does not cover this population — labels are guesses |
| Long left tail extending below 0.3 | A subpopulation is poorly matched (often outside the reference's training distribution) |

A practical cutoff: pick the score at the local minimum between the
high-confidence and low-confidence modes, OR use 0.5 if there's no
clear bimodality. Flag low-confidence cells but DO NOT relabel them as
"unknown" by default — the recipe's downstream uses can still benefit
from the best-guess label, with confidence as a separate column.

## Bimodal scores within one predicted label

A common pattern: `predicted.celltype.l1 == "T cell"` but the score
histogram for the T-cell subset is bimodal. Interpretation: some query
cells are confident T cells (high score) and some are ambiguous between
"T cell" and a nearby label (CD8 vs CD4, T vs NK).

To resolve:

```r
# Look at the next-best label per cell. The full prediction probability
# vector is in query[["prediction.score.celltype.l1"]]@data (a matrix of
# label x cell). The second-highest label per cell is the candidate
# alternative.
P <- GetAssayData(query, assay = "prediction.score.celltype.l1", layer = "data")
# For each cell, find the top-2 labels and their weights:
top2 <- apply(P, 2, function(x) {
  o <- order(x, decreasing = TRUE)
  c(top1 = rownames(P)[o[1]], w1 = x[o[1]],
    top2 = rownames(P)[o[2]], w2 = x[o[2]])
})
top2 <- as.data.frame(t(top2))
# Now top2$top2 carries the runner-up; mass of bimodal cells should land
# on a consistent runner-up.
```

When the runner-up is consistent (e.g. ambiguous T cells mostly
runner-up to "NK"), that population is biologically intermediate (γδT
or NKT) or the reference's label resolution doesn't distinguish them.

## The confident-wrong-label failure mode

Recipe's Decisions #1 calls this out: a tissue-mismatched reference
will give confident-looking labels that are biologically wrong. The
mechanism: the reference's label space is closed (every query cell gets
assigned to SOME reference label), and the anchor scoring is local —
if there's any cell-type signal in common (e.g. mitochondrial-high
cells in two tissues), the reference will pick the closest label even
if no real match exists.

Symptoms:

- High median scores across the board (>0.7) but biology that doesn't
  make sense (heart tissue mapping to PBMC populations, all cells
  labeled "T cell").
- Some labels never used (e.g. the reference has 10 labels but the
  query is concentrated in 2 of them) — common in mismatch.
- The `predicted_ADT` imputed proteins don't match any sensible
  lineage marker pattern on the reference UMAP (see
  `references/multimodal_mapping_internals.md` for what predicted_ADT
  represents).

Defenses:

1. **Anchor-fraction check.** `nrow(anchors@anchors) / ncol(query) <
   0.1` is a red flag — see `multimodal_mapping_internals.md`'s anchor
   heuristic.
2. **Marker-gene overlay.** Plot canonical lineage markers (CD3, CD19,
   MS4A1, …) via FeaturePlot on `ref.umap` and confirm they track the
   predicted labels. Divergence = wrong reference.
3. **Truth column.** If the query has any independent label (manual
   annotation, expected sample composition), cross-tabulate against
   `predicted.celltype.l1`. Disagreement at >5% is a red flag.

## Cells outside the reference's training distribution

Some query cells genuinely don't have a counterpart in the reference —
rare populations, doublets, contaminating cell types. These map to
their nearest reference neighbor with low score. They show up as:

- Clustered low-score regions on `ref.umap` (use
  `FeaturePlot(query, "predicted.celltype.l1.score", reduction = "ref.umap")`
  to visualize the spatial distribution of scores).
- A right-tail mode in the histogram concentrated at one label (the
  whole population is mapping to its nearest neighbor with mid-range
  score).

These are biologically meaningful — they're the populations the
reference doesn't cover. Document them. Don't relabel them silently as
"unknown" without checking; the best-guess label sometimes still
carries useful broad lineage information.

## Score-threshold workflow (recommended pattern)

A defensible pipeline:

```r
# 1. Examine score distribution per label.
df <- data.frame(label = query$predicted.celltype.l1,
                 score = query$predicted.celltype.l1.score)

# 2. Optionally compute a per-label cutoff (median - 2 IQRs floor).
cutoffs <- aggregate(score ~ label, df, function(x) {
  med <- median(x); iqr <- IQR(x)
  max(0.3, med - 2 * iqr)   # never go below 0.3
})
print(cutoffs)

# 3. Add a confidence column.
query$predicted.celltype.l1.confident <- FALSE
for (lbl in unique(df$label)) {
  cut <- cutoffs$score[cutoffs$label == lbl]
  query$predicted.celltype.l1.confident[df$label == lbl] <- df$score[df$label == lbl] >= cut
}
table(query$predicted.celltype.l1.confident, query$predicted.celltype.l1)
```

This produces both the best-guess label (kept for everyone) and a
confidence flag (TRUE/FALSE) per cell — downstream filtering can use
the flag without losing the label information.

## Combining l1 + l2 + l3 confidence

When the recipe transfers multiple granularities, a cell can be
confident at l1 but ambiguous at l2 or l3 — that's normal (broad
identity is easier to call than fine identity). A useful summary
statistic: how often does the l1 label "agree" with the l2 → l1
mapping?

```r
# The reference's l2 -> l1 mapping (extracted from the reference itself):
l2_to_l1 <- setNames(reference$celltype.l1, reference$celltype.l2)
l2_to_l1 <- l2_to_l1[!duplicated(names(l2_to_l1))]

query$expected_l1_from_l2 <- l2_to_l1[query$predicted.celltype.l2]
agreement <- mean(query$predicted.celltype.l1 == query$expected_l1_from_l2,
                  na.rm = TRUE)
cat(sprintf("l1-from-l2 agreement: %.1f%%\n", 100 * agreement))
```

Below 95% agreement indicates the per-cell l1/l2 transfers are
internally inconsistent — usually a sign of cells near label boundaries
where the prediction is ambiguous at both levels.
