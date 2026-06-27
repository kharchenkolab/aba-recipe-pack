# Advanced diagnostics — the broader follow-up set

The base recipe (heatmap + rank verdict + PCA-by-factor) answers "what's related,
what's fatal, what's real" for the common case. Reach for these when the user
wants more depth, when the base verdict is ambiguous, or when the design is
complex (many factors, nesting, continuous covariates). Each entry: what it adds,
when to use it, and the exact function in R and/or Python.

## Variance structure — which factor actually drives the data

- **variancePartition (preferred) / PVCA.** The quantitative version of "which
  factors drive expression variance." Fits a per-gene (or per-PC) mixed model
  with all factors as random effects and reports the fraction of variance each
  explains — a bar chart of "% variance from disease / batch / sex / residual."
  **Flag when a technical factor explains a large share, especially more than the
  primary biological factor.** Handles continuous covariates; faster and more
  modern than classic PVCA.
  - R: `variancePartition::fitExtractVarPartModel(exprObj, formula, metadata)`
    then `plotVarPart(...)`. Formula uses `(1|batch)` for categoricals,
    `age` for continuous. Classic alternative: `pvca::pvcaBatchAssess(eset, batch.factors, threshold)`.
  - Python: no canonical equivalent — call variancePartition via `rpy2`, or
    approximate with per-PC `statsmodels` mixed models.

- **gPCA — a batch-effect *test* (not just a correlation).** Gives a δ statistic +
  permutation p-value for whether a specified batch variable carries real variance
  in the data. Answers "is this batch effect actually present?" rather than "is it
  correlated with condition."
  - R: `gPCA::gPCA.batchdetect(x, batch)` → δ + p.

## Partial-confounding severity (the AMBER quantifier)

- **Generalized VIF.** Quantifies *partial* confounding as variance inflation per
  model term, handling multi-df factors correctly (Fox–Monette, reported as
  GVIF^(1/(2·df)) for comparability). High GVIF = estimable but unstable.
  Complements the binary rank check with a graded severity (>√5 ≈ 2.24 concerning,
  >√10 ≈ 3.16 serious).
  - R: `car::vif(lm(y ~ A + B + C))` (auto-switches to GVIF for >1 df terms);
    `car::alias(lm(...))` lists the exact aliased coefficients;
    `caret::findLinearCombos(model.matrix(...))` enumerates linear combinations.
  - Python: `statsmodels.stats.outliers_influence.variance_inflation_factor(X, i)`
    over the columns of `patsy.dmatrix(...)`.

## Multi-way structure (what pairwise heatmaps miss)

- **MCA (Multiple Correspondence Analysis).** The categorical analogue of PCA, run
  *on the metadata itself* — projects all categorical factors into a 2D factor map
  where co-occurring levels sit close. Reveals 3+-way confounding the pairwise
  heatmap cannot.
  - R: `FactoMineR::MCA(df)` + `factoextra::fviz_mca_var(...)`.
  - Python: `prince.MCA().fit(df)` then `.plot(...)`.

- **Mosaic plot with residual shading.** The visual chi-square for one important
  factor pair: tile areas ∝ cell counts, tiles colored by standardized residuals
  (blue = over-represented vs independence, red = under). Shows *which* level
  combinations drive the association and exposes empty cells (a confounding red
  flag).
  - R: `vcd::mosaic(table, shade=TRUE, gp=vcd::shading_max, legend=TRUE)`;
    `vcd::assoc(...)` for tall tables.
  - Python: `statsmodels.graphics.mosaicplot.mosaic(...)` (no residual shading —
    R is preferred for this one).

## Intuitive balance & data-quality views

- **Design balance bar plots.** Stacked/grouped sample counts per factor level,
  optionally split by a second factor — the most intuitive "is this balanced?"
  view for non-statisticians; exposes singleton levels and imbalance.
  - R: `ggplot2::geom_bar(position = "dodge")`. Python: `seaborn.countplot(...)`.

- **Missingness map.** Heatmap of which metadata cells are NA, ordered by
  missingness — flags factors too sparse to use and non-random missingness (itself
  a confounder).
  - R: `naniar::vis_miss(df)` / `visdat::vis_dat(df)`.
  - Python: `missingno.matrix(df)` / `missingno.heatmap(df)`.

## Alternative association measures

- **Theil's U (asymmetric) — nesting detection.** U(A|B) ≠ U(B|A); uniquely
  detects hierarchy/nesting — e.g. U(batch|run)=1 but U(run|batch)<1 means run is
  nested within batch. Drives the random-effects structure for variancePartition.
  - R: `DescTools::UncertCoef(x, y, direction = "row"/"column")`.
  - Python: `dython.nominal.theils_u(x, y)`; `dython.nominal.conditional_entropy(...)`.

- **Normalized mutual information.** Symmetric cat×cat association robust to many
  levels; a non-chi-square alternative to Cramér's V.
  - R: `infotheo::mutinformation(...)` (with discretization).
  - Python: `sklearn.metrics.normalized_mutual_info_score(x, y)`.

## How to sequence a deeper investigation

1. Base recipe first (heatmap + rank verdict + PCA-by-factor).
2. If a verdict is AMBER and the user needs to model around it → **GVIF** for
   severity, **mosaic** to see which cells drive it.
3. If the user asks "which factor matters most for the data" → **variancePartition**
   (quantitative) and/or **gPCA** (a batch test).
4. If the design has many factors or you suspect 3-way confounding → **MCA**.
5. Always sanity-check data quality → **missingness map**, and **Theil's U** if
   you suspect nested technical factors.

## References

- Li, Boedigheimer et al. (2009). Principal Variance Components Analysis. (PVCA)
- Reese et al. (2013), *Bioinformatics*. gPCA — a guided-PCA batch statistic.
- Hoffman & Schadt (2016), *BMC Bioinformatics*. variancePartition.
- Friendly; Zeileis, Meyer & Hornik (2007). vcd strucplot + residual-based shadings.
- Fox & Monette (1992). Generalized collinearity diagnostics (GVIF).
