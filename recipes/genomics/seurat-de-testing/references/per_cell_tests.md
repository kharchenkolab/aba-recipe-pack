# Per-cell DE tests — every `test.use` value in `FindMarkers`

Every `test.use` value Seurat's `FindMarkers` / `FindAllMarkers` accepts,
what each statistic actually is, when it's the right choice, when it's the
wrong choice, and where to read the upstream paper. The Seurat v5 vignette
(satijalab.org/seurat/articles/de_vignette) demonstrates four of these in
passing; this reference is the full catalogue with the trade-off semantics.

Important framing: every test below operates **per cell** — each cell is one
observation. For CONDITION effects across BIOLOGICAL samples (stim vs ctrl,
disease vs healthy), no per-cell test is appropriate — go to pseudobulk
(Path C in SKILL.md, deep detail in `pseudobulk_workflow.md`). Squair et al.
2021 (Nat Commun, DOI 10.1038/s41467-021-25960-2) is the canonical reference
for why; the per-cell tests below are for WITHIN-sample comparisons (cluster
A vs cluster B, two cell-state groups within one donor).

## The catalogue at a glance

| `test.use` | What it is | Returns | Per-cell covariates? | Speed | Bioconductor? |
|---|---|---|---|---|---|
| `wilcox` (default) | Wilcoxon rank-sum (Mann-Whitney U) | `p_val`, `avg_log2FC`, `pct.1/2` | NO | Fast (presto: very fast) | No |
| `wilcox_limma` | Wilcoxon implemented via `limma::wilcoxTwoSampleTest` | same | NO | Fast | `limma` |
| `t` | Welch's t-test on log-normalized expression | same | NO | Fastest | No |
| `bimod` | Likelihood-ratio for a two-component mixture (McDavid 2013) | same | NO | Medium | No |
| `roc` | Per-gene classifier AUC | `myAUC`, `avg_log2FC`, `pct.1/2` (NO p_val) | NO | Medium | No |
| `LR` | Logistic regression of group on expression | same | **YES** (`latent.vars`) | Medium | No |
| `negbinom` | Negative-binomial GLM on UMI counts | same | **YES** | Slow | No (uses MASS) |
| `poisson` | Poisson GLM on UMI counts | same | **YES** | Slow | No |
| `MAST` | Two-part hurdle model (Finak 2015) | same | **YES** | Slow | `MAST` |
| `DESeq2` | NB GLM with size factors — **PSEUDOBULK ONLY** | same | (via design — see pseudobulk_workflow.md) | n/a | `DESeq2` |

The default values for `FindMarkers` in Seurat 5.5.0 (verified via
`formals(Seurat:::FindMarkers.default)`):

```
test.use        = "wilcox"
logfc.threshold = 0.1
min.pct         = 0.01
slot            = "data"        # uses log-normalized expression by default
latent.vars     = NULL
min.cells.feature = 3
min.cells.group   = 3
```

`FindAllMarkers` shares those defaults. The PBMC3k tutorial and most of our
recipes pin `logfc.threshold = 0.25` + `min.pct = 0.25` for stability across
Seurat patch releases (5.x has loosened the defaults more than once).

## Wilcoxon (`wilcox`) — the right default

Non-parametric rank-sum test between the two cell groups, per gene. Doesn't
assume any distribution on expression; robust to the mixed zero-inflated /
log-normal shape of scRNA-seq counts. Single-cell papers default to Wilcoxon
because it's:

- **Calibrated** under the per-cell null. For WITHIN-sample comparisons
  (cluster vs cluster) the p-values are usable as-is; Soneson & Robinson
  2018 (Nat Methods) ranked Wilcoxon among the best-performing per-cell DE
  methods in their benchmark of 36 approaches.
- **Fast** with `presto` (the C++ rank-test backend, immunogenomics/presto).
  Seurat v5 auto-uses `presto` when installed — 10–100× faster, identical
  p-values. The SKILL.md `## Install` block pulls it.
- **Symmetric** about zero log-fold-change — no winner/loser asymmetry the
  way a one-sided LR test has.

### When Wilcoxon WINS

- "What marks each cluster?" (Path A). The de facto standard.
- Two cell-states within ONE sample, no covariates to adjust for.
- Quick scan before deciding whether a deeper test is warranted.

### When Wilcoxon FAILS

- **Per-cell covariates exist** (cell-cycle phase, percent.mt, batch effect).
  Wilcoxon doesn't accept `latent.vars` — passing it errors with
  `unused argument`. Switch to MAST or LR.
- **Condition effects across BIOLOGICAL samples** (stim vs ctrl, disease vs
  healthy). Wilcoxon p-values are anti-conservative for cross-sample
  comparisons — they treat each cell as an independent replicate when
  cells from the same donor are correlated. Squair 2021 documents
  FPR inflation; pseudobulk is the fix. Path C.
- **Tiny groups** (<10 cells per side). Wilcoxon has limited resolution;
  consider `t` for speed or just don't run the test.

## Wilcoxon-limma (`wilcox_limma`) — same statistic, different backend

Wilcoxon implemented through `limma::wilcoxTwoSampleTest`, with slightly
different ties handling. Identical decisions for the typical dataset.
Useful only if `presto` is unavailable and `limma` is already loaded for
something else. Rarely needed.

## Welch's t-test (`t`)

Two-sample t-test on log-normalized expression, with unequal variances.

### When `t` WINS

- **Tiny groups (5–20 cells per side)** where the Wilcoxon's discrete
  p-value resolution gets bumpy. The t-test is parametric and gives finer
  p-values.
- **Very fast** for ad-hoc exploration.

### When `t` FAILS

- Expression isn't actually normal — single-cell data is zero-inflated +
  heavy-tailed. The t-test's distributional assumption is violated. For a
  formal report use Wilcoxon (calibrated under the per-cell null) or MAST
  (models the zero-inflation explicitly).
- Per-cell covariates: doesn't accept `latent.vars`.

## Bimod (`bimod`) — the original Seurat method

Likelihood-ratio test for a two-component mixture model — splits expression
into a "drop-out" (zero-inflated) component and a continuous component, fits
each side per group, LRT for the joint shift. McDavid et al. 2013
(Bioinformatics). Pre-MAST predecessor of the zero-inflated approach;
historical but functional.

### When `bimod` WINS

- Reproducing original Seurat-v1/v2-era analyses verbatim.

### When `bimod` FAILS

- Almost never the right choice today. MAST is the modern descendant; if
  you want a zero-inflated model, use MAST.

## ROC (`roc`) — per-gene classifier AUC

For each gene, treats the expression value as a 1-D classifier between
ident.1 and ident.2 and computes the area-under-ROC. Returns `myAUC`
(0.5 = no separation, 1.0 = perfect ident.1 marker, 0.0 = perfect ident.2
marker) and `avg_log2FC` — **NO `p_val`**.

### When `roc` WINS

- "What's the cleanest single marker for ident.1 vs ident.2?" When you
  want a ranking by **how cleanly** a gene separates the groups, not by
  significance under a null.
- Reporting markers for an antibody panel — AUC > 0.85 is a useful
  "this gene cleanly distinguishes the populations" threshold.

### When `roc` FAILS

- "Is gene X significantly different?" — there's no p-value. AUC = 0.55
  could be noise or a tiny real effect; the test doesn't tell you which.
- Use AUC as a RANKING criterion, not a significance criterion.

## Logistic regression (`LR`)

Per-gene logistic regression of group membership on expression, optionally
adjusted for `latent.vars`. Lighter-weight covariate-adjusted alternative
to MAST.

```r
de <- FindMarkers(obj, ident.1 = "A", ident.2 = "B",
                  test.use = "LR",
                  latent.vars = c("percent.mt", "S.Score", "G2M.Score"))
```

### When `LR` WINS

- Per-cell covariates needed but MAST is unavailable / too slow.
- The categorical covariates are already one-hot encoded as numeric in
  `obj@meta.data` (LR handles `latent.vars` as numeric design columns).

### When `LR` FAILS

- Expression isn't modeled as continuous — LR predicts group from
  expression as a univariate predictor, NOT a count model. For UMI-shaped
  count data with strong zero inflation, MAST is a better fit.

## Negative binomial (`negbinom`) and Poisson (`poisson`) — count GLMs

GLMs on the UMI count data (NOT the log-normalized layer). Requires the
counts layer to be available. Both accept `latent.vars`.

```r
de <- FindMarkers(obj, ident.1 = "A", ident.2 = "B",
                  test.use = "negbinom",
                  latent.vars = c("percent.mt"))
```

- `poisson` assumes mean = variance — usually violated (overdispersion).
- `negbinom` models overdispersion with a per-gene dispersion parameter.

### When `negbinom` / `poisson` WIN

- UMI counts with strong overdispersion AND a need for covariate
  adjustment, AND MAST is unsuitable (e.g. proper count-level inference
  preferred over log-data + hurdle).

### When they FAIL

- **Slow.** Per-gene GLM fitting over many cells. For >10k cells and
  >5k genes, MAST is faster and similarly principled.
- The same per-cell vs pseudobulk caveat: still per-cell, still wrong
  for cross-sample condition effects.

## MAST (`MAST`) — the modern zero-inflated default for per-cell with covariates

Finak et al. 2015 (Genome Biology, DOI 10.1186/s13059-015-0844-5). Two-part
hurdle model: a Bernoulli for the zero/nonzero (drop-out) component plus a
Gaussian on the continuous-positive component, both modeled per gene with
shared covariate design. Returns the joint LRT p-value for the group effect.

Seurat's `MAST` wrapper (verified via `Seurat:::MASTDETest`):

```r
fmla <- as.formula(paste0(" ~ ", paste(c("condition", colnames(latent.vars)), collapse = "+")))
zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
summaryCond <- MAST::summary(object = zlmCond, doLRT = "conditionGroup2")
```

`latent.vars` MUST be column names in `obj@meta.data`. The Seurat wrapper
treats them as-is in the `zlm` design — categorical variables need one-hot
encoding before being passed (or treat them via the `formula` to MAST
directly, out of Seurat's wrapper).

### When MAST WINS

- Two-group within-sample comparison with per-cell covariates that need
  conditioning (cell-cycle phase, percent.mt, percent.ribo).
- The zero-inflation matters for your gene set — MAST jointly models the
  "is this gene detected?" question and the "what level when detected?"
  question, which Wilcoxon and t-test conflate.

### When MAST FAILS

- **Condition effects across samples.** Even with `latent.vars = "donor"`,
  MAST treats donor as a per-cell covariate, NOT a replicate unit. The
  effective sample size is still N_cells, so p-values stay
  anti-conservative for cross-sample inference. This is one of the most
  common live-session mistakes; see `anti_patterns.md`.
- **Slow.** On large objects MAST is the slowest per-cell option. Wilcoxon
  with presto is 50–500× faster.

## DESeq2 (`DESeq2`) — pseudobulk only

The `test.use = "DESeq2"` slot is **only correct on a pseudobulk object** —
one where columns are sample × cell-type aggregates, not individual cells.
The Seurat wrapper builds a `~group` design (`Seurat:::DESeq2DETest`):

```r
group.info[cells.1, "group"] <- "Group1"
group.info[cells.2, "group"] <- "Group2"
dds <- DESeq2::DESeqDataSetFromMatrix(countData = data.use, colData = group.info, design = ~group)
dds <- DESeq2::estimateSizeFactors(dds)
dds <- DESeq2::estimateDispersions(dds, fitType = "local")
dds <- DESeq2::nbinomWaldTest(dds)
```

This is fine when `data.use` is a pseudobulk matrix (rows × pseudobulk
columns). It is **catastrophically wrong** when `data.use` is a per-cell
matrix — the dispersion fit is meaningless across cells from the same
biological sample, and the false-positive rate inflates 10–100× (Squair
2021 Fig. 1).

### When DESeq2 WINS

- Path C — pseudobulk condition effect with a simple two-level design
  (ident.1 = celltype_STIM, ident.2 = celltype_CTRL).

### When DESeq2 FAILS

- Per-cell input (always wrong).
- Multi-factor designs — the Seurat wrapper hardcodes `~group`. For
  `~donor + condition`, paired designs, or interaction terms, drop to
  native DESeq2 / edgeR on the `GetAssayData(pseudo, layer="counts")`
  matrix. Details in `pseudobulk_workflow.md` §"Multi-factor designs".

## Decision tree — pick the right `test.use`

```
WITHIN-SAMPLE comparison?
├── YES, no per-cell covariates              → "wilcox" (default; presto-accelerated)
├── YES, per-cell covariates (MT, cell-cycle) → "MAST"   (or "LR" if MAST is too slow)
├── YES, want classifier ranking, not p-val   → "roc"    (myAUC)
└── YES, tiny groups (<20/side), exploratory → "t"      (faster than wilcox; report as exploratory)

CROSS-SAMPLE condition effect (stim vs ctrl across donors)?
└── ALWAYS use Path C: AggregateExpression → FindMarkers(test.use = "DESeq2")
    Per-cell tests with `latent.vars = "donor"` are NOT a substitute (see anti_patterns.md).
```

## Sources

- Seurat v5 DE vignette — satijalab.org/seurat/articles/de_vignette
- `Seurat:::FindMarkers.default` formals (Seurat 5.5.0) — for argument defaults
- `Seurat:::DESeq2DETest`, `Seurat:::MASTDETest` (Seurat 5.5.0) — for design choices
- Finak et al. 2015 — MAST, Genome Biology (DOI 10.1186/s13059-015-0844-5)
- McDavid et al. 2013 — bimod LR mixture, Bioinformatics
- Soneson & Robinson 2018 — Bias, robustness and scalability in single-cell DE analysis, Nat Methods (DOI 10.1038/nmeth.4612)
- Squair et al. 2021 — Confronting false discoveries in single-cell DE, Nat Commun (DOI 10.1038/s41467-021-25960-2)
- presto — github.com/immunogenomics/presto
