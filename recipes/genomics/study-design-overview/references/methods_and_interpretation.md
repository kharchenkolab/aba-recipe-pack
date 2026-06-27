# Methods & interpretation — what each measure means, and how to read it

## The one distinction that drives everything: three questions, three tools

The most common analyst error is collapsing these into a single "correlation":

1. **Association** — do two factors co-vary? → Cramér's V / correlation ratio η /
   Spearman ρ. A graded 0→1 quantity. *Some* association is normal and survivable.
2. **Aliasing / perfect confounding** — is one factor a deterministic function of
   another, so their effects are mathematically inseparable? → a yes/no property
   of the **design-matrix rank**, NOT a correlation threshold. This is the killer:
   it makes DE coefficients *not estimable*. Two factors can have V = 1.0 (perfect
   alias) or V = 0.6 (strong but separable) — the same heatmap cell magnitude
   range can mean either, which is why association alone is not a verdict.
3. **Variance structure** — which factors actually drive the high-dimensional
   data (not just each other)? → PVCA / variancePartition / PCA-vs-factor on the
   *expression matrix*, not the metadata table (see `advanced_diagnostics.md`).

Step 2 answers (1), Step 3 answers (2), Step 4 (and the advanced tier) answers
(3). Never report (1) as if it settled (2).

## The measures

- **Bias-corrected Cramér's V (cat × cat).** Cramér's V scales chi-square to
  [0,1]. Raw V is **upward-biased** at small n / many levels — routinely reports
  phantom ~1.0 associations on the small, unbalanced cohorts typical of GEO/SRA.
  The **Bergsma (2013)** correction subtracts a df-proportional term
  (φ̃² = max(0, φ² − (r−1)(c−1)/(N−1)), with adjusted r̃, c̃), shrinking V toward
  0 in proportion to degrees of freedom. Always keep `cramers_v_bias_correction=True`.
- **Correlation ratio η (cat × continuous).** η² = between-group SS / total SS —
  i.e. √(eta-squared from a one-way ANOVA of the continuous variable across the
  category's levels). In [0,1]; the categorical analogue of |correlation|.
- **Spearman ρ (continuous × continuous).** Rank correlation — robust to
  nonlinearity and outliers, the right default for metadata like age / RIN /
  date. Signed (−1..1); the heatmap shows magnitude.

## The aliasing verdict (Step 3) — why rank deficiency

Build the additive design `model.matrix(~ A + B)` (patsy `dmatrix("Q('A')+Q('B')")`).
If its **rank < its number of columns**, the columns are linearly dependent — one
factor is a linear function of the other, so a model containing both cannot
estimate their separate effects. This is exactly the *"Design matrix not of full
rank / the following coefficients are not estimable"* error limma, edgeR, and
DESeq2 raise — surfaced here before the user wastes a DE run. The cross-tab makes
it human-readable (block-diagonal counts = aliased); the rank check makes it
definitive. A strong-but-<1 association is rank-full → estimable, just lower-power.

## Interpretation bands to show the scientist

**Association magnitude** (same bands for V, η, and |ρ|):

| value | reading |
|---|---|
| < 0.10 | negligible |
| 0.10–0.30 | weak |
| 0.30–0.50 | moderate |
| ≥ 0.50 | strong — **flag if it's a biological × technical pair** |
| ≈ 1.0 | near-perfect — candidate aliasing → escalate to the Step 3 rank check |

Always state the per-cell **n** alongside the value: with N < ~100 and many
levels, even bias-corrected V is unstable. Show sample sizes, not just coefficients.

**Confounding verdict** (the action-driving flag):

- **RED — aliased:** design matrix rank-deficient (block-diagonal cross-tab).
  *Message:* effects inseparable; no batch-correction method rescues downstream
  DE; only re-design or excluding samples helps.
- **AMBER — partial:** estimable but empty design cells / strong-but-<1
  association (in the broader tier, GVIF^(1/2·df) > √5 ≈ 2.24 concerning,
  > √10 ≈ 3.16 serious). *Message:* include the technical factor as a covariate;
  expect reduced power/precision.
- **GREEN — balanced:** roughly orthogonal cross-tab; standard covariate
  adjustment suffices.

**η² (cat × continuous), e.g. condition vs RIN/age:** 0.01 small, 0.06 medium,
≥ 0.14 large. A large η² between a biological factor and a technical continuous
covariate (disease strongly predicts RIN) is a confounding warning.

## Housekeeping flags to always surface

- **Constant / redundant factors.** One-level factors carry no design info (drop).
  A factor that is a 1:1 relabeling of another (V = 1 with equal level counts, or
  Theil's U = 1 both directions) is redundant — keep one.
- **Nesting.** Asymmetric Theil's U (U(A|B)=1 but U(B|A)<1) means A is nested
  within B (e.g. run nested in batch) — matters for the random-effects structure
  in variancePartition. See `advanced_diagnostics.md`.
- **Missing metadata.** Exclude high-NA factors from the measures; never let
  NaN-as-category fabricate associations (the recipe uses `nan_strategy="drop_samples"`).
  Non-random missingness is itself a confounder — show the missingness map.
- **Pairwise ≠ joint.** The heatmap is pairwise; a clean pairwise picture can
  still hide 3-way confounding. MCA and the full-model rank check are the
  safeguards (`advanced_diagnostics.md`).

## References

- Bergsma, W. (2013). A bias-correction for Cramér's V and Tschuprow's T. *J. Korean Stat. Soc.*
- Fox, J. & Monette, G. (1992). Generalized collinearity diagnostics (GVIF). *JASA.*
- limma & edgeR user guides — non-full-rank design matrices (the "not estimable" condition).
- dython `nominal` module docs — Cramér's V w/ bias correction, correlation ratio, `associations`.
