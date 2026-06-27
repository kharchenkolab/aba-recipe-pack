---
name: study-design-overview
description: Survey the experimental design of a study's sample metadata — surface the design structure, factor-to-factor associations (clustered Cramér's V heatmap), and CONFOUNDING / aliasing between biological and technical factors, before any downstream analysis. Outputs 2-3 scientist-facing figures plus a RED/AMBER/GREEN confounding verdict.
when_to_use: You have a sample × metadata table (e.g. a GEO / SRA / ArrayExpress study just pulled in) mixing biological factors (disease, genotype, treatment, sex, age) and technical factors (batch, run, platform, date), and you want to understand the design, spot batch/condition confounding, and learn whether downstream differential analysis is even estimable — BEFORE running it. Also for "is X confounded with batch?", "test this study design", "what's confounded here".
avoid_when: Not for running the downstream analysis itself (DE, clustering, integration) — this is the design-QC step upstream of those. Not for single-factor power/sample-size calculations.
invocation: interactive+batch
requires_tools: [run_python]
capabilities_needed: [dython, scipy, statsmodels, scikit-learn, pandas, numpy, matplotlib, seaborn]
keywords: [experimental design, study design, confounding, confounder, batch effect, aliasing, Cramer's V, CramersV, association heatmap, metadata, covariate, design matrix, collinearity, orthogonality, balance, factor interaction, GEO, SRA, sample metadata, identifiability, estimable, eta squared, correlation ratio]
produces: [association_heatmap.png, design_balance.png, confounding_report.csv, pca_by_factor.png]
domain: genomics
source: "Methods synthesis: Bergsma (2013) bias-corrected Cramér's V; Fox & Monette (1992) GVIF; design-matrix rank/aliasing per the limma/edgeR non-full-rank diagnostics; dython.nominal.associations. Function signatures source-verified against dython 0.7.12. Validated on GSE192391."
---

# Study design overview — structure, association, and confounding

Take a study's **sample × metadata table** (one row per sample, columns are
candidate factors) and produce a scientist-facing overview of its experimental
design: which factors exist and how balanced they are, how strongly factors
co-vary (a clustered association heatmap), and — the part that matters most —
**which biological and technical factors are confounded**, with a verdict on
whether downstream differential analysis is even estimable.

The non-negotiable idea this recipe encodes: **association is graded and often
survivable; aliasing (perfect confounding) is binary and fatal.** Two factors
can have Cramér's V = 1.0 (perfectly aliased — their effects are mathematically
inseparable) or V = 0.6 (strong but still separable with enough samples). A
single correlation number cannot tell those apart — so this recipe reports both
the association magnitude (Step 2) AND the design-matrix rank verdict (Step 3),
and never collapses them into one number.

Base overview = 2-3 figures: a **clustered association heatmap** ("what's
related"), a **design-balance cross-tab** for the worst pair ("what's fatal"),
and — when an expression matrix is available — **PCA colored by factor** ("what
actually structures the data"). For deeper diagnostics on follow-up, see the
references below.

## Bundled references — load on demand

1. `references/methods_and_interpretation.md` — for the statistics behind each
   measure (why bias-corrected Cramér's V, the association-vs-aliasing
   distinction, the rank verdict) and the exact interpretation thresholds to
   show the scientist. Load before reporting numbers.
2. `references/advanced_diagnostics.md` — for the **broader follow-up tests**
   when the user wants more than the base overview: PVCA / variancePartition
   (rank variance sources), gPCA (a batch-effect test), MCA (multi-way
   structure), residual-shaded mosaic plots, the missingness map, Theil's U
   nesting detection, and GVIF (partial-confounding severity). R and Python
   functions for each.

## Install

```python
# dython supplies the unified mixed-type association matrix (bias-corrected
# Cramér's V / correlation ratio / Spearman) and the clustered heatmap. patsy
# (ships with statsmodels) builds the design matrix for the rank/aliasing check.
# scipy, scikit-learn, pandas, numpy, matplotlib, seaborn are in the base stack.
try:
    import dython, patsy  # noqa: F401
except Exception:
    pass  # call ensure_capability(name="dython") once; patsy comes via statsmodels
```

If `import dython` fails, run `ensure_capability(name="dython")` (one call,
idempotent) before Step 2.

## Decisions to surface up front

Tell the user these before running — they define what the overview means:

1. **Which columns are factors vs identifiers.** Sample IDs, accession numbers,
   and free-text descriptions are not design factors. Step 1 flags `drop:id-like`
   (one value per sample) and `drop:constant` (one value total); confirm the
   auto-classification, especially numeric-coded categoricals (a `batch` column
   coded 1/2/3 must be treated as categorical, not continuous).
2. **Which factors are biological vs technical.** The confounding verdict
   (Step 3) is computed for every **biological × technical** pair, because those
   are the ones that invalidate downstream DE. Step 3 uses a name heuristic
   (batch/run/lane/plate/platform/date/library/flowcell → technical); correct it
   if the study uses non-obvious names.
3. **Missing metadata.** Factors with high missingness are excluded from the
   measures (NaN is NOT treated as a real category — that would fabricate
   associations). If a key factor is sparse, say so.
4. **Bias correction is on.** Cramér's V is bias-corrected (Bergsma) by default —
   essential for the small, unbalanced cohorts typical of GEO; raw V over-reports
   near 1.0 at small n.

Figures shown as the analysis proceeds: `association_heatmap.png` (Step 2),
`design_balance.png` (Step 3), and `pca_by_factor.png` (Step 4, if expression
data is available).

---

## Step 1 — Load and type the sample metadata

Get the metadata into a DataFrame `meta` (rows = samples, columns = candidate
factors), then classify every column so the analysis uses sensible types and
drops non-informative ones.

```python
import pandas as pd, numpy as np

# Rows = samples, columns = candidate factors. If you already have this in memory
# (e.g. from a GEO fetch), use that DataFrame directly and skip the read.
meta = pd.read_csv("/path/to/sample_metadata.csv")

# Classify each column. GEO/SRA metadata arrives as STRINGS, so numeric factors
# (age, RIN) look like objects — test whether a column PARSES as numeric, don't
# trust the dtype; and an ID-NAMED column (patient id, donor) is a label, never a
# continuous measurement, even when coded with numbers. Heuristics (override per
# study where wrong):
#   n_unique <= 1                  -> drop:constant   (no design information)
#   n_unique == n_samples          -> drop:id-like    (a sample identifier)
#   pct_missing > 50%              -> drop:sparse      (too incomplete to use)
#   numeric (and not id-named) AND n_unique>10 -> continuous   (age, RIN, date)
#   categorical with n_unique > n/2 -> drop:high-cardinality   (near-id grouping,
#       e.g. patient/subject — every patient-level attribute is trivially nested
#       within it, which would flood the verdict with meaningless aliases; set it
#       aside by default, re-include only for an explicit repeated-measures model)
#   otherwise                      -> categorical
import re
ID_NAME = re.compile(r"\b(id|ids|patient|subject|donor|sample|accession|gsm|srr|barcode|uuid)\b", re.I)
def parses_numeric(s):
    return pd.to_numeric(s, errors="coerce").notna().mean() >= 0.8   # ≥80% are numbers
n = len(meta)
rows = []
for c in meta.columns:
    nun  = int(meta[c].nunique(dropna=True))
    miss = float(meta[c].isna().mean())
    is_num = (not ID_NAME.search(c)) and (pd.api.types.is_numeric_dtype(meta[c]) or parses_numeric(meta[c]))
    if   nun <= 1:              role = "drop:constant"
    elif nun >= n:             role = "drop:id-like"
    elif miss > 0.5:           role = "drop:sparse"
    elif is_num and nun > 10:  role = "continuous"
    elif nun > n / 2:          role = "drop:high-cardinality"
    else:                      role = "categorical"
    rows.append(dict(column=c, dtype=str(meta[c].dtype), n_unique=nun,
                     pct_missing=round(100*miss, 1), suggested=role))
summary = pd.DataFrame(rows)
print(summary.to_string(index=False))

# Keep informative factors; coerce continuous to numeric, categoricals to plain
# str. Use astype(str), NOT astype("string") — dython's NaN handling rejects the
# pandas nullable 'string' dtype (it silently returns 0 associations).
factors  = summary.loc[summary.suggested.isin(["categorical", "continuous"]), "column"].tolist()
num_cols = summary.loc[summary.suggested == "continuous", "column"].tolist()
cat_cols = [c for c in factors if c not in num_cols]
design = meta[factors].copy()
design[num_cols] = design[num_cols].apply(pd.to_numeric, errors="coerce")
for c in cat_cols:
    design[c] = design[c].astype(str)
print(f"\nUsing {len(factors)} factors ({len(cat_cols)} categorical, "
      f"{len(num_cols)} continuous) across {n} samples")
```

**Report:** print the `summary` table (column, dtype, n_unique, pct_missing,
suggested) and the kept-factor line. Tell the user which columns are dropped and
why, and ask them to correct any misclassification before continuing.

---

## Step 2 — Clustered mixed-type association heatmap

One symmetric factor×factor matrix where each cell uses the right measure for
that pair's types, hierarchically clustered so confounded blocks sit together.

```python
from dython.nominal import associations
import seaborn as sns, matplotlib.pyplot as plt

# cat×cat = bias-corrected Cramér's V; cat×continuous = correlation ratio η;
# continuous×continuous = Spearman ρ. Compute the matrix with dython, then render
# with seaborn.clustermap so correlated factor groups form visible blocks. (Use
# clustermap, NOT dython's clustering=True — that path calls sch.distance.pdist
# and is broken against recent scipy.) dython RETURNS A DICT — the matrix is
# res["corr"] (NOT res.corr).
res = associations(
    design,
    nominal_columns=cat_cols,            # be explicit; 'auto' can mis-detect numeric-coded categories
    nom_nom_assoc="cramer",
    nom_num_assoc="correlation_ratio",
    num_num_assoc="spearman",            # robust default for metadata (age, RIN, date)
    cramers_v_bias_correction=True,      # essential for small / unbalanced GEO cohorts
    nan_strategy="drop_samples",         # don't let NaN become a phantom category
    compute_only=True,                   # we render the heatmap ourselves (next)
)
assoc = res["corr"].astype(float)        # the N×N association matrix (DataFrame)
mag = assoc.abs().fillna(0)              # display magnitude (|ρ| for continuous pairs), in [0,1]

g = sns.clustermap(mag, cmap="Reds", vmin=0, vmax=1, annot=True, fmt=".2f",
                   figsize=(9, 8), dendrogram_ratio=0.12,
                   cbar_kws={"label": "association (V / η / |ρ|)"})
g.figure.suptitle("Factor–factor association (clustered)", y=1.02)
g.savefig("association_heatmap.png", dpi=120); plt.close("all")

# Rank the off-diagonal pairs — the ones to scrutinize.
m = mag.copy(); np.fill_diagonal(m.values, np.nan)
pairs = (m.stack().reset_index()
           .rename(columns={"level_0": "factor_a", "level_1": "factor_b", 0: "assoc"}))
pairs = pairs[pairs.factor_a < pairs.factor_b].sort_values("assoc", ascending=False)
print(pairs.head(10).to_string(index=False))
```

**Report:** save `association_heatmap.png` and print the top-10 pair table. Tell
the user: an association ≥0.5 between a **biological and a technical** factor is
the kind to flag; ≈1.0 is candidate aliasing that Step 3 will adjudicate. Read
these alongside per-cell sample sizes — at small n even bias-corrected V is
unstable. For the magnitude bands and what each measure means, read
`references/methods_and_interpretation.md`.

**Pitfall:** dython returns a **dict** — use `res["corr"]`, not `res.corr`.
Spearman cells (continuous×continuous) are signed (−1..1); the heatmap shows
magnitude.

---

## Step 3 — Confounding / aliasing verdict (the action-driving step)

For each factor pair worth scrutinizing, decide whether the two are merely
associated (survivable) or **aliased** (their effects cannot be separated). The
backbone is design-matrix rank: build `~ A + B` and check whether its rank
equals its column count — rank-deficient ⇔ one factor is a linear function of
the other ⇔ the exact "coefficients not estimable" condition limma/edgeR/DESeq2
throw, surfaced here BEFORE the user runs DE. **biological × technical** pairs
are the most dangerous (a nuisance confounded with the effect of interest), but
biological×biological aliasing (e.g. disease↔status) is just as inseparable — so
the verdict is computed for every pair and annotated by type, then ordered by
danger.

```python
import numpy as np, pandas as pd, patsy, re
from itertools import combinations

# Name heuristic for technical factors — EDIT `technical` directly if the study
# uses non-obvious names; it is the fallback, not a guarantee. (Many GEO studies
# expose no technical column at all — then every pair is bio×bio, which is fine.)
TECH = re.compile(r"batch|run|lane|plate|platform|date|librar|flow.?cell|chip|machine|center|site|seq", re.I)
technical = [c for c in factors if TECH.search(c)]
def role(c): return "tech" if c in technical else "bio"
print("technical factors detected:", technical or "(none — edit `technical` if any apply)")

# Association magnitude from Step 2, keyed by the sorted pair.
assoc_of = {tuple(sorted([r.factor_a, r.factor_b])): float(r.assoc) for r in pairs.itertuples()}

def verdict(a, b):
    d = design[[a, b]].dropna()
    # rank < ncol  ⇔  one factor is a deterministic function of the other (aliased)
    X = patsy.dmatrix(f"Q('{a}') + Q('{b}')", d, return_type="dataframe")
    rank, ncol = int(np.linalg.matrix_rank(X.values)), X.shape[1]
    strength = assoc_of.get(tuple(sorted([a, b])), float("nan"))
    # empty design cells are only meaningful between two categoricals (a crosstab
    # over a continuous factor is degenerate) — reported as info, not the trigger.
    cat_pair = (a in cat_cols) and (b in cat_cols)
    empty = int((pd.crosstab(d[a], d[b]).values == 0).sum()) if cat_pair else None
    if   rank < ncol:        flag = "RED"      # aliased — effects inseparable
    elif strength >= 0.5:    flag = "AMBER"    # strong but separable — partial confounding
    else:                    flag = "GREEN"    # estimable, weak association
    return dict(factor_a=a, factor_b=b, pair_type="×".join(sorted([role(a), role(b)])),
                assoc=round(strength, 2), n=len(d), rank=rank, ncol=ncol,
                empty_cells=empty, flag=flag)

report = pd.DataFrame([verdict(a, b) for a, b in combinations(factors, 2)])
# Danger order: RED→AMBER→GREEN; biological×technical ahead of same-type; then assoc.
report["_k"] = report.flag.map({"RED": 0, "AMBER": 1, "GREEN": 2}) * 10 \
             + (report.pair_type != "bio×tech").astype(int)
report = (report.sort_values(["_k", "assoc"], ascending=[True, False])
                .drop(columns="_k").reset_index(drop=True))
report.to_csv("confounding_report.csv", index=False)

# Show the pairs worth attention: anything not GREEN, or at least moderately associated.
print(report[(report.flag != "GREEN") | (report.assoc >= 0.3)].to_string(index=False))
```

**Report:** save `confounding_report.csv`. State the headline verdicts plainly:
- **RED:** "`<factor_a>` is perfectly confounded with `<factor_b>` — their effects
  are mathematically inseparable; no batch-correction method can rescue downstream
  DE. Only re-design or dropping samples helps." (Most severe when the pair is
  biological×technical.)
- **AMBER:** "`<factor_a>` and `<factor_b>` are strongly but not perfectly
  associated (V/η ≥ 0.5) — estimable but with lost power; include the nuisance
  factor as a covariate and interpret with caution."
- **GREEN:** well-balanced; standard covariate adjustment suffices.

Then draw the balance cross-tab for the worst pair:

```python
import matplotlib.pyplot as plt, seaborn as sns

worst = report.iloc[0]                          # most dangerous pair after sorting
ct = pd.crosstab(design[worst.factor_a], design[worst.factor_b])
plt.figure(figsize=(6, 4.5))
sns.heatmap(ct, annot=True, fmt="d", cmap="Reds", cbar_kws={"label": "samples"})
plt.title(f"Design balance: {worst.factor_a} × {worst.factor_b}  [{worst.flag}]")
plt.tight_layout(); plt.savefig("design_balance.png", dpi=120, facecolor="white"); plt.close()
```

**Report:** save `design_balance.png` — samples per design cell for the flagged
pair. Off-diagonal zeros are the visual signature of confounding.

**Pitfall:** `Q('col name')` quotes columns with spaces/special characters in the
patsy formula. The `rank < ncol` test is the definitive verdict — do not infer
aliasing from the association magnitude alone (a strong-but-<1 V is still
separable).

---

## Step 4 — PCA colored by factor (only if expression data is available)

The only view that touches the actual high-dimensional data: it answers "which
factor actually structures the data," which the metadata-only measures cannot.
**Skip this step if you only have the metadata table** — say so and stop at
Step 3.

```python
import numpy as np, matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# X = samples × features (genes), rows aligned to design.index. For single-cell,
# aggregate to pseudobulk per sample first; for bulk, use log-CPM.
pcs = PCA(n_components=2).fit_transform(np.log1p(np.asarray(X)))
color_by = factors[:4]                          # the first few factors
fig, axes = plt.subplots(1, len(color_by), figsize=(4 * len(color_by), 4), squeeze=False)
for ax, c in zip(axes[0], color_by):
    for lev, sub in design.groupby(c, observed=True):
        idx = design.index.get_indexer(sub.index)
        ax.scatter(pcs[idx, 0], pcs[idx, 1], s=25, alpha=0.8, label=str(lev))
    ax.set_title(c); ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    ax.legend(fontsize=6, frameon=False)
fig.suptitle("PCA of expression, colored by factor")
fig.tight_layout(); fig.savefig("pca_by_factor.png", dpi=120, facecolor="white"); plt.close()
```

**Report:** save `pca_by_factor.png`. If a technical factor cleanly separates
PC1/PC2 while the biological factor is scattered, the dominant variation is
technical → correct or include batch before DE. Two factors that color the plot
identically are confounded (cross-check Step 3). To quantify rather than eyeball
this, promote to PVCA / variancePartition — see
`references/advanced_diagnostics.md`.

---

## Batch variant

If invoked as `Skill(skill="study-design-overview", args="batch")` — e.g. to
screen every study in a list — switch to a quiet, one-line-per-study path:

- Run Step 1 (typing) and Step 3 (verdict) only; **skip the figures** (Steps 2/4
  and the `design_balance.png` save) unless a RED verdict is found.
- Still write `confounding_report.csv` per study (into a per-study subfolder).
- Print exactly one summary line per study:

```python
worst = report.iloc[0]
n_red   = int((report.flag == "RED").sum())
n_amber = int((report.flag == "AMBER").sum())
print(f"{acc}: {len(factors)} factors, {n} samples | RED={n_red} AMBER={n_amber} | "
      f"worst: {worst.factor_a}×{worst.factor_b} {worst.flag}")
```

(`acc` = the study accession the orchestrator passed.) No per-step narration, no
advisor calls.

---

## Final response checklist

At completion, summarize for the user:

- **Input:** how many samples × how many factors used, and which columns were
  dropped (id-like / constant / too-sparse) and why.
- **Structure:** the strongest factor associations from the heatmap (top pairs +
  values), called out as biological×biological, technical×technical, or the
  dangerous biological×technical.
- **Confounding verdict:** every RED pair stated plainly (effects inseparable —
  redesign/drop only), every AMBER pair (covariate-adjust, reduced power), and
  confirmation the rest are GREEN. This is the headline result.
- **Variance (if Step 4 ran):** which factor dominates the top PCs.
- **Outputs produced:** `association_heatmap.png`, `confounding_report.csv`,
  `design_balance.png` (+ `pca_by_factor.png` if expression was available).
- **Caveats:** small-n instability of V; pairwise measures miss 3-way
  confounding (point to MCA / the full-model rank check in
  `references/advanced_diagnostics.md`); association ≠ causation.
