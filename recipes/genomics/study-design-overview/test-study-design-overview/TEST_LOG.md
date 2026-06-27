# Test log — study-design-overview

**Validation status:** PASSED_WITH_NOTES
**Recipe content hash at validation:** 6279c36e6b4ef52e5a2346dda03ee11c42074d9ce98e6a35e3a6c2119cd2c847 (SKILL.md)
**Validated on:** 2026-06-27
**Validator:** aba-skill-authoring validation_loop (Phase 2)

## Fixture
- Dataset: GSE192391 sample metadata (`sample_metadata.csv`), pre-fetched into the test dir.
- Version pin: deterministic CSV checked into `test-study-design-overview/`; not re-downloaded.
- Shape: 30 samples × 11 columns (GSM, title, patient id, tissue, cell type, disease state,
  status, age, Sex, ethnicity, time). Metadata-only — exercises Steps 1–3.
- **Step 4 (PCA-by-factor) NOT validated** — requires an expression matrix; this metadata-only
  fixture does not provide one. `pca_by_factor.png` is therefore (correctly) not produced.
  Step 4 needs a separate expression fixture to validate.

## Environment
- Language: Python 3.12.13 (`/home/pkharchenko/aba/aba/.venv`)
- Key packages: dython 0.7.12, patsy 1.0.2, pandas 2.3.3, numpy 2.4.6, scipy 1.18.0,
  scikit-learn 1.9.0, seaborn 0.13.2, matplotlib 3.11.0, statsmodels 0.14.6
- Kernel: `python3` (ipykernel); nbconvert 7.17.1

## Notebooks
- `test_interactive.ipynb` (rendered: `test_interactive.html`, 432 KB) — Steps 1–3 verbatim
  (fixture path substituted for `"/path/to/sample_metadata.csv"`), inline embeds after each
  figure save, chdir into `outputs/` in setup, Step 4 a markdown-only "not exercised" cell.
  Executed with `--ExecutePreprocessor.kernel_name=python3`, NO `--allow-errors`. 0 errored cells.
- `test_batch.ipynb` (rendered: `test_batch.html`, 306 KB) — `## Batch variant` path with
  `acc="GSE192391"`: Step 1 + Step 3 only, no figures, per-study `confounding_report.csv`,
  one summary line + format/non-degeneracy assertion. 0 errored cells.
- Wall time: each notebook executes in a few seconds.

## Outputs
- `outputs/association_heatmap.png` — 62,156 bytes (1067×984 px)
- `outputs/design_balance.png` — 27,288 bytes (720×540 px)
- `outputs/confounding_report.csv` — 723 bytes (15 pairs)
- `outputs/GSE192391/confounding_report.csv` — 723 bytes (batch per-study subfolder)
- (`pca_by_factor.png` correctly absent — Step 4 not exercised, no expression data)

## Examination

### Figures (n=2)
- `association_heatmap.png`: CHECK existence/non-empty (62 KB ≥ 10 KB) | CHECK white background
  (corner pixel #ffffff) | CHECK dimensions (9×8 in @ 120 dpi → 1067×984 px) | CHECK Reds
  sequential palette, colorbar labeled "association (V / η / |ρ|)" | CHECK clustered (dendrograms
  on both axes) | CHECK **strong disease/status/time block present** — status, disease state, time
  cluster bottom-right with cells 0.87 / 0.98 / 0.60 | CHECK annotations readable (2-dp), diagonal
  = 1.00 | CHECK title "Factor–factor association (clustered)" and all axis labels readable.
  No FAILs.
- `design_balance.png`: CHECK existence/non-empty (27 KB) | CHECK white background | CHECK
  dimensions (6×4.5 in @ 120 dpi → 720×540 px) | CHECK title "Design balance: disease state ×
  time  [RED]" — correct worst pair (`report.iloc[0]`) with the RED flag in the title | CHECK
  content: COVID-19 = 12 (day 0) / 12 (day 7) / 0 (nan); control = 0 / 0 / 6 — off-diagonal
  zeros are the visual signature of confounding. No FAILs.
  - NOTE (not a FAIL): the `time` axis includes a literal `"nan"` column for the 6 control
    samples (which carry no timepoint). This is the honest, expected consequence of Step 1's
    `design[c].astype(str)` typing (NaN → the string `"nan"`) flowing into the verbatim
    `pd.crosstab(design[a], design[b])`; it does not gate any verdict (the RED flag comes from
    the rank check on dropna'd data) and it actually makes the control-vs-COVID design split
    clearer. Recorded as a cosmetic observation, not a deficiency.

### Values
- Step 1 typing: GSM/title → drop:id-like; patient id → drop:high-cardinality; tissue/cell type
  → drop:constant; age → continuous; disease state/status/Sex/ethnicity/time → categorical.
  Kept "6 factors (5 categorical, 1 continuous) across 30 samples". → PASS (matches expected)
- Step 2 top pairs: disease state×time = 0.98, disease state×status = 0.87, status×time = 0.60,
  Sex×status = 0.40. → PASS (within ±rounding of expected ranges)
- Step 3 verdict: disease state×time = RED (rank 2 < ncol 3, aliased); status×time = RED
  (rank 3 < 4); disease state×status = AMBER (assoc 0.87, full rank); status×Sex = GREEN (0.40).
  Counts: 12 GREEN / 2 RED / 1 AMBER. → PASS (matches expected exactly)
- technical factors detected: (none) — every pair bio×bio, as expected for this study. → PASS
- Worst pair (`report.iloc[0]`) = disease state × time, RED. → PASS
- Batch summary line (verbatim printed):
  `GSE192391: 6 factors, 30 samples | RED=2 AMBER=1 | worst: disease state×time RED`
  → PASS (matches documented format; assertion on format + `n_red > 0` passed)
- `produces:` files present in outputs/: 3 / 4 — association_heatmap.png, design_balance.png,
  confounding_report.csv present; pca_by_factor.png path-conditional on Step 4 (expression data),
  correctly absent. → PASS (path-conditional per recipe branching)

## Patches applied
- Round 1: (none — recipe ran verbatim with zero errored cells and all figure/value checks passing.)

## Sweep candidates (deferred to serial sweep pass)
- (none — this is a brand-new standalone recipe with no siblings to sweep, and no defect surfaced
  that would warrant a catalogue-wide grep.)

## Unresolved
- Step 4 (PCA colored by factor) is **not validated** — it requires an expression-matrix fixture
  (samples × genes) that this metadata-only fixture does not provide. To close this gap, validate
  Step 4 against an expression fixture (e.g. a small bulk/pseudobulk log-CPM matrix aligned to
  this metadata) in a follow-up run. The recipe itself instructs skipping Step 4 when only the
  metadata table is available, so this is a fixture-coverage gap, not a recipe defect.

## pip freeze (key subset)
python 3.12.13; dython 0.7.12; patsy 1.0.2; pandas 2.3.3; numpy 2.4.6; scipy 1.18.0;
scikit-learn 1.9.0; seaborn 0.13.2; matplotlib 3.11.0; statsmodels 0.14.6; nbconvert 7.17.1
