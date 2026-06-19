# Output interpretation — column semantics and biological sanity checks

What every column in a `FindMarkers` / `FindAllMarkers` result actually
means, what numerical ranges to expect, and the biological-coherence
sanity-checks that distinguish a successful DE run from a plausible-looking
failure.

Load this reference when:
- The user asks "what does `pct.1` / `pct.2` mean?" or similar column-semantics
  questions.
- A DE result has hits but you're not sure whether they're biologically
  plausible — you need a sanity-check checklist.
- You're writing the report and need to know what to highlight vs caveat.

## The result data.frame

Every per-cell DE test in Seurat (Wilcoxon, MAST, LR, t, bimod, negbinom,
poisson, DESeq2) returns a data.frame with these columns (`FindMarkers`)
plus a `cluster` column (`FindAllMarkers`):

| Column | Meaning | Range / units |
|---|---|---|
| `p_val` | Raw p-value from the test (gene-by-gene null hypothesis: no difference between groups) | `[0, 1]`, often very small |
| `avg_log2FC` | Mean log2 fold-change of ident.1 vs ident.2 on log-normalized expression | typically `[-3, 3]`; values >5 in either direction usually mean very-low-expression genes |
| `pct.1` | Fraction of cells in ident.1 expressing the gene (any UMI > 0) | `[0, 1]` |
| `pct.2` | Fraction of cells in ident.2 expressing the gene | `[0, 1]` |
| `p_val_adj` | Bonferroni-adjusted p-value, multiplied by `nrow(obj)` (total genes), clipped at 1 | `[0, 1]` |
| `cluster` (FindAllMarkers only) | The ident.1 of the (cluster vs rest) comparison | factor over cluster IDs |
| `gene` (added by recipe via `de$gene <- rownames(de)`) | The gene symbol | character |

The ROC test (`test.use = "roc"`) replaces `p_val` / `p_val_adj` with:

| Column | Meaning | Range |
|---|---|---|
| `myAUC` | Area-under-ROC of expression as a 1-D classifier (ident.1 = positive) | `[0, 1]`, 0.5 = no separation, >0.85 = clean marker |
| `power` | `abs(myAUC - 0.5) * 2` — re-centered AUC magnitude | `[0, 1]` |

## What each column tells you

### `p_val` and `p_val_adj`

The raw `p_val` is the test's p-value under its model assumptions (Wilcoxon
under no-difference null, MAST under the joint discrete + continuous null,
etc.). The Bonferroni-adjusted `p_val_adj = pmin(p_val * nrow(obj), 1)` is
what Seurat uses for significance reporting. `p_val_adj < 0.05` is the
default cutoff.

> **Why Bonferroni and not BH?** Seurat hardcoded Bonferroni in
> `FindMarkers` / `FindAllMarkers`. It is more conservative than
> Benjamini-Hochberg (`p.adjust(p, method = "BH")` for false-discovery
> rate). If the user is comparing to a BH-adjusted result from another
> tool, the counts will differ — Seurat's call is the more conservative
> one. Recompute BH from `p_val` if needed:
> `de$p_val_adj_BH <- p.adjust(de$p_val, method = "BH")`.

`p_val == 0` can appear when the test statistic is so extreme that the
backend reports zero. Seurat clips these for display via
`-log10(pmax(p_val_adj, .Machine$double.xmin))` — that's what the volcano
plot in Path B / C does.

### `avg_log2FC`

The log2 ratio of mean expression between ident.1 and ident.2. Sign
convention: **positive = UP in ident.1**, negative = UP in ident.2. The
"avg" here is per-gene across the group; the base-2 transform is on the
log-normalized expression layer (NOT the counts).

For the Seurat-default LogNormalize (CP10k + log1p), the relationship to
linear-scale fold-change is approximate; for the typical scRNA expression
range, `avg_log2FC = 1` corresponds to roughly 2× more expression in
ident.1, `avg_log2FC = 2` to ~4× more, etc. Don't read tight quantitative
fold-changes off `avg_log2FC` — it's a rank-and-direction quantity.

**Typical ranges:**
- |avg_log2FC| < 0.25 — small effect; Seurat's default `logfc.threshold`
  filters these out unless explicitly relaxed.
- 0.25 ≤ |avg_log2FC| < 0.5 — modest effect.
- |avg_log2FC| ≥ 0.5 — biologically meaningful (the recipe's "stricter"
  cutoff).
- |avg_log2FC| ≥ 2 — strong effect; often lineage markers.
- |avg_log2FC| > 5 — usually low-expression gene with one side near zero;
  inspect `pct.1` / `pct.2` to check.

### `pct.1` and `pct.2` — the expression-coverage check

The fraction of cells in each group with any UMI for the gene. These are
**not significance** but they are the BIOLOGICAL-PLAUSIBILITY anchor:

- **Lineage markers** typically show high `pct.1` (expressed in most cells
  of the cluster) AND low `pct.2` (rare in other clusters). Example: CD79A
  in B cells, `pct.1 ≈ 0.9`, `pct.2 ≈ 0.05`.
- **Transcriptional-state markers** (cell-cycle, stress, IFN response)
  show moderate-to-high `pct.1` in BOTH groups but a fold-change shift:
  ISG15 in stim vs ctrl is expressed in most cells either way; the
  difference is the magnitude.
- **`pct.1 < 0.1` AND `pct.2 < 0.1`** — the gene is rare in BOTH groups.
  `avg_log2FC` may be huge but the result is fragile (driven by a few
  cells). Treat as a hypothesis at best.

A useful diagnostic: `pct.1 - pct.2`. The Seurat `min.diff.pct` argument
filters on this (default `-Inf`, so no filtering). For clean markers,
`pct.1 - pct.2 > 0.3` is a stronger criterion than `avg_log2FC` alone.

### `cluster` (FindAllMarkers)

The ident.1 of each one-vs-rest comparison. Loop over `unique(markers$cluster)`
to inspect per-cluster behavior. Empty clusters (no significant markers)
are usually **not biologically distinct** from a neighboring cluster —
consider merging.

## Biological-coherence sanity checks

A DE result that runs cleanly but reports the wrong biology is the failure
mode this section exists to catch. Run these BEFORE writing the report.

### Path A — cluster markers

For each cluster, the top markers (by `avg_log2FC`, with `p_val_adj < 0.05`)
should match known lineage genes for the labelled cell type. The canonical
PBMC table:

| Cluster label | Expected top markers (high `pct.1`, low `pct.2`) |
|---|---|
| B (naive / memory) | CD79A, CD79B, MS4A1, CD19, IGHM |
| CD14 Mono | CD14, LYZ, S100A8, S100A9, FCN1 |
| CD16 Mono | FCGR3A, MS4A7, VMO1 |
| CD4 T (naive) | CCR7, LEF1, SELL, IL7R, CD3D/E |
| CD4 T (memory) | IL7R, CCR6, S100A4 |
| CD8 T | CD8A, CD8B, CD3D, GZMA, GZMK |
| NK | GNLY, NKG7, KLRD1, FCGR3A, GZMB, PRF1 |
| DC (myeloid) | FCER1A, CST3, CLEC10A |
| pDC | LILRA4, IL3RA, CLEC4C |
| Platelet / megakaryocyte | PPBP, PF4, ITGA2B, GP9 |
| Erythroid | HBB, HBA1, HBA2 |

If a cluster labelled "CD14 Mono" has top markers MS4A1 / CD19 / CD79A,
the labels are wrong (label mismatch) OR the clustering merged distinct
populations.

A cluster with **no markers at `p_val_adj < 0.05`** is either too small
(<10 cells) or not transcriptionally distinct — merge with a neighbor.

### Path B — two-group comparison

Examples that should hold up:

- **CD8 T vs NK** (cytotoxic-vs-cytotoxic): UP in CD8 — CD3D, CD3E, CD3G,
  CD8A, CD8B (T-cell receptor + co-receptor); UP in NK — GNLY, NKG7,
  KLRC1, KLRD1, TYROBP, FCER1G (innate-cytotoxicity components).
- **CD14 Mono vs CD16 Mono**: UP in CD14 — CD14, LYZ, S100A8, S100A9; UP
  in CD16 — FCGR3A (CD16 itself), MS4A7, VMO1.
- **B vs T**: UP in B — CD79A, CD79B, MS4A1; UP in T — CD3D, CD3E, CD3G.

If the comparison's top hits are batch-effect genes (mitochondrial,
ribosomal, MALAT1, NEAT1) or housekeeping genes (ACTB, GAPDH), the result
is **suspect**: most often a confound (cell-cycle, percent.mt, batch) is
driving the signal. Switch to MAST with `latent.vars = c("percent.mt",
"S.Score", "G2M.Score")` and re-test.

### Path C — pseudobulk condition effect (stim vs ctrl on ifnb)

The IFN-β stimulation in `ifnb` should put **interferon-stimulated genes
(ISGs) UP in STIM** vs CTRL. The canonical ISG set:

| Gene family | Members (top hits expected) |
|---|---|
| ISG core | ISG15, ISG20, IFI6, IFI27, IFI44L |
| IFIT family | IFIT1, IFIT2, IFIT3 |
| MX family | MX1, MX2 |
| OAS family | OAS1, OAS2, OAS3 |
| CXCL chemokines | CXCL10 (very strong), CXCL11 |
| Apolipoproteins | APOBEC3A, APOBEC3B |
| Other ISGs | RSAD2, IFITM3, TNFSF10 (TRAIL), GBP1, GBP5 |

A Path C run on ifnb that reports CXCL10, ISG15, IFIT1/2/3, RSAD2 UP in
STIM is methodologically sound. A run that reports housekeeping genes or
batch markers (orig.ident proxies) UP is broken — investigate.

Generic principle: every condition-effect Path C should have a
**predictable biology** that you can sanity-check the top hits against. If
you can't predict the biology before running the test, you can't sanity-
check the result; document that limitation in the report.

## Reporting the result — what to highlight

For Path A:
- Number of significant markers per cluster (range, mean).
- The top 3 markers per cluster (by `avg_log2FC`, filtered `p_val_adj < 0.05`).
- Clusters with <5 significant markers (weak — consider merging).
- Comparison to canonical lineage markers above (sanity-check).

For Path B:
- Total tested genes, # significant up / down.
- Top 10 by `avg_log2FC` (filtered `p_val_adj < 0.05`).
- Whether the comparison makes biological sense (the cytotoxic / B / T
  examples above).
- Covariates passed via `latent.vars` if MAST / LR was used.

For Path C:
- Design: # donors per condition, # cell-types tested, total pseudobulk
  columns, cell-types excluded for low replication.
- Per cell-type: # significant genes, top 10 up / down.
- Biological sanity: did the predicted up / down directions match? (For
  the IFN example: ISGs UP in STIM is the load-bearing check.)
- Limitations: any mocked donor_id (note explicitly — not a biological
  claim), batch effects not modeled, etc.

## Sources

- Seurat v5 DE vignette — satijalab.org/seurat/articles/de_vignette
- Seurat v5 `?FindMarkers` — column documentation
- PBMC canonical markers — Hao et al. 2021 (Seurat v4 paper, Cell, DOI 10.1016/j.cell.2021.04.048); Stuart et al. 2019 (Seurat v3 paper, Cell, DOI 10.1016/j.cell.2019.05.031)
- Interferon-stimulated gene set — Schoggins 2019 (Annu Rev Virol, DOI 10.1146/annurev-virology-092818-015756); the Interferome database (interferome.org)
- Squair et al. 2021 — context for why the biological-coherence check matters (Nat Commun, DOI 10.1038/s41467-021-25960-2)
