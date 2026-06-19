---
name: scrna-doublet-detection
description: Decision guide for scRNA-seq doublet calling. Compares scDblFinder, DoubletFinder, scrublet, Solo, scds (cxds/bcds/hybrid), and vendor (10x) doublet flags. Use when the user asks "do I need doublet detection?", "which caller should I use?", or "where in the pipeline does doublet calling go?".
when_to_use: User is preparing a single-cell RNA-seq pipeline and wants to know whether to add a doublet caller, which one to use for their dataset/loading concentration/protocol, and where it sits relative to empty-droplet filtering and QC. Triggered by phrases like "doublet detection", "remove doublets", "scrublet vs scDblFinder", "DoubletFinder threshold".
avoid_when: User asks for a runnable step-by-step doublet recipe (use `bp-quality-control` directly). User's experiment is a designed multiplet recovery (e.g. cell-cell interaction by deliberate co-encapsulation) — then doublets are signal, not noise.
invocation: interactive
kind: knowhow_draft
requires_tools: [WebFetch, Read]
keywords: [doublet detection, scDblFinder, DoubletFinder, scrublet, Solo, scds, multiplet rate, 10x Chromium, heterotypic, homotypic, single-cell QC]
domain: genomics
source: "Xi & Li 2021 Cell Systems benchmark; Germain et al. 2021 F1000Research scDblFinder paper; scDblFinder Bioconductor vignette (Germain, current); Heumos et al. 2023 Single-cell Best Practices book; 10x Genomics Chromium Next GEM Single Cell 3' v3.1 user guide multiplet table."
audience: both
produces: []
capabilities_needed: []
---

# scRNA-seq doublet detection — which caller, and when do I even need one?

The question this knowhow answers: **"I ran a 10x scRNA-seq experiment; should
I call doublets, and if so, with which tool?"**

Audience: both naive (first scRNA-seq experiment) and experienced (picking
between callers for a specific protocol). The decision is multi-axis (loading
concentration × pipeline language × protocol × heterotypic-vs-homotypic
expectation). The executable workflow that runs doublet detection is
`bp-quality-control` (Python/scrublet by default; calls out scDblFinder for the
R path). Standalone per-caller recipes (scDblFinder, DoubletFinder, Solo) do
not yet exist in the catalogue — flagged in §7.

## 1. Quick decision

For the typical 10x scRNA-seq experiment (≥3000 cells recovered per sample,
standard 0.8%-per-1000-cells loading, mixed cell types):

| Your situation | Use |
|---|---|
| Default first pass, Python/scanpy pipeline | **scrublet** via `sc.pp.scrublet`, per sample (recipe: `bp-quality-control`) |
| You want the current best accuracy, R/Bioconductor pipeline | **scDblFinder** (Germain 2021 — beats scrublet/DoubletFinder on the same 16-dataset benchmark) |
| Already using `scvi-tools` / have a GPU latent ready | **Solo** (`scvi.external.SOLO`) — fits naturally if you already trained an scVI model |
| You want to ensemble two callers for high-stakes work | **scDblFinder + DoubletFinder** averaged scores (the scDblFinder vignette explicitly suggests this) |
| Very large dataset, only need a fast screen | **scrublet** or **scds/cxds** — fastest at the cost of some accuracy (Xi & Li 2021) |
| Cells recovered < ~1000 per sample, loading well below recommended | **Skip — read §3** |

If none of these fit (multi-species mixing experiment, deliberate doublet
recovery, plate-based protocol), see §3 (when NOT to do this) before §4.

The single most defensible default in 2026: **scDblFinder if you're in R,
scrublet if you're in Python** — both run per-sample, both produce a
`predicted_doublet` flag, both should run BEFORE strict QC (see §5).

## 2. Background — what's a doublet and what does the rate depend on?

A **doublet** is a single droplet/well that captured ≥2 cells; its
transcriptome is the sum of those cells' expression. Doublets distort
clustering (they form bogus "intermediate" clusters between real types) and
inflate downstream DE / trajectory artifacts.

The doublet rate is fundamentally a **loading-concentration** property, not a
biology property:

| Cells loaded (10x Chromium Next GEM 3' v3.1) | Cells recovered | Expected multiplet rate |
|---|---|---|
| ~800 | ~500 | ~0.4% |
| ~1,600 | ~1,000 | ~0.8% |
| ~3,200 | ~2,000 | ~1.6% |
| ~8,000 | ~5,000 | ~3.9% |
| ~16,000 | ~10,000 | ~7.6% |

(10x Genomics CG000315 Rev E user guide; ~0.8% per 1000 cells recovered is the
canonical rule-of-thumb.) HT (high-throughput) kits run at roughly half this
rate per recovered cell (scDblFinder vignette, Germain).

Crucially, **only heterotypic doublets** (two different cell types) are
detectable from expression alone — a doublet of two T cells looks like a
slightly fatter T cell. Every algorithm in §4 targets heterotypic doublets
(Germain 2021: "we do not attempt to identify homotypic doublets, which are
anyway virtually unidentifiable and relatively innocuous").

## 3. When NOT to call doublets

The "you don't need this" escape valve is high-value here — doublet calling
is over-applied and sometimes does harm.

- **Very small datasets / very low loading concentration.** At ~500 recovered
  cells you have ~0.4% × 500 ≈ 2 expected doublets. The caller's false-positive
  burden will dwarf the true positives; you lose more real cells than doublets.
  Just flag the high-`n_genes` tail in QC and move on.
- **The doublets ARE the biology.** Designed cell-cell interaction assays
  (e.g. deliberate T-cell : APC co-encapsulation, PIC-seq, RNA-seq of
  co-cultures designed to recover physical doublets) treat doublets as
  signal. Calling and removing them deletes the experiment.
  `# REVIEW(doublet-detection-methods): cite a representative PIC-seq /
  designed-doublet paper here so the reader has a concrete anchor (Giladi 2020?).`
- **Plate-based protocols (Smart-seq2/3, MARS-seq) at one cell per well.**
  Doublets are rare by construction (visual FACS QC at sort time). Calling
  computationally adds noise without removing meaningful signal.
  `# REVIEW(doublet-detection-methods): confirm — does the field actually skip
  doublet calling for Smart-seq2? Or is it just less common?`
- **You've already aggressively QC'd by `n_genes_by_counts`.** A 99th-percentile
  cap on `n_genes` (as `scrna-qc-clustering` does) removes the bulk of
  heterotypic doublets implicitly. Adding scrublet after that gives diminishing
  returns; the score histogram becomes unimodal and the caller's auto-threshold
  becomes unreliable.
- **Single-cell-type sorted populations** (e.g. FACS-purified naive CD4 T
  cells). Almost all doublets are homotypic and undetectable; the caller's
  output is mostly noise.
- **The biological question is "are there transitional/intermediate states?"**
  Some "doublet" calls in this case are real intermediate-state cells (e.g.
  doublet-positive thymocytes, hybrid EMT states). Run the caller for the flag
  but **don't auto-delete** — inspect cluster-by-cluster.

## 4. Alternatives — full matrix

| Method | Characterization | When it wins | When it fails | Cost | Recipe |
|---|---|---|---|---|---|
| **scDblFinder** (R/Bioconductor) | Simulates artificial doublets from real cells, trains gradient-boosted classifier on PCA features; cluster-based simulation captures rare types (Germain 2021) | Current best accuracy on the 16-dataset Germain benchmark; fastest of the accurate methods; sensible defaults; auto-thresholds | Homotypic doublets (by design); R-only — awkward to call from a scanpy pipeline | Seconds to a few minutes per sample; lightweight | TBD — no per-caller recipe; runs as the R branch of `bp-quality-control` |
| **DoubletFinder** (R/Seurat) | Simulates artificial doublets, projects with PCA, kNN density score (McGinnis 2019); requires a manual pK sweep | Top-ranked in Xi & Li 2021 benchmark (note: scDblFinder was not in that benchmark); ubiquitous in older Seurat papers | pK sweep is fiddly + slow; needs Seurat object; less convenient than scDblFinder for the same R workflow | Minutes-to-tens-of-minutes per sample (pK sweep dominates) | TBD — `# TODO(recipe): doubletfinder-r` |
| **scrublet** (Python) | Simulates artificial doublets from observed counts, kNN density on PCA-reduced space, k-modal threshold (Wolock 2019) | The scanpy-native default; one-line call `sc.pp.scrublet(adata)`; runs per sample; well-validated | Mid-pack accuracy in Xi & Li 2021 (below DoubletFinder/Solo); auto-threshold sometimes fails on unimodal score distributions and needs visual override | Seconds to a minute per sample | `bp-quality-control` (Python path) |
| **Solo** (`scvi.external.SOLO`) | Trains an scVI VAE, then a supervised classifier head on simulated doublets (Bernstein 2020) | If you already trained an scVI model for integration, doublet calling comes "for free"; tied for #2 in Xi & Li 2021 (AUROC) | GPU-class compute even for a small dataset; opaque failure mode; per-sample training is recommended → undoes the cost amortization | Minutes-to-tens-of-minutes on a GPU; hours on CPU | TBD — `# TODO(recipe): solo-doublet-scvi` |
| **scds** (cxds / bcds / hybrid) | cxds: co-expression score on gene pairs that should be mutually exclusive. bcds: binary classifier on simulated doublets. hybrid: average (Bais & Kostka 2020) | cxds is the fastest of all methods (Xi & Li 2021); good for a "screen 100 samples in minutes" budget | Mid-low accuracy on the Germain 2021 benchmark; less commonly used in 2024-2026 pipelines | Seconds | TBD — `# TODO(recipe): scds-doublet-r` |
| **10x Cell Ranger built-in / GEM-X multiplet flags** | Cell Ranger reports the expected multiplet rate from loading; some kits (Cell Multiplexing, CITE-seq with hashing) emit per-droplet multiplet calls directly from the assay | Hashing/multiplexing protocols: the vendor flag IS the ground truth — use it | Standard 3' / 5' v3.1 Cell Ranger does NOT call per-droplet doublets — it only reports the expected rate; you still need an algorithmic caller | Free; from existing Cell Ranger output | Read the Cell Ranger `summary.html`; combine with one of the above |
| **Hashing-based (HTODemux / hashedDrops)** | If samples were multiplexed with hashtags (HTO, MULTI-seq), the cross-hashtag droplets ARE confirmed doublets | The closest thing to ground truth in a standard workflow; should be used whenever hashing data is present | Requires hashed samples; doesn't catch intra-sample doublets (a doublet of two cells from the same hashtag still looks like a singlet) | Trivial | `# TODO(recipe): htodemux-r` (or covered by a future demultiplex recipe) |

**Key tension between Xi & Li 2021 and Germain 2021.** Xi & Li 2021
benchmarked nine methods and ranked DoubletFinder #1 — but **scDblFinder was
not in their benchmark.** Germain 2021 added scDblFinder to the same kind of
benchmark and reported it beats DoubletFinder and scrublet on heterotypic
detection at a fraction of the runtime. The Heumos et al. 2023 Single-cell
Best Practices book reconciles this by naming scDblFinder as the top
performer and recommending it as the default; scrublet is named as the
scanpy-native fallback. We adopt that consensus.

`# REVIEW(doublet-detection-methods): the Xi & Li vs Germain tension is the
single most important judgment in this doc. Confirm the current
(2026) consensus is still "scDblFinder > DoubletFinder", or surface a more
recent benchmark.`

## 5. Anti-patterns and common mistakes

- **Calling doublets on a batch-aggregated AnnData (multi-sample concatenated
  matrix).** Both scrublet and scDblFinder build a kNN graph in expression
  space; if samples are batch-confounded, the graph collapses into per-sample
  clouds and the caller flags whole samples as "doublets". **Always run
  per-sample.** scDblFinder's `samples=` argument is the explicit fix; for
  scrublet, loop over samples. (Heumos et al. 2023, "Per-sample doublet calls
  — aggregating batches breaks the detector"; scDblFinder vignette: "The most
  common cause for an unexpectedly large proportion of doublets is if you have
  a multi-sample dataset and did not split by samples.")

- **Calling doublets AFTER strict QC.** Counter-intuitive but well-documented:
  doublet rate is a function of cells **input to the caller**. If you've
  already deleted half the cells via MAD outlier QC, the caller's expected
  doublet rate is now wrong and the auto-threshold misbehaves. Worse, doublets
  formed between a good cell and a low-quality cell are no longer
  representable (the low-quality partner has been removed). **Run doublet
  detection BEFORE strict QC; let the doublet flag join the QC vector; then
  filter.** (scDblFinder vignette, explicit; Germain 2021: input "should not
  otherwise have undergone very stringent filtering".)

- **Calling doublets BEFORE empty-droplet filtering, though.** This is the
  asymmetric flip of the above: empty droplets must be removed first (via
  Cell Ranger's `filtered_feature_bc_matrix` or `emptyDrops`). The caller
  expects cell-containing barcodes. Order: empty-drop filter → doublet call
  → QC filter. (scDblFinder vignette implies; community consensus.)

- **Treating the caller's `predicted_doublet` boolean as the final word.**
  The Best Practices book explicitly recommends KEEPING flagged doublets
  initially, inspecting where they sit in the UMAP / which clusters they
  enrich, and then deciding. A cluster that is 80% flagged doublets is
  doublet-derived; a cluster with 5% flagged is fine.

- **Applying one caller's threshold blindly across protocols.** Drop-seq,
  10x v2, 10x v3, BD Rhapsody, and ParseBio all have different doublet rates
  per recovered cell. Don't carry over a `dbr=0.04` from a 10x run to a Parse
  experiment without recomputing the expected rate from the loading numbers.

- **Choosing DoubletFinder solely because Xi & Li 2021 ranked it #1.**
  That benchmark didn't include scDblFinder. The same author group's
  follow-up plus the broader 2023-2024 literature has shifted to scDblFinder
  as the default. Recency-by-paper-date isn't a good reason either — but
  the more recent benchmark with the broader method set is.

- **Skipping doublet calling because "we'll just filter high-`n_genes`
  cells".** That removes the obvious heterotypic doublets but misses doublets
  where one of the partners is small (e.g. RBC + T cell — the count depth
  is barely elevated). Use a dedicated caller and add the score to the QC
  vector.

- **Ensembling 5+ callers naively.** The scDblFinder vignette suggests
  averaging scDblFinder + DoubletFinder for modest improvement. Stacking
  every caller you can install is not a free win — false positives
  compound on cells where any method is sensitive to outliers.
  `# REVIEW(doublet-detection-methods): is there a published consensus
  ensemble strategy beyond Germain's "average two methods"?`

## 6. Sanity checks — how to know your choice was right

After running the caller, before deleting cells:

- **Doublet flag rate matches expectation.** If you loaded 8000 cells (expect
  ~3.9% multiplet rate per 10x table) and the caller flagged 18%, something
  is wrong — most likely batch-aggregation (you forgot per-sample). If the
  caller flagged 0.1%, the auto-threshold likely failed on a unimodal score
  distribution — inspect the histogram and threshold manually.

- **Doublet score histogram is bimodal.** scDblFinder + scrublet both produce
  a per-cell score. Plot the score distribution; you should see a real-cell
  mode near 0 and a doublet mode near the threshold. A flat or unimodal
  distribution means the caller didn't find a clear separation — either the
  doublet rate is genuinely tiny or the kNN landscape is too noisy.

- **Doublet-flagged cells localize in UMAP.** Plot the UMAP coloured by
  `predicted_doublet`. Genuine doublets cluster **between** singlet clusters
  (intermediate position in low-dim space) or form their own small
  high-`n_genes` clusters. A scatter of doublet calls uniformly distributed
  through every cluster suggests the caller is over-firing.

- **Cluster doublet fraction informs cluster identity.** Run clustering on
  the unfiltered object, then tally `predicted_doublet` fraction per
  cluster. Clusters with >40-50% doublet fraction are doublet-derived
  (delete the cluster); clusters with <10% are real (keep the cluster,
  delete the per-cell flagged cells if desired).

- **Top-doublet cells have elevated `n_genes_by_counts` AND a mix of
  marker genes from ≥2 cell types.** Spot-check the top-10 highest-score
  cells: if you see, say, both T-cell markers (CD3D) and myeloid markers
  (LYZ) in the same cell, that's a real heterotypic doublet. If the top
  hits look like normal cells of one type with slightly high counts, the
  caller is partly firing on count-depth outliers, not biology.

- **Don't trust the result if you ran one caller on a mixed-protocol
  pool.** Each protocol's doublet rate is different. Run per-protocol-batch
  or set the `dbr` parameter per-batch (scDblFinder supports this; for
  scrublet, loop).

## 7. See also

**Benchmarks + reviews:**
- Xi & Li 2021, Cell Systems — nine-method doublet caller benchmark,
  ranks DoubletFinder #1 (note: scDblFinder not included)
  (https://doi.org/10.1016/j.cels.2020.11.008)
- Germain et al. 2021, F1000Research — scDblFinder paper and 16-dataset
  benchmark, scDblFinder beats DoubletFinder/scrublet
  (https://doi.org/10.12688/f1000research.73600.2)
- Heumos et al. 2023, Nat Rev Genet — Single-cell Best Practices book,
  doublet-detection section names scDblFinder as the top performer
  (https://sc-best-practices.org/preprocessing_visualization/quality_control.html)

**Method papers:**
- McGinnis et al. 2019, Cell Systems — DoubletFinder
  (https://doi.org/10.1016/j.cels.2019.03.003)
- Wolock et al. 2019, Cell Systems — scrublet
  (https://doi.org/10.1016/j.cels.2018.11.005)
- Bernstein et al. 2020, Cell Systems — Solo
  (https://doi.org/10.1016/j.cels.2020.05.010)
- Bais & Kostka 2020, Bioinformatics — scds (cxds, bcds, hybrid)
  (https://doi.org/10.1093/bioinformatics/btz698)

**Vendor / protocol references:**
- 10x Genomics Chromium Next GEM Single Cell 3' v3.1 user guide (CG000315
  Rev E) — multiplet rate table by cells loaded/recovered
  (https://www.10xgenomics.com/support/single-cell-gene-expression)

**Recipes that execute paths from §1:**
- `bp-quality-control` — best-practice scRNA-seq QC; calls scrublet (Python)
  per sample and notes scDblFinder as the R alternative; the executable
  point of entry for any path in §1.
- `# TODO(recipe): scdblfinder-r` — dedicated R/Bioconductor scDblFinder
  recipe (currently only mentioned inline in `bp-quality-control`).
- `# TODO(recipe): doubletfinder-r` — dedicated DoubletFinder recipe.
- `# TODO(recipe): solo-doublet-scvi` — Solo via `scvi.external.SOLO`.
- `# TODO(recipe): htodemux-r` — hashtag-based doublet calling for
  multiplexed experiments.

**Adjacent knowhow:**
- `scrna_pipeline.md` — upstream framing: keep samples separate, run QC
  per-sample (and therefore call doublets per-sample).
- `scrna-de-methodology` — downstream framing: doublet flags should be
  excluded before pseudobulking; pseudobulk aggregation hides the
  contaminant otherwise.

**When to leave this knowhow altogether:**
- If the experiment is designed-doublet (PIC-seq etc.), this knowhow does
  not apply — see the protocol's own analysis guidance.
- If the data is plate-based (Smart-seq2/3, MARS-seq), doublets are not the
  dominant noise source — focus on per-well visual QC at sort time.
