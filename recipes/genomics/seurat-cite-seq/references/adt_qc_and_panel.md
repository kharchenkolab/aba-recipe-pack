# ADT panel QC — suffix cleanup, isotype controls, per-protein gates

Deep detail behind Step 1 (panel-suffix cleanup) and Step 4 (per-protein
VlnPlot diagnosis) of `seurat-cite-seq`. Loaded when the agent needs to
audit a panel — clean up vendor name suffixes, identify isotype/hashtag
rows that shouldn't enter analysis, or diagnose a VlnPlot that shows
flat/diffuse signal.

## Panel-suffix cleanup — what the regex menu covers

Antibody panels from commercial vendors ship with technology suffixes
glued onto the symbol. They look like `CD3_TotalSeqB` or
`MouseIgG2a_control_TotalSeqB`. Strip them so the recipe's downstream
`adt_<protein>` keys are clean:

```r
# Standard 10x + BioLegend TotalSeq-B
rownames(adt_counts) <- gsub("_[A-Za-z_]*TotalSeq[A-C]?$", "",
                             rownames(adt_counts))
```

The character class `[A-Za-z_]*` covers the `_control_` insert that
isotype antibodies carry (e.g. `MouseIgG1_control_TotalSeqB`). The
`TotalSeq[A-C]?$` tail matches the three vendor variants (TotalSeq-A,
TotalSeq-B, TotalSeq-C — all three are shipped by BioLegend with the
same suffix shape).

Other common suffixes you might encounter:

| Vendor / chemistry | Suffix shape | Regex |
|---|---|---|
| BioLegend TotalSeq-A/B/C | `_TotalSeqB`, `_control_TotalSeqB` | `_[A-Za-z_]*TotalSeq[A-C]?$` |
| BD AbSeq | `_AbSeq` | `_AbSeq$` |
| Custom in-house | varies — inspect first | `<custom>` |
| Hashtag oligos (HTO) | usually a separate "Multiplexing Capture" feature_type, NOT under `Antibody Capture`. If you see `HTO_*` in `Antibody Capture`, route them out to a separate assay before normalization (`obj[["HTO"]]` rather than `obj[["ADT"]]`). |

**Audit the rownames BEFORE downstream code.** A suffix you missed will
silently bleed into every plot title and every Step 6 protein-to-RNA
pair lookup:

```r
cat("ADT rownames sample (first 10):\n")
print(head(rownames(adt_counts), 10))
```

If suffixes remain, extend the `gsub` and re-run.

## Isotype controls — what they are and what to do with them

Isotype controls are antibodies raised against a target the cells
don't express (typically mouse IgG isotypes on a human sample). They
measure non-specific binding — the noise floor each cell carries.
Panels designed for DSB normalization include 2–4 isotypes; many
published panels include none.

How to find them in a panel:

```r
adt_features <- rownames(obj[["ADT"]])
isotype_patterns <- c("IgG[0-9]*[ab]?$", "isotype", "ISO_", "[Cc]ontrol",
                      "Rat[-_]?IgG", "Mouse[-_]?IgG", "Armenian[-_]?Hamster")
isotype_idx <- which(Reduce(`|`, lapply(isotype_patterns, grepl, adt_features)))
cat("Suspected isotypes:\n")
print(adt_features[isotype_idx])
```

What to do with them:

1. **Keep them in the assay during CLR normalization** — they participate
   in the per-feature geometric mean (with `margin=2`) but on a panel
   of 25 markers + 3 isotypes the effect is minor.
2. **Exclude them from clustering** — `VariableFeatures(obj) <-
   setdiff(rownames(obj[["ADT"]]), isotypes)` before any ADT-side
   `RunPCA`. Including isotypes adds rank to the PCA without adding
   biology.
3. **Plot them on the QC VlnPlot (Step 4)** — they should be uniformly
   low across all clusters. If an isotype spikes in one cluster,
   that's evidence of non-specific binding bias (cell-size, viability,
   FcR expression).
4. **Use them for DSB** — `isotype.control.name.vec` argument; see
   `references/adt_normalization.md`.

## Hashtag oligos (HTO) — not ADT

If the panel multiplexes donors via cell hashtags, those usually arrive
as a third feature_type:

- 10x with cell-multiplexing: `"Multiplexing Capture"`
- Some pipelines park them under `"Antibody Capture"` as `HTO_*`
  rownames.

The downstream call `HTODemux` (Seurat) expects them in their **own**
assay, not the ADT assay:

```r
# Pull HTO rows OUT of ADT before normalization.
hto_idx <- grep("^HTO_", rownames(adt_counts))
if (length(hto_idx) > 0) {
  hto_counts <- adt_counts[hto_idx, , drop = FALSE]
  adt_counts <- adt_counts[-hto_idx, , drop = FALSE]
  obj[["HTO"]] <- CreateAssayObject(counts = hto_counts[, colnames(obj)])
}
```

This recipe assumes HTO de-multiplexing happened upstream — if the
input matrix still carries HTO rows, route them out before Step 3.

## Per-protein QC patterns and what they mean

Step 4 plots a per-protein VlnPlot grouped by RNA cluster. Three
patterns to recognize:

| Pattern | What it looks like | Likely cause | What to do |
|---|---|---|---|
| **Lineage marker — clean** | Tall violins on 1–2 clusters, flat-zero on the rest | Antibody works, lineage is captured | Keep — use for cluster annotation |
| **Pan-marker — clean** | Tall violins on every cluster (e.g. CD45 on all hematopoietic cells) | Antibody works, target is broadly expressed | Keep — useful as a viability check, not for cluster separation |
| **Diffuse / no signal** | Flat distributions across all clusters near the floor | Failed antibody, dropped during staining, OR symbol misread (e.g. `CD8` vs `CD8a` panel naming) | Audit symbol — `rownames(obj[["ADT"]])` — first. If symbol is right, drop the protein from analysis |
| **Background spike** | Flat distributions but with a long tail in one cluster | Non-specific binding bias (FcR, cell-size, viability) | If panel has isotypes, compare to isotype's same-cluster pattern. If isotypes spike too, the cluster has an FcR / size confound — exclude from cluster annotation, not from QC |
| **Bimodal across clusters** | Each cluster shows a clear two-peak distribution | Cell-state heterogeneity within the cluster — protein is dynamic | Don't "fix" — this is biology. The protein splits the cluster; consider WNN (`seurat-wnn-multimodal`) for the joint workup |

## Common symbol mismatches

ADT panels label antibodies by the protein epitope (often the surface
glycoform); RNA names follow HGNC. The mapping is rarely 1:1:

| ADT panel name | RNA HGNC symbol | Note |
|---|---|---|
| `CD3` | `CD3E`, `CD3D`, or `CD3G` | T-cell receptor complex — usually `CD3E` for the canonical T-cell marker |
| `CD4` | `CD4` | matches |
| `CD8`, `CD8a` | `CD8A` | bmcite uses `CD8a` (lowercase a); 10x panels typically `CD8`. Case-sensitive |
| `CD19` | `CD19` | matches; for B-cell ID, the RNA `MS4A1` is often a better marker than `CD19` itself |
| `CD56`, `NCAM1` | `NCAM1` | older panels use the CD name, newer use the gene symbol |
| `CD16` | `FCGR3A`, `FCGR3B` | the protein FcγRIII has two genes |
| `CD20` | `MS4A1` | the antibody binds MS4A1's surface form |
| `CD25` | `IL2RA` | the IL-2 receptor alpha chain |
| `HLA-DR` | `HLA-DRA`, `HLA-DRB1`, etc. | the antibody recognizes the heterodimer; multiple genes |

For Step 6's ADT-vs-RNA pair plot, the protein-to-gene map matters.
**Build the map per panel, don't assume it.** When in doubt, ask the
user or look up the antibody clone in the vendor datasheet.

## Panel-size selection criteria — when to skip ADT-only steps

The SKILL.md's Decision #2 says "small panels (<10) should skip
ADT-only clustering entirely". Concrete criteria:

- **<10 markers** — overlay only. Don't `RunPCA` on ADT. Don't try WNN.
  The protein contribution to a joint graph is noise-dominated.
- **10–20 markers, lineage-focused** — overlay + biaxials are
  diagnostic. WNN is possible but rarely wins over RNA-only.
- **20–50 markers, lineage + activation** — full workup; this is the
  WNN regime. CITE-seq vignette uses 13 markers (cbmc) and 25 markers
  (bmcite) and gets meaningful joint clusters at the 25-marker end.
- **>50 markers (custom panels)** — DSB if isotypes are in. Long-tail
  panels often have protein-only sub-populations the RNA misses.

## Audit script — run after Step 1

```r
adt <- rownames(obj[["ADT"]])
cat(sprintf("ADT panel: %d markers\n", length(adt)))
cat("Suspected isotypes:\n")
iso <- grep("IgG|isotype|control", adt, ignore.case = TRUE, value = TRUE)
print(iso)
cat("Marker names sample:\n")
print(head(setdiff(adt, iso), 15))
if (length(adt) - length(iso) < 10) {
  cat("WARNING: <10 informative markers; treat ADT as overlay only.\n")
}
```

## References

- Stoeckius M., et al. (2018). *Cell Hashing with barcoded antibodies
  enables multiplexing and doublet detection for single cell genomics.*
  Genome Biol. 19, 224 — the HTO demultiplexing paper, source of the
  HTO-as-separate-assay pattern.
- Hao Y., et al. (2021). *Integrated analysis of multimodal single-cell
  data.* Cell 184, 3573–3587 — bmcite + the 25-marker hematopoietic
  panel used in the WNN vignette.
- BioLegend TotalSeq antibody documentation
  (<https://www.biolegend.com/en-us/totalseq>) — authoritative source
  for the `_TotalSeq[A-C]` suffix conventions.
