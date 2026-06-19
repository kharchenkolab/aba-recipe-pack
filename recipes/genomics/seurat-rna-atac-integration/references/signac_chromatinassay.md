# Signac `ChromatinAssay` — internals

What a Signac `ChromatinAssay` is, how `CreateChromatinAssay` consumes
peaks + fragments + annotation, the GRanges interplay, the soft-dep on
`biovizBase`, and the Cell Ranger ARC multimodal-H5 layout the recipe
reads from. Load this when the agent needs to debug an empty
GeneActivity matrix, a missing-feature error from CoveragePlot, or any
unexpected behavior tied to the fragments file.

Sources consulted: Stuart et al. *Nat Methods* 18, 1333–1341 (2021) —
"Single-cell chromatin state analysis with Signac" (the ChromatinAssay
design paper); Signac source `R/objects.R` (constructor) and
`R/fragments.R` (fragment-lookup machinery), tested against Signac
1.16.0; Signac vignettes at `stuartlab.org/signac/articles/`.

## Class structure

A `ChromatinAssay` is an S4 class extending Seurat's `Assay`. The
class adds:

- `@ranges` — a `GRanges` object: one interval per row (peak) of the
  counts matrix. Parsed from the rownames of the input counts via the
  `sep` argument.
- `@fragments` — a list of `Fragment` objects, each wrapping a path to
  a bgzipped fragments file (`.tsv.gz`) plus its tabix index (`.tbi`).
- `@annotation` — a `GRanges` of gene/transcript intervals (typically
  from `GetGRangesFromEnsDb`), used by `GeneActivity`, `CoveragePlot`,
  and `LinkPeaks` to relate peaks to genes.
- `@genome` — a `Seqinfo` (or character genome name) declaring the
  reference build. `seqlevelsStyle` matters here — see below.
- `@motifs` — optional, populated by `AddMotifs`.
- `@positionEnrichment` — optional, populated by `TSSEnrichment`.
- Standard `Assay` slots (counts, data, scale.data, var.features,
  meta.features) for the peak × cell matrix itself.

You don't manipulate these slots directly in normal use; access via
`granges(obj[["ATAC"]])`, `Fragments(obj[["ATAC"]])`,
`Annotation(obj[["ATAC"]])`.

## `CreateChromatinAssay` — argument-by-argument

```r
obj[["ATAC"]] <- CreateChromatinAssay(
  counts     = counts$Peaks,                       # peak x cell sparse matrix (dgCMatrix)
  sep        = c(":", "-"),                        # parses "chr1:100-200" -> GRanges("chr1", 100-200)
  fragments  = "/path/to/atac_fragments.tsv.gz",   # PATH (not the .tbi)
  annotation = annotation,                         # GRanges from GetGRangesFromEnsDb
  genome     = "hg38"                              # OR Seqinfo() with seqlengths
)
```

Each argument:

- **`counts`** — `dgCMatrix` (peak × cell), integer counts.
  `Read10X_h5` on a multiome H5 returns a list; the ATAC slot is named
  `Peaks` (with a capital P), NOT `ATAC`.
- **`sep`** — the two delimiters in the peak rowname. Cell Ranger ARC
  writes `chr1:100-200`, so `sep = c(":", "-")`. Older outputs may use
  `chr1-100-200`, in which case `sep = c("-", "-")`.
- **`fragments`** — path to the bgzipped fragments file (`.tsv.gz`),
  NOT the tabix index. The index (`.tbi`) MUST live next to the
  `.tsv.gz`; Signac reads it implicitly via the `Rsamtools` tabix
  bindings. If you symlink, symlink both.
- **`annotation`** — a `GRanges` of gene/transcript intervals. The
  canonical way to build it is
  `GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)` then
  `seqlevelsStyle(annotation) <- "UCSC"`. See "seqlevels style"
  below.
- **`genome`** — character (`"hg38"`) or a `Seqinfo` object. The
  string form is shorthand; Signac uses it to label outputs but
  doesn't fetch sequence (you need a separate `BSgenome` for that).
- **`min.cells` / `min.features`** — defaults are 0 (keep
  everything). You can pre-filter rare peaks here, but the canonical
  filtering is `FindTopFeatures` after construction.

## Why `biovizBase` matters (and why it's not auto-pulled)

`GetGRangesFromEnsDb()` internally calls `biovizBase::crunch()` —
biovizBase is the package that converts EnsDb objects into the
GRanges Signac expects. But:

- `biovizBase` is a Signac `Suggests`, NOT `Imports`. The package
  installs fine without it.
- The error message ("Please install biovizBase") only appears at the
  first call to `GetGRangesFromEnsDb` — Signac's loader doesn't check
  for it at startup.

Always run `BiocManager::install("biovizBase")` once before using
ChromatinAssay. The recipe's Install block does this; if the agent
hits "Please install biovizBase", the install step was skipped.

## seqlevels style — UCSC vs Ensembl

EnsDb annotations come with Ensembl-style chromosome names: `1`, `2`,
…, `X`, `Y`, `MT`. Cell Ranger ARC outputs peaks with UCSC-style
names: `chr1`, `chr2`, …, `chrX`, `chrY`, `chrM`. The two styles do
not match without translation.

```r
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
print(head(seqlevels(annotation)))             # "1" "2" ...

seqlevelsStyle(annotation) <- "UCSC"
print(head(seqlevels(annotation)))             # "chr1" "chr2" ...
```

This is the silent-killer step. If you forget it:

- `granges(obj[["ATAC"]])` will report `chr1`/`chr2` (UCSC) but
  `Annotation(obj[["ATAC"]])` will report `1`/`2` (Ensembl).
- `GeneActivity()` returns an empty matrix (no annotation intervals
  overlap any peak).
- `CoveragePlot()` fails or shows blank panels.

The recipe documents this in Step 1's pitfalls. Set
`seqlevelsStyle(annotation) <- "UCSC"` BEFORE passing to
`CreateChromatinAssay`.

For organism / build pinning details, read
`references/genome_pinning.md`.

## Fragment lookup — what happens when CoveragePlot reads

When you call `CoveragePlot(obj, region = "MS4A1", ...)`:

1. Signac resolves `"MS4A1"` to a GRanges interval using the
   `Annotation` slot.
2. For each cluster (or grouping), Signac calls
   `Rsamtools::scanTabix` on the fragments file with that interval as
   the query, parses the fragments, and counts them.
3. The result is a per-cluster pileup track + the annotated gene
   model.

So the `.tbi` index is the load-bearing dependency for plotting; if
it's missing or stale (older than the `.tsv.gz`), Signac will rebuild
it on first call (which is slow on large samples).

If the fragments path moved (sample moved between directories), update
it with `Fragments(obj[["ATAC"]]) <- ...` — the path stored in the
ChromatinAssay is absolute.

## Multimodal H5 layout (Cell Ranger ARC)

Cell Ranger ARC writes a single H5 per sample containing BOTH RNA and
ATAC count matrices:

```
filtered_feature_bc_matrix.h5
├── matrix
│   ├── data           # nonzero values, concatenated across all features
│   ├── indices        # row indices per nonzero value
│   ├── indptr         # column-wise offsets
│   └── shape          # (n_features, n_barcodes)
├── features
│   ├── id             # feature IDs (Ensembl IDs for RNA, "chr1:100-200" for peaks)
│   ├── name           # feature symbols (gene symbols for RNA; same as id for peaks)
│   ├── feature_type   # "Gene Expression" or "Peaks"
│   ├── genome         # genome name per feature
│   └── interval       # genomic interval per feature (relevant for Peaks)
└── barcodes           # cell barcodes
```

`Read10X_h5` splits on `feature_type` and returns:

```r
counts <- Read10X_h5("/path/to/filtered_feature_bc_matrix.h5")
# counts$`Gene Expression`   # dgCMatrix (genes x cells), gene symbols as rownames
# counts$Peaks               # dgCMatrix (peaks x cells), "chr1:100-200" as rownames
```

The barcodes are SHARED — same cells across both modalities (that's
what makes it a multiome). If `ncol(counts$`Gene Expression`) !=
ncol(counts$Peaks)`, the H5 is malformed.

A single-modality H5 (RNA-only or ATAC-only) returns a single matrix,
not a list. Test with `is.list(counts)` if you're unsure.

## Inspecting a ChromatinAssay

```r
atac <- obj[["ATAC"]]
class(atac)                        # "ChromatinAssay"
nrow(atac); ncol(atac)             # peaks x cells
head(granges(atac))                # GRanges of peaks
Fragments(atac)                    # list of Fragment objects with paths
Annotation(atac)                   # GRanges of gene/transcript intervals
genome(atac)                       # "hg38" (or named Seqinfo)

# Sanity: do peaks and annotation share seqlevels?
common <- intersect(seqlevels(atac), seqlevels(Annotation(atac)))
length(common)                     # should be ~24 (autosomes + X/Y/M)
```

If `common` is 0, the seqlevels-style mismatch happened — fix the
annotation and re-`CreateChromatinAssay` (cheaper than re-running the
whole pipeline).
