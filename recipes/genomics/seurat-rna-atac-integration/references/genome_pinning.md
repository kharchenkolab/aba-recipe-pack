# Genome / build pinning — EnsDb, BSgenome, seqlevels style

Picking the right EnsDb / BSgenome combination per organism + Cell
Ranger ARC reference build (hg38 vs hg19 vs mm10), seqlevels style
mapping (Ensembl `1` vs UCSC `chr1`), and why a mismatch silently
breaks downstream calls. Load this when the agent needs to swap to a
non-human organism, debug an empty GeneActivity / failed
CoveragePlot, or confirm the annotation matches the FASTQs' alignment
reference.

Sources consulted: Bioconductor `EnsDb.*` and `BSgenome.*` package
descriptions; Cell Ranger ARC documentation
(`support.10xgenomics.com/single-cell-multiome-atac-gex/software/`);
GenomeInfoDb / GenomicRanges vignettes; Signac's
`GetGRangesFromEnsDb` source.

## What needs to match

Three things must agree on the SAME genome build:

1. **The FASTQs' alignment.** Set when Cell Ranger ARC was run — the
   reference package used. For human, this is one of:
   - `refdata-cellranger-arc-GRCh38-2020-A-2.0.0` (hg38, Ensembl 98)
   - `refdata-cellranger-arc-GRCh38-2020-A` (hg38, Ensembl 98 — older
     packaging)
   Cell Ranger ARC has not (as of the recipe's source vignette)
   released an hg19 reference; hg19 multiome is rare.
2. **The EnsDb annotation** used in `GetGRangesFromEnsDb`. Must match
   the same build (hg38 → `EnsDb.Hsapiens.v86`).
3. **The BSgenome** used for `RegionStats` / `LinkPeaks`. Must match
   (`BSgenome.Hsapiens.UCSC.hg38`).

If any of the three is on a different build, peak coordinates,
annotation intervals, and sequence lookups will silently disagree.

## The canonical pairings

### Human

| FASTQ alignment | EnsDb | BSgenome | seqlevels style |
|---|---|---|---|
| hg38 (GRCh38) — default | `EnsDb.Hsapiens.v86` | `BSgenome.Hsapiens.UCSC.hg38` | UCSC (`chr1`) |
| hg19 (GRCh37) — legacy | `EnsDb.Hsapiens.v75` | `BSgenome.Hsapiens.UCSC.hg19` | UCSC (`chr1`) |

`EnsDb.Hsapiens.v86` is Ensembl release 86 = GRCh38. Newer EnsDb
versions exist (v97, v109, …) but v86 is the one Signac vignettes
test against; later versions are usually compatible (gene model
updates, not coordinate changes).

### Mouse

| FASTQ alignment | EnsDb | BSgenome | seqlevels style |
|---|---|---|---|
| mm10 (GRCm38) | `EnsDb.Mmusculus.v79` | `BSgenome.Mmusculus.UCSC.mm10` | UCSC (`chr1`) |
| mm9 (NCBIm37) | `EnsDb.Mmusculus.v75` | `BSgenome.Mmusculus.UCSC.mm9` | UCSC (`chr1`) |

Cell Ranger ARC's default mouse reference is mm10
(`refdata-cellranger-arc-mm10-2020-A-2.0.0`).

The recipe's Step 2 mitochondrial pattern also changes by organism:

```r
# Human
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Mouse
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
```

Mouse genes are typically lower-case (e.g. `mt-Co1`), human genes
upper-case (`MT-CO1`).

### Other organisms

Bioconductor has EnsDb + BSgenome packages for many organisms. Find
the canonical names:

```r
# All EnsDb packages installed (or available via BiocManager::available())
BiocManager::available("EnsDb")
BiocManager::available("BSgenome")

# Or browse the catalogues:
# https://bioconductor.org/packages/release/data/annotation/
# (filter by "EnsDb" or "BSgenome")
```

Common ones:

| Organism | EnsDb | BSgenome |
|---|---|---|
| Rat | `EnsDb.Rnorvegicus.v79` | `BSgenome.Rnorvegicus.UCSC.rn6` |
| Zebrafish | `EnsDb.Drerio.v79` | `BSgenome.Drerio.UCSC.danRer11` |
| Drosophila | `EnsDb.Dmelanogaster.v79` | `BSgenome.Dmelanogaster.UCSC.dm6` |

For non-model organisms (the EnsDb / BSgenome pair doesn't exist),
you'll need to build the annotation manually — `ensembldb` has tools
for this; see `?makeEnsembldbPackage`.

## Seqlevels style — the silent killer

Ensembl chromosome names: `1`, `2`, …, `X`, `Y`, `MT`.
UCSC chromosome names: `chr1`, `chr2`, …, `chrX`, `chrY`, `chrM`.

EnsDb returns Ensembl style by default. Cell Ranger ARC writes peaks
in UCSC style. They DO NOT MATCH WITHOUT TRANSLATION.

```r
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
print(head(seqlevels(annotation)))            # "1" "2" ... "Y" "MT"

seqlevelsStyle(annotation) <- "UCSC"
print(head(seqlevels(annotation)))            # "chr1" "chr2" ... "chrY" "chrM"
```

The translation:

| Ensembl | UCSC |
|---|---|
| `1` | `chr1` |
| `2` | `chr2` |
| `X` | `chrX` |
| `Y` | `chrY` |
| `MT` | `chrM` |

Note `MT` ↔ `chrM` — the mitochondrial chromosome name changes letter,
not just prefix. The `seqlevelsStyle` setter handles this
automatically.

### How to spot a seqlevels mismatch

If `GeneActivity` returns an empty matrix:

```r
ga <- GeneActivity(obj)
dim(ga)                     # 0 x N — empty
```

The most likely cause: peak GRanges (UCSC `chr1`) doesn't overlap
annotation GRanges (Ensembl `1`). Verify:

```r
peaks_seq <- seqlevels(granges(obj[["ATAC"]]))
anno_seq  <- seqlevels(Annotation(obj[["ATAC"]]))
intersect(peaks_seq, anno_seq)            # should be ~24; if 0, mismatch
```

Fix: rebuild the annotation with UCSC style and re-set
`Annotation(obj[["ATAC"]]) <- annotation`.

### How to spot a build mismatch (hg38 vs hg19)

Even with matching seqlevels style, the COORDINATES change between
builds. A peak at `chr11:60,455,000-60,455,500` in hg38 lands in a
different gene in hg19. Symptoms:

- CoveragePlot at a known marker (e.g. `MS4A1`) shows a peak in an
  unexpected location.
- Markers don't match between RNA and GeneActivity in CoveragePlot
  panels.

Verify the build:

```r
unique(genome(annotation))           # should be "hg38" or "GRCh38"
unique(genome(granges(obj[["ATAC"]])))   # should match
```

If they disagree, you used the wrong EnsDb. Fix is to start over
with the right one.

## EnsDb installation

EnsDb packages are on Bioconductor's "AnnotationHub" tier — install
with BiocManager:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("EnsDb.Hsapiens.v86")
```

For lab-shared environments, prefer the bioconda channel:

```bash
mamba install -c bioconda bioconductor-ensdb.hsapiens.v86
```

The package is ~150 MB; install once per environment.

## BSgenome installation

BSgenome packages are much larger (genome sequence is bundled — ~1
GB for human). Skip if you don't need `RegionStats` / `LinkPeaks`.

```r
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

bioconda:

```bash
mamba install -c bioconda bioconductor-bsgenome.hsapiens.ucsc.hg38
```

If the only reason you'd install BSgenome is one CoveragePlot call,
skip it — CoveragePlot doesn't need it (it reads fragment positions
from the tabix-indexed `.tsv.gz`, not sequence). BSgenome is ONLY
needed for:

- `RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)` — computes
  GC content per peak (prerequisite for LinkPeaks).
- `LinkPeaks` (uses the GC content from RegionStats for the
  matched-control resampling).
- Motif analysis (`AddMotifs` + `RunChromVAR`) — extracts sequence
  under peaks.

## Cell Ranger ARC reference identification

The reference build is encoded in the Cell Ranger ARC run's outputs.
The truth source is the `summary.csv` or the parameter file:

```
sample_name/
├── outs/
│   ├── filtered_feature_bc_matrix.h5
│   ├── atac_fragments.tsv.gz
│   ├── atac_fragments.tsv.gz.tbi
│   └── summary.csv               ← grep for "Reference"
└── _outs/
```

`summary.csv` has a row like:

```
Reference,refdata-cellranger-arc-GRCh38-2020-A-2.0.0
```

The build name tells you everything:
- `GRCh38` → hg38 → `EnsDb.Hsapiens.v86` + `BSgenome.Hsapiens.UCSC.hg38`
- `GRCm38` → mm10 → `EnsDb.Mmusculus.v79` + `BSgenome.Mmusculus.UCSC.mm10`

If the user can't tell you, ask. Don't guess — wrong-build analyses
are the most common cause of "this multiome looks completely wrong"
reports.
