---
name: fetch-geo-processed-matrices
description: Download processed/supplementary files (count matrices) for a GEO series or sample — 10x mtx triplets, h5/h5ad, or per-sample count tables — and load them for analysis.
when_to_use: You have a KNOWN GEO accession (GSE… or GSM…) and want either (a) its already-processed expression data — count matrices, not raw FASTQ — or (b) just to LIST the study's samples (GSM accessions) + per-sample metadata/characteristics. For a known accession use THIS recipe (it has the GEOparse sample table), not query-geo (which only searches for studies and can't list a series' real per-sample roster). Fast path for scRNA-seq / bulk RNA-seq when authors deposited matrices; try it BEFORE the FASTQ/realignment path.
requires_tools: [run_python]
capabilities_needed: [GEOparse, pandas, scanpy, anndata]
keywords: [GEO, GSE, GSM, count matrix, supplementary files, 10x, mtx, barcodes, features, h5ad, h5, processed data, expression matrix, scRNA-seq, bulk RNA-seq, download]
produces: [downloaded supplementary files on disk, loaded AnnData / count matrix, sample metadata table]
domain: database
source: "ABA original (2026 acquisition research) — GEOparse, scanpy"
---
# Fetch GEO processed count matrices

Most depositors upload **processed** expression data as GEO *supplementary files*,
either at the **series** level (one bundle for the study) or the **sample** level
(one file per GSM). For scRNA-seq this is almost always there; downloading and
loading it is far faster and more reliable than re-fetching FASTQ and re-aligning.
**Try this recipe first.** Fall back to `fetch-sequencing-fastq` only if the
series has no usable processed matrices, or the user explicitly needs raw reads.

**Match the deliverable to the ask — downloading is NOT the same as converting.**
If the user asked to *download* (and/or *register*) the data, fetch the
supplementary files, register them as a dataset, and stop — the files are the
deliverable. Load into **AnnData only when the NEXT step is actual analysis**
(clustering/DE). And **don't merge multiple samples into one object unless asked**:
"register them together" means one dataset entity spanning the files, not a single
merged `.h5ad`. Converting/merging by default discards the raw files the user
wanted and bakes in choices they didn't make.

`run_python` only; `ensure_capability("GEOparse")` and (`scanpy`/`anndata` for
loading 10x/h5ad, `pandas` for plain tables).

## Decision point — what kind of processed file is here?

GEO supplementary layout is not standardized. Before loading, **list the files
first** and branch on what you see:

- **10x triplet** (`*matrix.mtx.gz` + `*barcodes.tsv.gz` + `*features.tsv.gz`/`*genes.tsv.gz`):
  files are usually loose and prefixed per-sample (`GSM..._matrix.mtx.gz`), so
  `sc.read_10x_mtx` (which needs a standard CellRanger dir) fails — read the three
  parts EXPLICITLY with `sc.read_mtx` + pandas (Step 3a).
- **HDF5** (`*.h5ad`, `*.h5`, `*filtered_feature_bc_matrix.h5`): `sc.read_h5ad`
  or `sc.read_10x_h5`. Single richest case.
- **Flat table** (`*.csv.gz`, `*.txt.gz`, `*.tsv.gz`, `*counts*`): `pandas.read_csv`
  (genes × samples or cells). Watch the separator and the orientation.
- **RDS / Seurat / loom**: harder; note it and prefer the FASTQ path or ask the user.

## Step 1 — list supplementary files (cheap, no big download)

**Check the GSM (sample) level FIRST.** A single sample almost always has its own
supplementary files — for scRNA-seq that's the 10x triplet (`*barcodes.tsv.gz`,
`*features.tsv.gz`/`*genes.tsv.gz`, `*matrix.mtx.gz`), typically only tens of MB.
**Download those directly. Do NOT jump to the series `GSE…_RAW.tar`** (often many
GB — it bundles every sample) unless the GSM genuinely lists no files.

Reliable, GEOparse-free listing — parses the SOFT text directly, so it works for a
GSM *or* a GSE even when GEOparse raises (it can, on some records):

```python
import urllib.request, re
def geo_supp_files(acc):
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
           f"?acc={acc}&targ=self&form=text&view=quick")
    txt = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "replace")
    # GSM lines are !Sample_supplementary_file_N ; GSE lines are !Series_supplementary_file_N.
    # Keep the native ftp:// scheme — see Step 2's comment for why we no longer
    # rewrite it to https://.
    return re.findall(r"^![A-Za-z]+_supplementary_file_\d+\s*=\s*(\S+)", txt, re.M)

acc = "GSM5746259"                       # GSM… (one sample) or GSE… (series)
files = geo_supp_files(acc)
print(f"{acc}: {len(files)} supplementary file(s)"); [print(" ", u) for u in files]
```

(GEOparse convenience alternative: `g = GEOparse.get_GEO("GSM…", destdir="./geo_meta",
silent=True); g.metadata.get("supplementary_file")`. Handy, but it raises on some
records — the urllib method above is the dependable primary.)

If the GSM lists files (the common case) → download them directly (Step 2). Only if
it lists **none** → resolve the parent series (`!Sample_series_id`) and grab the
series supplementary, falling back to `GSE…_RAW.tar` **only as a last resort** —
and then stream it and extract just this sample's members, not the whole archive.

## Step 2 — download the listed files to disk (parallel, ftp://)

```python
import os, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed
os.makedirs("./geo_data", exist_ok=True)

# Transport choice — measured 2026-06 on ftp.ncbi.nlm.nih.gov:
#   ftp://  6-way parallel:  ~13s for a 2-sample (83 MB) 10x triplet
#   https:// 6-way parallel: ~55s for the same payload
# NCBI throttles HTTPS to ~400 KB/s per connection under contention; the
# ftp:// frontend doesn't (~6 MB/s sustained). The historical "ftp is slow"
# advice is no longer true here. KEEP the native ftp:// URLs Step 1 returned;
# fall back to https:// PER FILE if ftp:// errors (covers FTP server outages
# without sinking the whole batch).
#
# No `-C -` (resume) on the first attempt: a partial file left by a cancelled
# prior run triggers `HTTP 416 → curl rc=33`. Reserve `-C -` for an explicit
# retry branch (or just rm -f the partial and re-download).

PER_FILE_TIMEOUT = 120   # a 50 MB file at NCBI's slowest measured rate is <60s;
                         # 120 is forgiving but fails fast on a genuine stall
                         # instead of sitting silent for 10 minutes.

def _fetch(url, timeout=PER_FILE_TIMEOUT):
    """Returns (dst, ok, dt, err). Tries ftp:// first; on any failure falls
    back to https:// for that one file."""
    dst = os.path.join("./geo_data", os.path.basename(url))
    t0 = time.time()
    r = subprocess.run(["curl", "-fsS", "--max-time", str(timeout), "-o", dst, url],
                       capture_output=True, text=True)
    if r.returncode == 0:
        return dst, True, time.time() - t0, ""
    # ftp:// failed — try https:// fallback for this one URL
    url_https = url.replace("ftp://ftp.ncbi.nlm.nih.gov",
                            "https://ftp.ncbi.nlm.nih.gov")
    if url_https != url:
        r2 = subprocess.run(["curl", "-fsSL", "--max-time", str(timeout * 2),
                             "-o", dst, url_https], capture_output=True, text=True)
        if r2.returncode == 0:
            return dst, True, time.time() - t0, "(https fallback)"
        err = r2.stderr.strip()
    else:
        err = r.stderr.strip()
    return dst, False, time.time() - t0, err[:200]

# Print the plan UP FRONT so a slow start is obvious immediately, rather than
# silence followed by a single dump at the end.
print(f"Fetching {len(files)} file(s) into ./geo_data — expect ~"
      f"{max(2, len(files)*2)}-15s over ftp://ftp.ncbi.nlm.nih.gov.", flush=True)

with ThreadPoolExecutor(max_workers=min(6, len(files))) as ex:
    futs = {ex.submit(_fetch, u): u for u in files}
    done, total_bytes = 0, 0
    for fut in as_completed(futs):
        dst, ok, dt, msg = fut.result()
        done += 1
        sz = os.path.getsize(dst) if os.path.exists(dst) else 0
        total_bytes += sz
        mark = "✓" if ok else "✗"
        suffix = f"  {msg}" if msg else ""
        print(f"  [{done}/{len(files)}] {mark} {os.path.basename(dst)}  "
              f"{sz/1e6:.1f} MB in {dt:.1f}s{suffix}", flush=True)
        if not ok:
            raise RuntimeError(f"download failed for {dst}: {msg}")

print(f"Done — {total_bytes/1e6:.1f} MB total.", flush=True)
```

Per-sample 10x triplets are 30–100 MB (`matrix.mtx.gz` dominates); the parallel
ftp:// fetch above completes in tens of seconds even on a household connection.
**For larger archives** — `GSE…_RAW.tar` is often multi-GB — run as a background
job, write to a durable path, stream (`requests.get(stream=True)` +
`iter_content`) rather than buffering in the kernel, and `tar.extract()` only
the members whose names start with your `GSM…` prefix. That's the *one*
legitimate place to use `curl -C -` (resume on an interrupted multi-GB pull). Verify each file exists and is
non-empty before loading.

## Step 3 — load by branch

### 3a. 10x mtx triplet
GEO supplementary triplets are **loose, GSM-prefixed** (`GSM..._matrix.mtx.gz` /
`...barcodes.tsv.gz` / `...features.tsv.gz` all in one dir) — `sc.read_10x_mtx`
will NOT find these (non-standard names), so read the three parts EXPLICITLY.
This is the usual GEO layout and the #1 source of "lost gene names" flailing:

```python
import scanpy as sc, pandas as pd, anndata as ad, os
D = "./geo_data"
def load_geo_10x(prefix):                       # one GEO loose, GSM-prefixed triplet
    a = sc.read_mtx(f"{D}/{prefix}.matrix.mtx.gz").T          # mtx is genes×cells → transpose
    a.obs_names = pd.read_csv(f"{D}/{prefix}.barcodes.tsv.gz", header=None)[0].values
    a.var_names = pd.read_csv(f"{D}/{prefix}.features.tsv.gz", header=None, sep='\t')[1].values  # col 2 = symbols
    a.var_names_make_unique(); return a         # do NOT skip — duplicate symbols are common

# single sample:
adata = load_geo_10x("GSM5354513_...")          # the shared file prefix (sans .matrix.mtx.gz etc.)

# MULTIPLE samples — concat with a batch key (one object spanning the samples,
# only when the next step is batch-aware analysis; see the merge caveat above):
prefixes = ["GSM5354513_...", "GSM5354514_..."]
sample_names = ["S1", "S2"]
adata = ad.concat([load_geo_10x(p) for p in prefixes],
                  label='sample', keys=sample_names, index_unique='-')
```

If the files happen to be in a **standard CellRanger dir** (canonical
`matrix.mtx.gz`/`barcodes.tsv.gz`/`features.tsv.gz` names), use the loader instead:
`adata = sc.read_10x_mtx(dir, var_names='gene_symbols')`.

### 3b. HDF5
```python
import scanpy as sc
adata = sc.read_10x_h5("./geo_data/GSM.../filtered_feature_bc_matrix.h5")  # cellranger h5
# or: adata = sc.read_h5ad("./geo_data/.../something.h5ad")
```

### 3c. flat table
```python
import pandas as pd
df = pd.read_csv("./geo_data/GSM.../GSE..._counts.csv.gz", index_col=0)
# Decide orientation: genes are usually rows. If columns look like ENSG/symbols, transpose.
import anndata as ad
adata = ad.AnnData(df.T)   # AnnData wants cells/samples × genes
```

## Step 4 — attach sample metadata

```python
import pandas as pd
meta = pd.DataFrame({g: s.metadata for g, s in gse.gsms.items()}).T
# characteristics_ch1 holds the useful condition/tissue/genotype annotations.
```

## Key decisions
- **Processed-first**: this recipe is the default for "get the data from GSE…".
  Only go to FASTQ if matrices are absent/unusable or reads are truly needed.
- `download_sra=False` always — keep SRA out of this path.
- Per-sample vs series-level files: scRNA-seq is usually per-GSM 10x triplets;
  bulk is often one series-level count table.

## Gotchas
- **Just need a cell count?** Download only `barcodes.tsv.gz` (one barcode = one
  cell) and count its lines — don't pull the (often 50–150 MB) `matrix.mtx.gz`.
- **No naming standard.** You must list files and branch; do not assume 10x.
- 10x triplets are commonly prefixed per sample and need regrouping (Step 3a).
- `features.tsv.gz` vs legacy `genes.tsv.gz`; CellRanger v2 vs v3 layouts differ.
- Flat tables: check separator (`\t` vs `,`) and orientation before trusting it.
- **Use `ftp://` for NCBI supplementary files, not `https://`.** Measured
  2026-06: ftp:// to `ftp.ncbi.nlm.nih.gov` runs ~6 MB/s sustained even on a
  household connection; HTTPS to the same host throttles each connection to
  ~400 KB/s under contention. A 6-way parallel 2-sample 10x fetch is ~13s
  over ftp vs ~55s over https. The old "ftp:// is slow/flaky" lore was
  written against a different NCBI; today it's inverted. Step 2 keeps ftp://
  and falls back to https:// per-file if a single ftp transfer errors.
- **Don't use `curl -C -` on a fresh download.** Triggers `HTTP 416 → curl rc=33`
  when a partial file from a cancelled prior run is on disk. Reserve `-C -` for
  an explicit retry on the multi-GB `_RAW.tar` path only.
- **Parallelize AND print as files complete.** Six silent curls behind a join
  look identical to a hang. Step 2's example uses `as_completed` + per-file
  `print(..., flush=True)` so the user sees motion immediately (the small
  `barcodes.tsv.gz` lands in <1s and gives a heartbeat while the big
  `matrix.mtx.gz` is still streaming). If you write your own loop, do the
  same — never bury a parallel fetch behind `ex.map` + a single end-of-batch
  log.
- **Per-file timeout should be tight (≤120s for a 50 MB file).** Long
  timeouts (`--max-time 600`) hide stalls behind silence; short ones surface
  them as fast failures the recipe can retry/fall back on.
- Some series only deposit RDS/Seurat objects → note it; FASTQ path may be cleaner.
- `GEOparse` caches SOFT files in `destdir`; reuse it to avoid refetching metadata.

## In ABA
`run_python`; `ensure_capability("GEOparse")`, `ensure_capability("scanpy")`
(pulls anndata), `ensure_capability("pandas")`. Big file pulls → background/streamed
job, durable disk path, verify size before load.
