---
name: fetch-geo-processed-matrices
description: Download processed/supplementary files (count matrices) for a GEO series or sample — 10x mtx triplets, h5/h5ad, or per-sample count tables — and load them for analysis.
when_to_use: You have a KNOWN GEO accession (GSE… or GSM…) and want either (a) its already-processed expression data — count matrices, not raw FASTQ — or (b) just to LIST the study's samples (GSM accessions) + per-sample metadata/characteristics. For a known accession use THIS recipe (it has the GEOparse sample table), not query-geo (which only searches for studies and can't list a series' real per-sample roster). Fast path for scRNA-seq / bulk RNA-seq when authors deposited matrices; try it BEFORE the FASTQ/realignment path.
requires_tools: [run_python, register_dataset]
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
(clustering/DE) — and on that path **still register first**: what you fetched is the
subject of the analysis, so it needs a dataset entity before you analyze, not instead
of it. Either way registration is `register_dataset`, one entity spanning the files.
And **don't merge multiple samples into one object unless asked**:
"register them together" means one dataset entity spanning the files, not a single
merged `.h5ad`. Converting/merging by default discards the raw files the user
wanted and bakes in choices they didn't make.

`run_python` for the fetch/load; `register_dataset` for the files you pulled;
`ensure_capability("GEOparse")` and (`scanpy`/`anndata` for loading 10x/h5ad,
`pandas` for plain tables).

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
    # GEO hands back ftp:// URLs. Do NOT pin a scheme here — Step 2 measures both
    # against this host and downloads over whichever is actually faster from where
    # you are running.
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

## Step 2 — download the listed files (verify every byte; probe the transport)

Two rules here, both learned the hard way.

**1. Never trust `curl`'s exit code alone.** A transfer cut mid-stream — a
`--max-time` cutoff, or an FTP data channel the server closes cleanly — leaves a
*partial* file on disk. `returncode == 0` does **not** mean "complete", and
`os.path.getsize() > 0` is not a completeness check. A truncated `matrix.mtx.gz`
sails through as ✓ and only detonates later inside `gzip` or `sc.read_mtx`, far from
the cause. Check the byte count against the server's advertised size, and
decompress-test every `.gz`. **A short file is a FAILED download** — delete it and
retry on the other scheme. This is the single most important thing in this recipe.

**2. Don't hardcode `ftp://` or `https://`.** Which is faster depends entirely on the
network you are on. Same 47 MB `matrix.mtx.gz` on `ftp.ncbi.nlm.nih.gov`, single stream:

| where | `ftp://` | `https://` |
|---|---|---|
| VBC compute node (2026-07) | 3.8 MB/s | **20.6 MB/s** |
| household connection (2026-06) | **~6 MB/s** | ~0.4 MB/s (HTTPS throttled) |

Both are real; neither generalizes. So **probe**: race the two schemes on the
smallest file (a ~37 KB `barcodes.tsv.gz` — under a second either way), pull the rest
over the winner, and fall back per-file to the loser on any failure or truncation.

```python
import os, re, gzip, time, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

DEST = "./geo_data"; os.makedirs(DEST, exist_ok=True)

# Generous: with completeness verified below, a timeout no longer silently truncates —
# it fails, deletes the partial, and retries the other scheme. Tight timeouts used to be
# a defence against silent stalls; verification is a better one.
PER_FILE_TIMEOUT = 300

def _as(url, scheme):
    """Same host + path, other scheme."""
    return re.sub(r"^(?:ftp|https)://", scheme + "://", url)

def _remote_size(url, timeout=30):
    """Server-advertised size: Content-Length (https) or SIZE (ftp). None if unknown."""
    r = subprocess.run(["curl", "-fsSIL", "--max-time", str(timeout), url],
                       capture_output=True, text=True)
    m = re.search(r"(?im)^content-length:\s*(\d+)", r.stdout) if r.returncode == 0 else None
    return int(m.group(1)) if m else None

def _complete(dst, expect):
    """Complete iff the byte count matches AND (for .gz) it inflates to the end.
    This is what catches silent truncation."""
    if not os.path.exists(dst) or os.path.getsize(dst) == 0:
        return False
    if expect is not None and os.path.getsize(dst) != expect:
        return False
    if dst.endswith(".gz"):
        try:
            with gzip.open(dst, "rb") as fh:
                while fh.read(1 << 20):
                    pass                     # raises on truncation / bad CRC
        except Exception:
            return False
    return True

def _curl(url, dst, timeout):
    # No `-C -` (resume) on a fresh download: a partial left by a cancelled prior run
    # triggers `HTTP 416 -> curl rc=33`. We overwrite, and verify afterwards.
    return subprocess.run(["curl", "-fsSL", "--max-time", str(timeout), "-o", dst, url],
                          capture_output=True, text=True).returncode == 0

def _fetch(url, scheme, timeout=PER_FILE_TIMEOUT):
    """Download over `scheme`; on failure OR truncation retry the other scheme.
    Never leaves a partial file behind. Returns (dst, ok, dt, note)."""
    dst = os.path.join(DEST, os.path.basename(url))
    t0, other = time.time(), ("https" if scheme == "ftp" else "ftp")
    for sch in (scheme, other):
        u = _as(url, sch)
        if _curl(u, dst, timeout) and _complete(dst, _remote_size(u)):
            return dst, True, time.time() - t0, ("" if sch == scheme else f"({sch} fallback)")
        if os.path.exists(dst):
            os.remove(dst)                   # a partial must never reach Step 3
    return dst, False, time.time() - t0, "failed or truncated on both ftp:// and https://"

# --- probe: race the schemes on the SMALLEST file, keep the winner ---------------
sizes = {u: (_remote_size(_as(u, "https")) or 1 << 60) for u in files}
probe  = min(sizes, key=sizes.get)
timing = {}
for sch in ("https", "ftp"):
    t0 = time.time()
    timing[sch] = (time.time() - t0) if _curl(_as(probe, sch), "/tmp/_geo_probe.bin", 60) else float("inf")
os.path.exists("/tmp/_geo_probe.bin") and os.remove("/tmp/_geo_probe.bin")
scheme = min(timing, key=timing.get)
if timing[scheme] == float("inf"):
    raise RuntimeError("neither ftp:// nor https:// reached ftp.ncbi.nlm.nih.gov")
print("transport probe on %s (%.0f KB): %s -> using %s://"
      % (os.path.basename(probe), sizes[probe] / 1e3,
         ", ".join("%s=%.2fs" % (s, t) for s, t in timing.items()), scheme), flush=True)

# Print as files land — six silent curls behind a join look identical to a hang.
print(f"Fetching {len(files)} file(s) into {DEST} over {scheme}://", flush=True)
with ThreadPoolExecutor(max_workers=min(6, len(files))) as ex:
    futs = {ex.submit(_fetch, u, scheme): u for u in files}
    done, total = 0, 0
    for fut in as_completed(futs):
        dst, ok, dt, note = fut.result()
        done += 1
        sz = os.path.getsize(dst) if os.path.exists(dst) else 0
        total += sz
        print(f"  [{done}/{len(files)}] {'✓' if ok else '✗'} {os.path.basename(dst)}  "
              f"{sz/1e6:.1f} MB in {dt:.1f}s {note}", flush=True)
        if not ok:
            raise RuntimeError(f"download failed for {dst}: {note}")

print(f"Done — {total/1e6:.1f} MB, every file size- and gzip-verified.", flush=True)
```

Per-sample 10x triplets are 30–100 MB (`matrix.mtx.gz` dominates); the parallel fetch
above lands in tens of seconds on either transport.
**For larger archives** — `GSE…_RAW.tar` is often multi-GB — run as a background
job, write to a durable path, stream (`requests.get(stream=True)` +
`iter_content`) rather than buffering in the kernel, and `tar.extract()` only
the members whose names start with your `GSM…` prefix. That's the *one*
legitimate place to use `curl -C -` (resume on an interrupted multi-GB pull) — and
still verify the final size against `Content-Length` before extracting.

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
- **A download is not "done" because curl exited 0 — verify it.** Truncated
  `.gz` files are the #1 way this recipe wastes an hour: `curl` returns 0 (or is cut
  by `--max-time`), the partial file looks plausible, and the failure surfaces much
  later as a `gzip`/`sc.read_mtx` error that reads like corrupt *data* rather than an
  incomplete *download*. Always compare against the server's `Content-Length`/FTP
  `SIZE` and inflate-test `.gz` files (Step 2's `_complete()`), delete partials, and
  retry on the other scheme.
- **Neither `ftp://` nor `https://` is universally faster — probe.** Same 47 MB file
  on `ftp.ncbi.nlm.nih.gov`: from a VBC compute node https is ~5x faster
  (20.6 vs 3.8 MB/s, 2026-07); from a household line ftp was ~15x faster
  (~6 MB/s vs ~0.4 MB/s throttled, 2026-06). Both measurements are real. Step 2
  races the two on the smallest file and uses the winner, so it self-tunes per site
  instead of encoding one network's result as universal advice. Symptom of getting
  this wrong: downloads that crawl or stall for minutes.
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
- **Per-file timeout is a backstop, not a correctness check.** A short `--max-time`
  used to be the defence against a silent stall, but it *causes* the truncation above:
  curl is killed mid-stream and leaves a plausible-looking partial. Now that Step 2
  verifies completeness and deletes partials, prefer a generous timeout (300s) and let
  verification catch the bad transfer. Keep the per-file progress prints for the
  heartbeat.
- Some series only deposit RDS/Seurat objects → note it; FASTQ path may be cleaner.
- `GEOparse` caches SOFT files in `destdir`; reuse it to avoid refetching metadata.

## In ABA
`run_python`; `ensure_capability("GEOparse")`, `ensure_capability("scanpy")`
(pulls anndata), `ensure_capability("pandas")`. Big file pulls → background/streamed
job, durable disk path, verify size before load.
