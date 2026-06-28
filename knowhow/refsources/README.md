# Reference-source catalog (`refsources/`)

The install-modifiable data layer behind `fetch_reference` (misc/refs.md §5.1).
Each `*.yaml` here is a **provider manifest** telling the resolver *where* a
standard reference (genome FASTA, annotation, transcriptome) or a **pre-built
aligner index** lives, and how to turn a request
(`organism` / `assembly` / `role`, or an `accession`) into a concrete
download URL or CLI invocation. This is *data*, not code — bump a release or add
a species by editing YAML; no code change.

The schema is whatever `aba/backend/core/data/refsources.py` reads — that file is
the source of truth; this README documents it.

## How resolution works

`resolve_asset(provider, organism=, assembly=, role=, accession=, filename=)`:

1. `load_providers()` reads every `*.yaml` across the search path and keys them
   by their `provider:` field. **First file to declare a given provider name
   wins** (see *Layering* below).
2. It dispatches on the provider's `kind` (default `manifest`).

### Facets are matched NORMALIZED

The resolver folds case + separators so the agent's natural inputs hit:

- **role** — case/separator-folded (`Genome` == `genome`, `bowtie2-index` ==
  `bowtie2_index`). Match is exact-after-fold.
- **assembly** — case/separator-folded (`GRCh38` == `grch38`; `R64-1-1` ==
  `r64_1_1`). Exact-after-fold; omit it in a query to take the first asset of
  that role+organism.
- **organism** — folded to a canonical slug **with common-name aliases** AND a
  substring match. So `human`, `Homo sapiens`, and `homo_sapiens` all resolve to
  the `homo_sapiens` assets; `drosophila` matches `drosophila_melanogaster`.
  Aliases known to the resolver: human, mouse, fly/fruit fly/drosophila,
  zebrafish, rat, yeast/budding yeast, worm/c elegans. (No alias for
  arabidopsis — query it as `arabidopsis_thaliana` or `arabidopsis`, both of
  which substring-match.)

Use the slugs in `knowhow/refs/NAMING.md`: organisms as the lowercase binomial
with `_` (`homo_sapiens`, `mus_musculus`, …); roles from the controlled
vocabulary.

## The two `kind`s and every field

### `kind: manifest` — an explicit asset list

Used by FASTA/annotation/index providers (ensembl, gencode, ucsc, aws-indexes,
tenx-cellranger, igenomes). The resolver walks `assets:` and returns the first
whose `role` (folded) matches, `organism` passes the alias/substring test, and
`assembly` (folded) matches — skipping any facet the caller left `None`.

```yaml
provider: <name>          # the lookup key; must be unique across the search path
description: <one line>
kind: manifest
homepage: <url>           # informational
assets:
  - role: <role-slug>     # genome, transcriptome, cds, protein, gtf, gff,
                          # bowtie2_index, hisat2_index, star_index, twobit, ...
    organism: <slug>      # homo_sapiens, mus_musculus, ...
    assembly: <assembly>  # GRCh38, GRCm39, hg38, R64-1-1, ...
    version: <tag>        # provenance: ensembl-116, gencode-v50, GRCh38_noalt_as
    url: <direct download URL>
    unpack: zip | tar.gz | gz   # OMIT when the file is used as-is (.fa.gz, .2bit)
```

The resolver returns `{url, unpack, version, role, organism, assembly}`. The
fetch executor pulls `url`, applies `unpack`, and `register_reference`s it with
`version` recorded for reproducibility.

`unpack` semantics: set it when the download is an **archive of a tree** (a
zip/tarball that expands to multiple files or a dir — index packages, GENCODE,
10x). Omit it for a single already-final file even if itself compressed: a bare
`.fa.gz` / `.gtf.gz` is registered as-is (the consumer gunzips on use), and a
`.2bit` is binary. (This matches the seed: the `.fa.gz` FASTAs carry no
`unpack`; the `.zip`/`.tar.gz` index bundles do.)

### `kind: template` — parametric by accession

Used by accession-driven providers (ncbi). The resolver looks up `role` in
`roles{}`, builds params `{accession, filename, **roles[role]}`, and `.format`s
them into `command` and/or `url_template`. **The template branch receives only
`accession` + `filename` + the per-role params — NOT organism/assembly** — so it
is for accession-keyed sources only; species/assembly/release-driven sources use
`manifest`.

```yaml
provider: ncbi
kind: template
command: "datasets download genome accession {accession} --include {include} --filename {filename}"
unpack: zip
url_template: "<optional curl-able URL with {accession}/{include}/...>"   # optional
roles:
  genome: {include: genome}     # per-role params merged into the format() dict
  gtf:    {include: gtf}
  ...
```

The resolver returns `{command, url (if url_template), unpack, version=accession,
role, accession}`. `filename` defaults to `<accession>.zip` if the caller omits
it. An unknown `role` (not in `roles{}`) or a missing `accession` raises.

## Layering / search path (first match wins)

`refsources.py` searches, in order:

1. **`$ABA_REFSOURCES_DIR`** — operator / test override (this pack).
2. **Built-in seed** — `aba/backend/content/bio/knowhow/refsources/` (ships
   `aws-indexes` + `ncbi` minimal).

Because providers are keyed by name and the first dir wins, a provider file
**here overlays the seed of the same name**: e.g. this pack's `aws-indexes.yaml`
(comprehensive) fully replaces the 3-asset seed, and this `ncbi.yaml` replaces
the seed `ncbi`. New names (ensembl, gencode, ucsc, tenx-cellranger, igenomes)
are simply added.

### Activate this pack

Point the env var at this directory:

```bash
export ABA_REFSOURCES_DIR=/path/to/aba-recipe-pack/knowhow/refsources
```

(or wire it via the install's environment). Then `fetch_reference` /
`resolve_asset` see all seven providers below. Verify:

```bash
ABA_REFSOURCES_DIR=$ABA_REFSOURCES_DIR python -c \
 "import sys; sys.path.insert(0,'/path/to/aba/backend'); \
  from core.data.refsources import load_providers; print(sorted(load_providers()))"
# -> ['aws-indexes','ensembl','gencode','igenomes','ncbi','tenx-cellranger','ucsc']
```

### Extend / override a provider

- **Override an asset:** drop a YAML here with an existing `provider:` name; it
  shadows the seed entirely (it's whole-file, first-match-wins — not a merge).
  To keep the seed's other entries, copy them in.
- **Add a species/role:** append an `assets:` entry with the right
  role/organism/assembly slugs and a verified `url`.
- **Add a provider:** new file, new `provider:` name.

### Bump a release

Edit the YAML in place:

- **ensembl** — change `release-N` in the URL path and the `version` tag; for
  GTF/GFF also change the `.N.` filename infix. Re-confirm the assembly string
  (Ensembl re-bases assemblies between releases — rat went mRatBN7.2 → GRCr8,
  fly BDGP6.32 → BDGP6.54). Note vertebrates use the vertebrate release number
  (116) while fly/worm/yeast carry the Ensembl Genomes number (63) in their
  annotation filenames; arabidopsis is on the Plants mirror at EBI.
- **gencode** — change `release_NN`/`release_MNN` in the path AND the
  `vNN`/`vMNN` filename infix, plus `version`.
- **ucsc** — UCSC db names are stable; add a db by following the
  `goldenPath/{db}/bigZips/{db}.fa.gz` (+ `.2bit`) scheme.
- **aws-indexes / tenx-cellranger / igenomes** — swap the published archive
  name / dated reference filename and the `version` tag.

After any edit, re-HEAD a couple of changed URLs and re-run the resolve battery
(see the task's §5 / the in-repo `test_refs_fetch`).

## Providers in this pack

| File | kind | provider | what it serves |
|---|---|---|---|
| `ensembl.yaml` | manifest | `ensembl` | genome/cdna/cds/pep FASTA + GTF/GFF3, 8 species (vertebrate r116, metazoa/fungi r63, plants r63) |
| `gencode.yaml` | manifest | `gencode` | human v50 + mouse vM39 GTF/GFF3 + transcripts + primary-assembly genome |
| `ucsc.yaml` | manifest | `ucsc` | goldenPath `.fa.gz` (`genome`) + `.2bit` (`twobit`) for 9 UCSC dbs |
| `aws-indexes.yaml` | manifest | `aws-indexes` | pre-built Bowtie2 + HISAT2 indices (Langmead lab S3) |
| `tenx-cellranger.yaml` | manifest | `tenx-cellranger` | Cell Ranger GEX (`star_index`) / VDJ (`vdj_ref`) / ATAC (`atac_ref`) reference dirs |
| `ncbi.yaml` | template | `ncbi` | `datasets download genome accession …` by RefSeq/GenBank accession |
| `igenomes.yaml` | manifest | `igenomes` | iGenomes single-file `genome.fa` + `genes.gtf` (index dirs excluded — they are directories) |

### Role vocabulary used here

`genome`, `transcriptome`, `cds`, `protein`, `gtf`, `gff`, `twobit`,
`bowtie2_index`, `hisat2_index`, `star_index` (Cell Ranger GEX dir), plus the
deliberate extensions `vdj_ref` and `atac_ref` (10x). The resolver matches role
by normalized string equality, so a caller must request these exact slugs.
Index roles available elsewhere in the platform's vocabulary
(`bowtie_index`/bowtie1, `salmon_index`, `bwa_index`, `bismark_index`, …) are
**not** in this pack because no provider here publishes them as a single fetch —
acquire those by fetching a FASTA (ensembl/ucsc/gencode) and building.

### Organism slugs covered

`homo_sapiens`, `mus_musculus`, `rattus_norvegicus`, `danio_rerio`,
`drosophila_melanogaster`, `caenorhabditis_elegans`, `saccharomyces_cerevisiae`,
`arabidopsis_thaliana` (coverage varies by what each provider actually hosts —
human + mouse are covered most thoroughly).
