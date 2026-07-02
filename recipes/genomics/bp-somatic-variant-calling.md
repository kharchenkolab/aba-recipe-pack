---
name: bp-somatic-variant-calling
description: Best-practice somatic (tumor) variant calling from DNA sequencing — tumor-normal Mutect2+Strelka
  for SNV/indel, ASCAT/Control-FREEC for CNV, Manta for SV, via nf-core/sarek.
when_to_use: 'Use to call SOMATIC (acquired, cancer/tumor) mutations from tumor DNA sequencing (WGS, WES,
  or targeted panel): SNV/indel, and optionally CNV, SV, purity/ploidy. Covers tumor-normal PAIRED and
  TUMOR-ONLY (no matched normal, use panel-of-normals). NOT for germline/inherited variant calling, NOT
  for RNA variants, NOT for microbial/pathogen variants.'
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- sarek
- gatk_mutect2
- strelka
- ascat
keywords:
- somatic variant calling
- tumor-normal
- tumor-only
- Mutect2
- Strelka2
- panel of normals
- gnomAD germline resource
- ASCAT purity ploidy
- Control-FREEC CNV
- Manta structural variant
- cancer WGS WES panel
- somatic SNV indel
produces:
- somatic_snv_indel.vcf.gz
- cnv_calls.txt
- sv.vcf.gz
- purity_ploidy_estimates
- annotated_variants.vcf.gz
- multiqc_report.html
domain: genomics
source: GATK Best Practices — Somatic short variant discovery (SNVs + Indels), Mutect2 with gnomAD germline
  resource + panel-of-normals (gatk.broadinstitute.org); nf-core/sarek (Hanssen et al., NAR Genom Bioinform
  2024, lqae031)
---

**Question:** How do I call somatic (cancer/tumor) mutations from DNA sequencing — which mode (tumor-normal vs tumor-only), which callers, and how do I get CNV/SV and purity/ploidy?

## Recommended approach (with why)

**Sequence a matched normal (tumor-normal PAIRED design) whenever possible, and call SNV/indel with Mutect2 as the primary caller.** This is the GATK Best-Practices somatic workflow and the field default. A matched normal (e.g. adjacent normal tissue or blood) lets the caller subtract that patient's germline variants directly, which is by far the strongest control: when lung tumors were sequenced tumor-only, ~94% of "somatic" calls were false positives, and even filtering common SNPs left ~48% false positives ([tumor-only filtration study, JMD 2019](https://pubmed.ncbi.nlm.nih.gov/30576869/)). Mutect2 does local haplotype reassembly and Bayesian somatic modeling, and in recent WES benchmarks achieved the highest recall among common callers while retaining ~99.9% precision ([Mutect2/Strelka2/FreeBayes benchmark, Biomolecules 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12650410/)).

Even with a matched normal, Mutect2 should use **both** a germline-resource (gnomAD allele frequencies) **and** a panel-of-normals (PoN). GATK is explicit that a germline resource such as gnomAD is a more refined tool for probabilistic germline filtering than any PoN, while the PoN captures recurrent technical/sequencing artifacts that germline databases miss ([GATK Mutect2 resources guide](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Mutect2)). Both are standard inputs to the best-practice workflow.

**Run a second SNV/indel caller and take the consensus for anything clinically or biologically load-bearing.** Callers disagree at low VAF: Mutect2 is strongest above ~10% VAF, while Strelka2's position-wise model detects lower-VAF variants (down to ~5%) ([benchmark, Biomolecules 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12650410/)). Ensemble/voting approaches (e.g. Mutect2 + Strelka2 + MuSE for SNVs; Mutect2 + Strelka2 + VarScan2 for indels) materially improve F1 over any single caller ([ensemble benchmark, Brief Bioinform 2025](https://academic.oup.com/bib/article/26/1/bbae697/7960049)). In practice: **Mutect2 + Strelka2** is the pragmatic two-caller default.

**Somatic calling is not just SNV/indel — add the variant classes your biology needs:**
- **CNV / copy number:** call **ASCAT** (allele-specific, and it estimates tumor **purity and ploidy** — essential for interpreting VAFs and copy states) for WGS/WES with a matched normal; **Control-FREEC** or **CNVkit** are the alternatives (CNVkit is the go-to for targeted/panel and for tumor-only CNV).
- **Structural variants (SV):** call **Manta** (fusions, large deletions, translocations) on WGS; **TIDDIT** as an alternative.

## Alternatives (competing approaches, and when each wins)

- **Strelka2 as primary SNV/indel caller** — when you care most about low-VAF/subclonal SNVs or want strict, high-precision indel calls; often paired with Manta for indel candidate seeding. Best combined with, not instead of, Mutect2.
- **Tumor-only mode (no matched normal)** — use only when a matched normal is genuinely unavailable (archival FFPE, cell lines, retrospective cohorts). Requires an **aggressive PoN plus population-frequency germline filtering**: filter on gnomAD/ExAC AF > 0.01 and rescue known somatic sites via COSMIC; ~10 unmatched normals in the PoN maintains ~94% sensitivity / ~99% specificity ([UNMASC, NAR Cancer 2021](https://academic.oup.com/narcancer/article/3/4/zcab040/6382329); [tumor-only filtration, JMD 2019](https://pubmed.ncbi.nlm.nih.gov/30576869/)). Expect residual germline contamination; do not report tumor-only germline-looking variants as somatic without orthogonal evidence.
- **FreeBayes / Lofreq / MuSE / VarScan2 / Sentieon TNScope** — supplementary callers for ensemble voting or specific niches (Lofreq/MuSE for very low VAF; TNScope for speed at scale). Not first-line alone.
- **Consensus/ensemble callers** (e.g. SomaticSeq-style voting) — when maximal accuracy justifies the extra compute and you can run ≥3 callers; overkill for exploratory work.

## Decision features (what your data dictates)

- **Matched normal available?** → the single biggest fork. Yes → tumor-normal paired (default). No → tumor-only + PoN + population filtering (weaker; flag uncertainty).
- **Assay: WGS vs WES vs targeted panel.** WGS enables genome-wide CNV/SV and ASCAT purity/ploidy. WES/panel need `--wes`/`--intervals` (target BED); on panels prefer CNVkit for CNV and expect limited SV power.
- **Sequencing depth / expected VAF.** Deep panels (500–1000×+) chase low-VAF/subclonal or ctDNA variants → favor low-VAF-sensitive callers (Strelka2, Lofreq, MuSE) and higher-quality PoN. Standard WGS ~30–60× tumor / ~30× normal is fine for clonal SNV/indel.
- **Purity (tumor content) and ploidy.** Low-purity samples suppress somatic VAF and confound CNV; ASCAT estimates both, and its purity/ploidy can be overridden if you have orthogonal estimates.
- **Sample type / FFPE.** FFPE deamination artifacts (C>T/G>A) inflate false positives → orientation-bias filtering (Mutect2 handles this) and a matched FFPE-derived PoN matter.
- **Variant classes needed.** Driver SNV/indel only → skip CNV/SV. Copy-number-driven cancers, LOH, WGD, or fusion questions → add ASCAT/Control-FREEC and Manta.
- **Cohort with multiple samples per patient** (primary + relapse/mets) — sarek supports multiple tumor samples against one normal.

## Pitfalls

- **Tumor-only without a PoN is not somatic calling** — it is mostly germline + artifact. If no normal, a PoN and gnomAD/COSMIC filtering are mandatory, and results remain provisional.
- **Skipping the germline resource** — Mutect2 without `--germline-resource` (gnomAD) leaves many germline variants uncalled-as-germline; the PoN alone does not fix this.
- **Interpreting VAF without purity/ploidy** — a "subclonal" VAF in a low-purity or copy-altered region can be a clonal variant; get ASCAT purity/ploidy before calling clonality.
- **Single-caller reporting at low VAF** — callers diverge most exactly where it matters (subclonal/ctDNA); require concordance for low-VAF calls.
- **Forgetting `--wes`/`--intervals` on exome/panel** — running WGS defaults on targeted data wrecks CNV and coverage stats.
- **Contaminated or swapped normal** — always run sample-concordance/contamination QC; a mismatched normal silently produces garbage somatic calls.
- **FFPE artifacts** — treat C>T/G>A low-VAF calls skeptically; enable orientation-bias filtering.

## In ABA (which pipeline to run)

Route somatic variant calling to **nf-core/sarek** via `run_nextflow`. Sarek is the nf-core cancer variant pipeline: pre-processing → somatic (and germline) variant calling → annotation, for WGS / WES / targeted, in tumor-normal paired and tumor-only modes ([nf-core/sarek](https://nf-co.re/sarek); [Hanssen et al. 2024](https://academic.oup.com/nargab/article/6/2/lqae031/7658070)).

**Launch:** `run_nextflow` with pipeline `nf-core/sarek`, a fixed `revision` (pin the current stable release, e.g. `3.8.1` — do not run `dev`/`master`), and `profile` matching the environment (`singularity`/`apptainer` on this cluster; `docker` elsewhere). Provide the sarek samplesheet (`--input`) and `--outdir`.

**Selection routing (set via params, details in `describe_pipeline`):**
- **Tumor-normal:** pair tumor and normal rows under the same `patient` in the samplesheet (Sarek infers pairing from the `status` column). **Tumor-only:** provide tumor rows only and supply a `--pon` (panel of normals).
- **Callers via `--tools`:** default somatic SNV/indel = `mutect2,strelka`; add CNV with `ascat` (WGS, needs matched normal; yields purity/ploidy) or `controlfreec`/`cnvkit` (cnvkit for panels/tumor-only); add SV with `manta` (or `tiddit`); add annotation with `snpeff`/`vep`.
- **Assay:** for WES/panel set `--wes` and pass the target BED via `--intervals`.
- **Germline filtering resources:** point Mutect2 at the gnomAD germline resource and `--pon`; use iGenomes/GATK bundle for the reference.

Call **`describe_pipeline`** for the exact sarek parameter names, samplesheet schema (`patient,sample,status,lane,fastq_1,fastq_2` etc.), `--step` values (start from FASTQ, from mapped BAM, or `variant_calling` from recalibrated BAM), and the full `--tools` vocabulary before launching. Keep per-param detail there — this knowhow only decides mode, callers, and variant classes.
