---
name: bp-germline-variant-calling
description: Best-practice germline SNV/indel discovery from WGS/WES/panel reads — GATK-style align+BQSR
  then HaplotypeCaller (GVCF/joint) or DeepVariant, via nf-core/sarek.
when_to_use: Use to call germline (inherited) SNVs and short indels from short-read DNA sequencing — human/model-organism
  WGS, whole-exome (WES), or targeted capture panels — starting from FASTQ/BAM/CRAM, including single-sample
  calling and multi-sample joint genotyping of a cohort/trio. NOT for somatic/tumor variant calling (tumor-normal),
  NOT for structural/copy-number variants, NOT for RNA-seq variant calling.
requires_tools:
- run_nextflow
capabilities_needed:
- nextflow
- nf-core/sarek
- gatk4
- deepvariant
- bwa-mem2
keywords:
- germline variant calling
- SNV indel
- HaplotypeCaller
- DeepVariant
- joint genotyping
- GVCF
- BQSR
- VQSR
- WGS WES panel
- GATK best practices
- sarek
produces:
- annotated.vcf.gz
- joint_germline.vcf.gz
- gvcf
- multiqc_report.html
domain: genomics
source: GATK Best Practices — Germline short variant discovery (SNPs+Indels), gatk.broadinstitute.org;
  DeepVariant/GLnexus cohort calling (Yun et al., Bioinformatics 2020)
---

## Question
How do I call germline (inherited) SNVs and short indels from short-read WGS/WES/panel data, and which caller/workflow should I use?

## Recommended approach (with why)
Follow the **GATK Best-Practices germline short-variant discovery** shape, and run it through **nf-core/sarek**:

1. **Align** (BWA-MEM2) → **mark duplicates** → **BQSR** (BaseRecalibrator + ApplyBQSR).
2. **Call per sample** with either **GATK HaplotypeCaller** (local de-novo haplotype assembly, calls SNVs+indels jointly) or **DeepVariant** (CNN over pileup images).
3. **For a cohort/trio/family**, call each sample in **GVCF mode**, then **joint-genotype** all GVCFs together (GATK GenomicsDBImport → GenotypeGVCFs), then filter the multi-sample callset with **VQSR** (or hard filters for small cohorts).

**Why this shape.** GATK Best Practices is the reference standard for germline discovery in humans and is what benchmark studies, clinical pipelines, and reference callsets are built against ([GATK docs](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels)). **Joint genotyping** is preferred over merging independently-called single-sample VCFs because calling all samples together lets the genotyper distinguish a true homozygous-reference site from a no-data site, rescues low-depth evidence that is corroborated across samples, and produces a coherent squared-off matrix suitable for VQSR and Mendelian/segregation analysis ([GATK: logic of joint calling](https://gatk.broadinstitute.org/hc/en-us/articles/360035890431-The-logic-of-joint-calling-for-germline-short-variants)). The GVCF intermediate decouples per-sample calling from genotyping so a cohort can grow without re-calling every sample ("N+1" problem).

**Caller choice.** Multiple independent benchmarks report **DeepVariant matches or exceeds HaplotypeCaller in precision/recall**, with the largest gains on **indels and difficult regions** (homopolymers, low-complexity) and lower trio Mendelian-error rates ([Comparison of GATK and DeepVariant by trio sequencing, Sci Rep 2022](https://www.nature.com/articles/s41598-022-05833-4); [Accuracy of germline pipelines, Sci Rep 2020](https://www.nature.com/articles/s41598-020-77218-4)). DeepVariant needs no BQSR/VQSR tuning and is a strong default for single samples. HaplotypeCaller's edge is the mature, GenomicsDB-backed **joint-genotyping** path — the standard for cohorts, families, and any downstream that expects a GATK VQSR callset. (DeepVariant scales to cohorts via GLnexus merging of GVCFs — [Yun et al. 2020](https://academic.oup.com/bioinformatics/article/36/24/5582/6064144) — but that path is less turnkey in sarek.)

## Alternatives (with when each is preferred)
- **DeepVariant (single-sample or per-sample GVCF).** Prefer when you want best raw accuracy on indels/hard regions with minimal filtering effort, for one or few samples, or non-model organisms without good VQSR truth resources. Choose `--tools deepvariant`.
- **HaplotypeCaller + joint genotyping + VQSR.** Prefer for **cohorts, trios, and family/population studies** where you need a jointly-genotyped multi-sample VCF and standard VQSR filtering. Choose `--tools haplotypecaller --joint_germline`.
- **HaplotypeCaller single-sample + hard filters / CNN.** Prefer for a single WGS/WES sample when you want the GATK ecosystem but have no cohort; VQSR is unreliable with few samples, so use hard filters.
- **Strelka2 / FreeBayes / bcftools mpileup.** Lighter/alternative callers; Strelka2 is fast and accurate, bcftools is common for **non-human/agricultural** germline work where GATK resource bundles don't exist. Available in sarek as additional `--tools`.
- **Sentieon Haplotyper.** Drop-in, much faster HaplotypeCaller-equivalent (licensed); prefer when throughput matters and a Sentieon license is available.
- **Ensemble / multi-caller.** Run 2+ callers and intersect when maximizing precision on a fixed cohort matters more than turnaround.

## Decision features (what drives the choice)
- **Number of samples / study design.** Single sample → DeepVariant or single-sample HaplotypeCaller. Cohort / trio / family / population → GVCF + `--joint_germline` (HaplotypeCaller).
- **Assay variant.** WGS → default intervals. **WES / targeted panel → set `--wes` and provide the capture `--intervals` BED**; skipping this inflates runtime and calls off-target noise.
- **Organism.** Human/well-resourced model organism → VQSR feasible (needs known-sites/truth resources). Non-model/agricultural or no truth set → prefer DeepVariant or hard filters, and note BQSR/VQSR known-sites may be unavailable.
- **Region difficulty / indel importance.** Indel-heavy or difficult-region questions favor DeepVariant.
- **Depth & read length.** Standard short-read Illumina (30x WGS / 100x+ WES) fits all callers; very low depth benefits most from joint calling's cross-sample corroboration.
- **Downstream expectation.** Need a VQSR-filtered joint callset / trio segregation → HaplotypeCaller joint path. Need a quick, accurate per-sample VCF → DeepVariant.
- **Throughput / licensing.** Large batches with Sentieon license → Sentieon Haplotyper.

## Pitfalls
- **Merging single-sample VCFs is not joint genotyping** — you lose hom-ref vs no-call information and get an inconsistent matrix. Use GVCF → joint genotyping for cohorts.
- **VQSR with too few samples/variants fails or misbehaves** — for a single exome or tiny cohort use GATK hard filters instead.
- **Forgetting `--wes`/`--intervals` for capture data** wastes compute and adds off-target false positives; conversely, applying WGS assumptions to panels distorts depth-based filters.
- **Reference/build mismatch** — mapping reference, intervals BED, and known-sites/VQSR resources must all be the same genome build (and contig naming) or the run breaks or silently miscalls.
- **Skipping BQSR** is defensible for DeepVariant (it does not require it) but expected for the classic GATK path.
- **Wrong problem class:** this is germline. For tumor variant calling use the somatic tumor-normal workflow (Mutect2), not this knowhow; for SV/CNV use dedicated callers.
- **DeepVariant model must match the data type** (WGS vs WES vs PacBio/ONT); using the wrong model degrades accuracy.

## In ABA (which pipeline to run)
Route to **nf-core/sarek** (whole-genome/targeted germline + somatic variant analysis; current release **3.9.0**).

- Launch with `run_nextflow` targeting `nf-core/sarek`; pin a stable revision (`-r 3.9.0`) for reproducibility.
- **Profile:** use the site container profile (`singularity`/`apptainer` on this cluster) plus any institutional config; do not run `docker` on HPC.
- **Caller/step selection is driven by `--tools`:** `deepvariant` (single-sample default), `haplotypecaller` (add `--joint_germline` for cohorts), optionally `strelka`, `freebayes`, `mpileup`, or `sentieon_haplotyper`. Add annotation with `snpeff`/`vep`.
- **Assay handling:** for WES/panel pass `--wes` and `--intervals <capture.bed>`; for WGS leave intervals at genome default.
- Call **`describe_pipeline` for nf-core/sarek** to get the exact parameter names, defaults, reference/`igenomes` options, BQSR/known-sites resources, and the VQSR settings — keep that param detail there, not here.
- Decision summary for the agent: **single sample → `--tools deepvariant`; cohort/trio/family → `--tools haplotypecaller --joint_germline` (+VQSR); WES/panel → always add `--wes --intervals`; non-human/no truth set → DeepVariant or hard filters.**

