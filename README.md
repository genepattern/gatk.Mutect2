# Mutect2 (v4.3.0.0)

**Description**: Calls somatic short mutations — single nucleotide variants (SNVs) and small insertions/deletions (indels) — from tumor sequencing data, with or without a matched normal sample, via local assembly of haplotypes using the GATK Mutect2 Bayesian somatic genotyping engine.
**Authors**: Broad Institute; GenePattern Team, UC San Diego
**Contact**: https://groups.google.com/forum/#!forum/genepattern-help
**Algorithm Version**: GATK 4.3.0.0

---

## Summary

**gatk.Mutect2** is a GenePattern module wrapping the GATK (Genome Analysis Toolkit) Mutect2 somatic variant caller. It detects somatic mutations — genetic changes present in tumor tissue but absent from the germline — by performing local *de novo* haplotype assembly over active genomic regions and scoring candidate variants with a Bayesian somatic genotyping model.

### What problem does it solve?

Identifying somatic mutations (mutations acquired during tumor development, as opposed to inherited germline variants) is a central task in cancer genomics. Standard germline variant callers are not appropriate for this task because somatic mutations are often present at low allele fractions, may be heterogeneous across a tumor, and arise against a background of technical noise and germline variation. Mutect2 is specifically designed to overcome these challenges.

### How does it work?

1. **Active region detection**: Mutect2 scans the input BAM file(s) for genomic regions showing evidence of variation. Only these "active regions" are processed in depth.
2. **Local assembly**: Within each active region, the tool assembles the reads into a set of candidate haplotypes using a De Bruijn graph.
3. **Likelihood scoring**: Each read is aligned to each candidate haplotype using a pair hidden Markov model (PairHMM), and the genotype likelihoods are computed.
4. **Somatic genotyping**: A Bayesian model evaluates the evidence for each variant being somatic rather than germline, using optional population allele frequency priors (from a germline resource) and matched normal data.
5. **Filtering** (optional): After calling, `FilterMutectCalls` can be run to annotate variants as `PASS` or assign one or more failure reasons, leveraging read orientation artifact models (e.g., for FFPE samples), contamination estimates, and statistical filters.

### When to use this module

- Whole-genome sequencing (WGS) or whole-exome sequencing (WES) of tumor-normal pairs
- Tumor-only somatic calling (no matched normal available)
- Mitochondrial variant calling / heteroplasmy detection
- Construction of a Panel of Normals (PoN) — run Mutect2 in single-sample mode on a set of normal samples, then use `GenomicsDBImport` and `CreateSomaticPanelOfNormals`

### Companion index file handling

The module wrapper automatically checks for required companion index files (`.bai` for BAM files; `.tbi` for VCF.gz files; `.fai` and `.dict` for the reference FASTA). If any index is missing, the wrapper generates it on-the-fly using the appropriate GATK/samtools utility (`gatk BuildBamIndex`, `gatk IndexFeatureFile`, `samtools faidx`, `gatk CreateSequenceDictionary`) before running the main analysis.

---

## References

1. Benjamin, D. et al. (2019). **Calling Somatic SNVs and Indels with Mutect2**. *bioRxiv*. https://doi.org/10.1101/861054
2. McKenna, A. et al. (2010). **The Genome Analysis Toolkit: A MapReduce framework for analyzing next-generation DNA sequencing data**. *Genome Research*, 20(9):1297–1303. https://doi.org/10.1101/gr.107524.110
3. Van der Auwera, G.A. & O'Connor, B.D. (2020). *Genomics in the Cloud: Using Docker, GATK, and WDL in Terra*. O'Reilly Media.
4. GATK Mutect2 Documentation: https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2
5. GATK Best Practices for Somatic SNVs/Indels: https://gatk.broadinstitute.org/hc/en-us/articles/360035894731

---

## Source Links

* [GATK Source Repository (GitHub)](https://github.com/broadinstitute/gatk)
* [GenePattern Mutect2 Docker Image](https://hub.docker.com/r/genepattern/mutect2)
* [GenePattern Mutect2 Dockerfile](https://github.com/genepattern/Mutect2/blob/main/Dockerfile)

---

## Parameters

| Name | Description | Default Value |
| :--- | :--- | :--- |
| tumor.bam * | Tumor sample BAM (or CRAM) file with index. The wrapper auto-generates a `.bai` index if absent. | — |
| reference.fasta * | Reference genome FASTA file. Companion `.fai` and `.dict` files are auto-generated if absent. | — |
| normal.bam | Matched normal sample BAM (or CRAM) file. If provided, enables tumor-normal calling mode. The wrapper auto-generates a `.bai` index if absent. | — |
| output.vcf | Base name for output VCF file(s). Final output will be `<output.vcf>_unfiltered.vcf.gz` (and `_filtered.vcf.gz` if filtering is enabled). | output.vcf.gz |
| germline.resource | Population germline allele frequency VCF.gz (e.g., `af-only-gnomad.hg38.vcf.gz`). The wrapper auto-generates a `.tbi` index if absent. | — |
| panel.of.normals | Panel of Normals (PoN) VCF.gz for artifact filtering (e.g., `1000g_pon.hg38.vcf.gz`). The wrapper auto-generates a `.tbi` index if absent. | — |
| intervals | Genomic intervals to restrict variant calling (`.intervals`, `.list`, `.bed`, or `.interval_list`). Strongly recommended for WES/targeted panels. | — |
| tumor.sample.name | SM tag of the tumor sample as recorded in the BAM `@RG` header. Auto-detected from BAM header if left blank. | — |
| normal.sample.name | SM tag of the normal sample as recorded in the BAM `@RG` header. Auto-detected from BAM header if left blank. Only used in tumor-normal mode. | — |
| af.of.alleles.not.in.resource | Prior allele fraction assigned to variants not found in the germline resource. Mode-dependent default: `1e-6` (tumor-normal), `5e-8` (tumor-only), `4e-3` (mitochondria). | 0.0000025 |
| mitochondria.mode | Enable mitochondrial variant calling mode with specialized LOD thresholds suited for heteroplasmy. | false |
| max.reads.per.alignment.start | Downsampling cap: maximum reads retained per alignment start position. Increase to ≥200 for amplicon or very high-depth data. | 50 |
| callable.depth | Minimum read depth required at a site to be considered callable during filtering. | 10 |
| extra.args | Additional GATK Mutect2 command-line arguments passed directly to the tool (e.g., `--genotype-germline-sites true`). | — |

\* required

---

## Input Files

1. **tumor.bam**
   The primary input: a coordinate-sorted BAM (Binary Alignment/Map) or CRAM file containing the sequencing reads from the **tumor sample**. The file must include at least one read group (`@RG`) header line with a sample name (`SM:`) tag. A companion `.bai` BAM index file is required for random access; if it does not exist alongside the BAM file, the wrapper will automatically generate it by running `gatk BuildBamIndex`. Accepted formats: `.bam`, `.cram`.

2. **reference.fasta**
   The reference genome in FASTA format (`.fa` or `.fasta`) against which the reads were aligned. Must match the reference used during alignment. Two companion files are required:
   - `.fai` — FASTA index (generated by `samtools faidx` if absent)
   - `.dict` — Sequence dictionary (generated by `gatk CreateSequenceDictionary` if absent)
   
   The wrapper will auto-generate either file if it is missing.

3. **normal.bam** *(optional)*
   A coordinate-sorted BAM or CRAM file for the **matched normal sample** from the same patient. Providing a matched normal dramatically improves somatic specificity by allowing Mutect2 to distinguish somatic mutations from germline variants and sequencing artifacts. As with `tumor.bam`, a `.bai` index is auto-generated if absent. Accepted formats: `.bam`, `.cram`.

4. **germline.resource** *(optional)*
   A block-compressed, tabix-indexed VCF (`.vcf.gz`) containing population germline allele frequencies — typically the gnomAD allele-frequency-only VCF (e.g., `af-only-gnomad.hg38.vcf.gz` from the GATK resource bundle). Must include an `AF` INFO field. Used to set germline prior probabilities during somatic genotyping. A `.tbi` index is auto-generated via `gatk IndexFeatureFile` if absent.

5. **panel.of.normals** *(optional)*
   A block-compressed, tabix-indexed VCF (`.vcf.gz`) representing a Panel of Normals (PoN): a collection of variant sites observed in multiple normal samples, representing recurrent sequencing and technical artifacts. Any variant found in the PoN is flagged during filtering. The GATK resource bundle provides `1000g_pon.hg38.vcf.gz` as a community PoN. A `.tbi` index is auto-generated via `gatk IndexFeatureFile` if absent.

6. **intervals** *(optional)*
   A file specifying genomic regions to restrict variant calling. **Strongly recommended** for whole-exome sequencing (WES) and targeted panel sequencing to reduce runtime and limit off-target calls. Accepted formats: `.intervals`, `.list` (one `chr:start-end` entry per line), `.bed`, or `.interval_list` (Picard-format). For WGS without specific targets, this parameter may be omitted.

---

## Output Files

1. **`<output.vcf>_unfiltered.vcf.gz`**
   A block-compressed, tabix-indexed VCF containing all candidate somatic variants called by Mutect2, including both true somatic mutations and likely artifacts. Each record includes extensive annotations: allele depths (`AD`), allele fractions (`AF`), genotype likelihoods, and supporting read counts. This file is the direct output of the Mutect2 calling step before any post-call filtering.

2. **`<output.vcf>_unfiltered.vcf.gz.tbi`**
   Tabix index for the unfiltered VCF, enabling fast random access by genomic coordinate.

3. **`<output.vcf>_filtered.vcf.gz`** *(produced when filtering is enabled)*
   The post-filtered VCF produced by `FilterMutectCalls`. Each variant is annotated with either `PASS` (considered a high-confidence somatic call) or one or more filter flags (e.g., `germline`, `panel_of_normals`, `strand_bias`, `weak_evidence`, `orientation_bias`). Downstream analysis should typically use this file.

4. **`<output.vcf>_filtered.vcf.gz.tbi`** *(produced when filtering is enabled)*
   Tabix index for the filtered VCF.

5. **`<output.vcf>_f1r2.tar.gz`** *(produced when orientation bias filtering is enabled)*
   Compressed archive of F1R2 read orientation counts per site, used as input to `LearnReadOrientationModel`. Present only when the orientation bias filter option is enabled.

6. **`<output.vcf>_read-orientation-model.tar.gz`** *(produced when orientation bias filtering is enabled)*
   The fitted read orientation artifact model (output of `LearnReadOrientationModel`), passed to `FilterMutectCalls`. Present only when the orientation bias filter option is enabled.

7. **`mutect2_run.log`**
   A plain-text log file capturing the standard output and error streams of all GATK commands executed, including any auto-generated index steps. Useful for diagnosing errors or auditing the exact commands run.

---

## Example Data

**Input:**
- Tumor BAM: HG008-T (NIST/GIAB HG008 tumor sample)
- Normal BAM: HG008-N-D (NIST/GIAB HG008 matched normal, DNA)
- Reference: hg38 (GRCh38) — available from the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)
- Germline resource: `af-only-gnomad.hg38.vcf.gz` — available from the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)
- Panel of Normals: `1000g_pon.hg38.vcf.gz` — available from the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811)

**Output:**
- Example filtered VCF: `HG008-T_filtered.vcf.gz` — somatic SNVs and indels called against HG008-N-D, filtered with `FilterMutectCalls`

> **Note**: GIAB HG008 tumor-normal data and truth sets are available from the [NIST Genome in a Bottle Consortium](https://www.nist.gov/programs-projects/genome-bottle).

---

## Requirements

- **Platform**: GenePattern server (version 3.9.11 or later recommended)
- **Language/Runtime**: Java 8+ (provided within the Docker container)
- **Docker Image**: `genepattern/mutect2:latest` (based on `broadinstitute/gatk:4.3.0.0`)
- **Memory**: Minimum 8 GB RAM recommended; 16–32 GB for WGS or high-depth samples. Adjust `java.heap.size` accordingly.
- **CPU**: Multi-core recommended; the PairHMM step scales with `--native-pair-hmm-threads`.
- **Disk**: Sufficient scratch space for intermediate BAM processing and VCF output. WGS runs may require 50–200 GB of free disk.
- **Input BAM requirements**:
  - Must be coordinate-sorted
  - Must contain `@RG` read group headers with `SM:` sample name tags
  - Duplicate marking is recommended (e.g., via `MarkDuplicates`) prior to running Mutect2
  - Base Quality Score Recalibration (BQSR) is recommended for best results
- **Reference requirements**:
  - Must match the genome build used for alignment (hg38 recommended)
  - `.fai` and `.dict` companion files required (auto-generated if absent)

---

## License

This GenePattern module wraps the open-source GATK toolkit.
- **GATK License**: BSD 3-Clause — see https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT
- **GenePattern Module**: MIT License — see https://github.com/genepattern/Mutect2/blob/main/LICENSE

---

## Version Comments

| Version | Release Date | Description |
| :--- | :--- | :--- |
| 4.3.0.0 | 2023-06-01 | Initial GenePattern module release wrapping GATK Mutect2 v4.3.0.0; includes auto-indexing of BAM, VCF, and FASTA companion files; tumor-normal and tumor-only modes; optional FilterMutectCalls and LearnReadOrientationModel steps. |
