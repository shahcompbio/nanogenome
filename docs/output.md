# shahcompbio/nanogenome: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Phasing](#phasing) - Variant calling, phasing, and BAM haplotagging
- [Somatic SV calling](#somatic-sv-calling) - Structural variant detection from tumor-normal pairs
- [Germline SV calling](#germline-sv-calling) - Structural variant detection from normal samples
- [Consensus SV calling](#consensus-sv-calling) - Merging results from multiple callers with MINDA
- [Copy number analysis](#copy-number-analysis) - Haplotype-resolved CNA detection
- [TE calling](#te-calling) - Transposable element insertion detection
- [SV annotation](#sv-annotation) - Gene and clinical annotation of SVs
- [Visualization](#visualization) - Circos plots and karyoplots
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Phasing

<details markdown="1">
<summary>Output files</summary>

- `longphase/`
  - `*.vcf.gz`: Phased VCF files from LongPhase.
  - `*.vcf.gz.tbi`: Tabix index files.
- `whatshap/`
  - `*.tsv`: WhatsHap phasing statistics per sample.

</details>

The phasing workflow calls SNPs/indels using [Clair3](https://github.com/HKU-BAL/Clair3), phases variants using [LongPhase](https://github.com/twolinin/LongPhase), and haplotags BAM files using [WhatsHap](https://whatshap.readthedocs.io/). Phasing statistics are collected for MultiQC reporting.

### Somatic SV calling

<details markdown="1">
<summary>Output files</summary>

- `severus/`
  - `*.vcf`: Severus SV calls.
- `savana/`
  - `*.vcf.gz`: SAVANA classified SV calls.
- `nanomonsv/`
  - `*.vcf`: NanoMonSV assembly-based SV calls.

</details>

Somatic SV calling uses an ensemble approach with up to three callers: [Severus](https://github.com/KolmogorovLab/Severus), [SAVANA](https://github.com/cortes-ciriano-lab/savana), and [NanoMonSV](https://github.com/friend1ws/nanomonsv). Callers can be selected with `--somatic_callers` (default: `severus,savana,nanomonsv`).

### Germline SV calling

<details markdown="1">
<summary>Output files</summary>

- `severus/`
  - `*.vcf`: Severus germline SV calls.
- `longcalld/`
  - `*.vcf`: LongcallD germline SV calls.

</details>

Germline SV calling (enabled with `--germline`) uses [Severus](https://github.com/KolmogorovLab/Severus), [Sniffles](https://github.com/fritzsedlazeck/Sniffles), [CuteSV](https://github.com/tjiangHIT/cuteSV), and [LongcallD](https://github.com/ydLiu-HIT/LongcallD). Callers can be selected with `--germline_callers` (default: `severus,longcallD,cutesv,sniffles`).

### Consensus SV calling

<details markdown="1">
<summary>Output files</summary>

- `minda/`
  - `*_minda_union.vcf`: Union of all SV calls from every caller (always produced).
  - `*_minda_consensus.vcf`: Consensus SV calls supported by at least `--min_callers` callers (default: 2).
  - `*_min_callers_<N>/`: Full MINDA output directory containing per-caller breakdowns and ensemble results.

</details>

[MINDA](https://github.com/shahcompbio/minda) merges SV calls from multiple callers into both union and consensus call sets. MINDA is always run twice — once with `min_callers=1` to produce the union set (all calls from any caller), and once with the user-specified `--min_callers` (default: 2) to produce the consensus set. The `--tolerance` parameter (default: 100bp) controls the breakpoint distance for merging, and `--min_size` (default: 50bp) filters out small variants.

### Copy number analysis

<details markdown="1">
<summary>Output files</summary>

- `wakhan/`
  - Copy number segments and haplotype-resolved CNA profiles from Wakhan.
- `savana/`
  - SAVANA CNA output when enabled via `--cna_tools`.

</details>

Haplotype-resolved copy number analysis is performed by [Wakhan](https://github.com/shahcompbio/wakhan) and optionally [SAVANA](https://github.com/cortes-ciriano-lab/savana). This step integrates SV haplotype information for improved phasing of copy number segments when `--use_sv_haplotypes` is enabled.

### TE calling

<details markdown="1">
<summary>Output files</summary>

- `longcalld/`
  - `*.all_variants.vcf`: All variant calls from LongcallD TE mode including TE-mediated insertions.
  - `*.sv_only.vcf.gz`: Filtered structural variants for annotation.
  - `*_phased_refined.sorted.cram`: Refined read alignments (when `--longcalld_realign` is enabled).
  - `*_phased_refined.sorted.cram.crai`: CRAM index file.
- `annotsv/`
  - AnnotSV annotation of LongcallD TE calls.

</details>

Transposable element calling (enabled with `--te_calling`) detects TE-mediated insertions using [LongcallD](https://github.com/ydLiu-HIT/LongcallD) and/or [tldr](https://github.com/adamewing/tldr). LongcallD calls are annotated with [AnnotSV](https://lbgi.fr/AnnotSV/) for comprehensive TE annotation. The workflow supports both somatic (tumor-normal with `--mosaic` mode) and germline calling modes (via `--skip_somatic_te`). Callers can be selected with `--te_calling_tools` (default: `longcalld,tldr`).

### SV annotation

<details markdown="1">
<summary>Output files</summary>

- `annotated_sv/`
  - `*.tsv`: Annotated SVs with gene, OncoKB, and strand information.
- `annotsv/`
  - `*.tsv`: AnnotSV annotation output with regulatory and clinical database annotations.
- `oncokb/`
  - OncoKB cancer gene list (downloaded automatically if not provided).

</details>

SV annotation integrates gene information from [BioMart](https://www.ensembl.org/info/data/biomart/index.html), cancer gene classification from [OncoKB](https://www.oncokb.org/), and comprehensive structural variant annotation from [AnnotSV](https://lbgi.fr/AnnotSV/).

### Visualization

<details markdown="1">
<summary>Output files</summary>

- `plotcircos/`
  - `*.png`: Circos plots showing somatic SVs and CNAs.
- `svkaryoplot/`
  - `*.png`: Karyoplot visualizations for germline SVs.

</details>

Somatic results are visualized as Circos plots integrating SVs and CNAs. Germline SVs are displayed using karyoplots with the [karyoploteR](https://bernatgel.github.io/karyoploter_tutorial/) R package.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. The report includes WhatsHap phasing statistics (from both the phasing workflow and TE calling) and software version tracking. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
