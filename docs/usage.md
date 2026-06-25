# shahcompbio/nanogenome: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**nanogenome** is a Nextflow pipeline for variant and mobile element analysis from long-read DNA sequencing data. It supports phasing, somatic and germline structural variant calling, copy number aberration analysis, transposable element calling, and comprehensive annotation. The pipeline is designed to work with Oxford Nanopore Technologies (ONT) aligned BAM files.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location:

```bash
--input '[path to samplesheet file]'
```

The samplesheet must be a comma-separated file with the following columns:

| Column                 | Description                                                                                                                               | Required |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- | -------- |
| `sample`               | Sample identifier. Must be the same for tumor-normal pairs.                                                                               | Yes      |
| `condition`            | Either `tumor` or `normal`.                                                                                                               | Yes      |
| `bam`                  | Full path to aligned BAM file.                                                                                                            | Yes      |
| `bai`                  | Full path to BAM index file (.bai).                                                                                                       | Yes      |
| `snp_vcf`              | Path to pre-phased SNP VCF file (.vcf or .vcf.gz). Required if `--skip_phasing` is used.                                                  | No       |
| `snp_tbi`              | Path to VCF index file (.tbi). Required when `snp_vcf` is provided.                                                                       | No       |
| `severus_vcf`          | Path to pre-computed Severus SV VCF (.vcf or .vcf.gz).                                                                                    | No       |
| `annotated_sv_tsv`     | Path to a pre-computed annotated SV table (.tsv). Required for standalone insertion classification (`--classify_inserts --skip_somatic`). | No       |
| `nanomonsv_result_txt` | Path to a pre-computed NanoMonSV result table (.txt). Required for standalone insertion classification.                                   | No       |
| `sbnd_result_txt`      | Path to a pre-computed NanoMonSV single-breakend result file (.txt). Required for standalone single-breakend classification.              | No       |

### Somatic analysis (tumor-normal pairs)

For somatic SV calling, provide paired tumor and normal samples with matching `sample` identifiers:

```csv title="samplesheet.csv"
sample,condition,bam,bai,snp_vcf,snp_tbi,severus_vcf
SAMPLE1,tumor,/path/to/tumor.bam,/path/to/tumor.bam.bai,,,
SAMPLE1,normal,/path/to/normal.bam,/path/to/normal.bam.bai,,,
```

### Germline analysis

For germline-only analysis (using `--germline`), provide normal samples:

```csv title="samplesheet.csv"
sample,condition,bam,bai,snp_vcf,snp_tbi,severus_vcf
SAMPLE1,normal,/path/to/normal.bam,/path/to/normal.bam.bai,,,
```

### Pre-phased input

If you have already phased your samples and want to skip the phasing step (`--skip_phasing`), provide the phased VCF and its index:

```csv title="samplesheet.csv"
sample,condition,bam,bai,snp_vcf,snp_tbi,severus_vcf
SAMPLE1,tumor,/path/to/tumor.bam,/path/to/tumor.bam.bai,/path/to/phased.vcf.gz,/path/to/phased.vcf.gz.tbi,
SAMPLE1,normal,/path/to/normal.bam,/path/to/normal.bam.bai,/path/to/phased.vcf.gz,/path/to/phased.vcf.gz.tbi,
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run shahcompbio/nanogenome \
   -profile docker \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --fasta <REFERENCE_FASTA> \
   --fai <REFERENCE_FAI>
```

This will launch the default somatic analysis pipeline with phasing, SV calling, and CNA analysis.

### Common workflow modes

#### Full somatic analysis (default)

```bash
nextflow run shahcompbio/nanogenome \
   -profile docker \
   --input samplesheet.csv \
   --outdir results \
   --fasta reference.fa \
   --fai reference.fa.fai
```

#### Germline SV calling

```bash
nextflow run shahcompbio/nanogenome \
   -profile docker \
   --input samplesheet.csv \
   --outdir results \
   --fasta reference.fa \
   --fai reference.fa.fai \
   --germline \
   --skip_somatic
```

#### Transposable element calling

Enable TE calling alongside the standard workflow:

```bash
nextflow run shahcompbio/nanogenome \
   -profile docker \
   --input samplesheet.csv \
   --outdir results \
   --fasta reference.fa \
   --fai reference.fa.fai \
   --te_calling \
   --tldr_te_fasta /path/to/te_reference.fa \
   --longcalld_te_fasta /path/to/te_reference.fa
```

By default, TE calling runs in somatic mode on tumor-normal pairs. To run germline TE calling on normal samples only:

```bash
--te_calling --skip_somatic_te
```

You can select which TE callers to run (default: `longcalld,tldr`):

```bash
--te_calling --te_calling_tools "longcalld"
--te_calling --te_calling_tools "tldr"
--te_calling --te_calling_tools "longcalld,tldr"
```

#### Insertion classification

Classify somatic insertions (L1, Alu, SVA, processed pseudogene, VNTR) using `nanomonsv insert_classify`:

```bash
nextflow run shahcompbio/nanogenome \
   -profile docker \
   --input samplesheet.csv \
   --outdir results \
   --fasta reference.fa \
   --fai reference.fa.fai \
   --classify_inserts \
   --ref_gtf /path/to/gencode.annotation.gtf \
   --line1_db /path/to/LINE1.hg38.bed.gz \
   --vntr_bed /path/to/human_GRCh38_no_alt_analysis_set.trf.bed \
   --bwa_index /path/to/bwa_index_dir
```

> **Note:** Only insertions called by nanomonsv or severus are classified, because these callers provide resolved consensus insertion sequences. SAVANA insertions are excluded because as of SAVANA v1.3.7, it emits multiple per-read supporting sequences rather than a single consensus insertion sequence, which is incompatible with the `nanomonsv insert_classify` input format.

Setting `--classify_inserts` also runs single-breakend (SBND) classification automatically alongside insertion classification. NanoMonSV single-breakend contigs are aligned with BWA, annotated with RepeatMasker, and classified to identify mobile-element-derived single breakends. Per-contig PDF visualizations are generated and merged by default; disable with `--skip_sbnd_vis`.

##### Standalone classification (`--skip_somatic`)

If you already have somatic SV calling results and only want to run classification, set `--skip_somatic` alongside `--classify_inserts` and provide the pre-computed inputs via the samplesheet:

```csv title="samplesheet.csv"
sample,condition,bam,bai,snp_vcf,snp_tbi,severus_vcf,annotated_sv_tsv,nanomonsv_result_txt,sbnd_result_txt
SAMPLE1,tumor,/path/to/tumor.bam,/path/to/tumor.bam.bai,,,,/path/to/annotated_sv.tsv,/path/to/nanomonsv_result.txt,/path/to/sbnd_result.txt
```

`annotated_sv_tsv` and `nanomonsv_result_txt` are required for insertion classification; `sbnd_result_txt` is required for single-breakend classification. Either can be omitted if you only want to run the other.

#### Somatic SNV/indel calling

Enable somatic SNV and indel calling with [ClairS](https://github.com/HKU-BAL/ClairS) or [DeepSomatic](https://github.com/google/deepsomatic):

```bash
nextflow run shahcompbio/nanogenome \
   -profile docker \
   --input samplesheet.csv \
   --outdir results \
   --fasta reference.fa \
   --fai reference.fa.fai \
   --somatic_snv_calling \
   --somatic_snv_caller clairs
```

By default, VCF output is annotated with [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) via [vcf2maf](https://github.com/mskcc/vcf2maf) and converted to MAF format. To skip VEP annotation (e.g., for testing or when a VEP cache is unavailable):

```bash
--somatic_snv_calling --inhibit_vep
```

You can select between callers:

```bash
--somatic_snv_calling --somatic_snv_caller "clairs"
--somatic_snv_calling --somatic_snv_caller "deepsomatic"
```

If you have a local VEP cache, provide it with `--vep_cache`:

```bash
--somatic_snv_calling --vep_cache /path/to/vep_cache
```

#### T2T-CHM13v2.0 genome build

The pipeline supports T2T-CHM13v2.0 for SV annotation and karyoplot visualization via `--genome_build t2t`:

```bash
nextflow run shahcompbio/nanogenome \
   -profile docker \
   --input samplesheet.csv \
   --outdir results \
   --fasta /path/to/chm13v2.0.fa \
   --fai /path/to/chm13v2.0.fa.fai \
   --genome_build t2t
```

When `--genome_build t2t` is set and `--gene_annotations` is not provided, the pipeline generates a T2T gene annotation table using the `T2TGENETABLE` module (powered by the [BiocT2T](https://bioconductor.org/) R package) instead of querying Ensembl BioMart, which does not support T2T-CHM13. Supported values for `--genome_build` are `hg38` (default), `hg19`, and `t2t`.

#### Skip phasing (use pre-phased data)

```bash
nextflow run shahcompbio/nanogenome \
   -profile docker \
   --input samplesheet.csv \
   --outdir results \
   --fasta reference.fa \
   --fai reference.fa.fai \
   --skip_phasing
```

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run shahcompbio/nanogenome -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: "./samplesheet.csv"
outdir: "./results/"
fasta: "/path/to/reference.fa"
fai: "/path/to/reference.fa.fai"
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull shahcompbio/nanogenome
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [shahcompbio/nanogenome releases page](https://github.com/shahcompbio/nanogenome/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
