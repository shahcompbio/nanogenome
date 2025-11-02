# CLAUDE.md

## Project Overview

**nanogenome** is a Nextflow-based bioinformatics pipeline for variant and mobile element analysis from long-read DNA sequencing data. The pipeline supports both somatic and germline structural variant (SV) calling, phasing, copy number aberration (CNA) analysis, and comprehensive annotation of genomic variants.

## Repository Structure

### Core Files

- **[main.nf](main.nf)** - Main pipeline entry point
- **[workflows/nanogenome.nf](workflows/nanogenome.nf)** - Primary workflow orchestration
- **[nextflow.config](nextflow.config)** - Pipeline configuration with default parameters
- **[nextflow_schema.json](nextflow_schema.json)** - Schema validation for pipeline parameters

### Key Directories

#### Subworkflows ([subworkflows/local/](subworkflows/local/))

- **phasing/** - Variant phasing and BAM haplotagging
- **sv_calling_somatic/** - Somatic structural variant calling
- **sv_calling_germline/** - Germline structural variant calling
- **annotate_sv/** - SV merging and annotation

#### Modules ([modules/](modules/))

**Local Modules:**

- **wakhan/** - Copy number aberration analysis (cna, hapcorrect, rephase_cna)
- **whatshap/** - Phasing operations (haplotag, stats)
- **nanomonsv/** - SV calling (get, parse)
- **savana/** - SV classification
- **longcalld/** - Long-read variant calling
- **minda/** - Consensus SV calling
- **plotcircos/** - Circular genome visualization
- **annotategenes/** - Gene annotation
- **biomart/** - BioMart integration
- **svtypes/** - SV type classification
- **vcf2tsv/** - VCF to TSV conversion

**nf-core Modules:**

- **clair3/** - SNV/indel calling
- **longphase/** - Variant phasing
- **severus/** - SV calling
- **sniffles/** - SV detection
- **cutesv/** - SV calling
- **bcftools/** - VCF manipulation
- **samtools/** - BAM operations
- **multiqc/** - Quality control reporting

### Configuration ([conf/](conf/))

- **base.config** - Base process configuration
- **modules.config** - Module-specific settings
- **igenomes.config** - Reference genome configurations
- **test\*.config** - Test profile configurations

### Testing ([tests/](tests/))

- nf-test based unit tests for pipeline components

### Documentation ([docs/](docs/))

- Pipeline documentation and usage guides

## Pipeline Workflow

### 1. Phasing Workflow (Optional)

If `--skip_phasing` is not set:

- Variant calling with Clair3
- Phasing with LongPhase
- BAM haplotagging with WhatsHap
- Phasing statistics generation

### 2. Somatic SV Calling (Default)

If `--skip_somatic` is not set:

- Multiple SV callers: Severus, SAVANA, NanoMonSV
- Consensus calling with MINDA
- CNA analysis with Wakhan
- Rephasing of CNAs based on SV haplotypes

### 3. Germline SV Calling (Optional)

If `--germline` is set:

- Multiple SV callers: Severus, LongCallD, CuteSV, Sniffles
- Consensus calling with MINDA

### 4. Annotation and Visualization

- SV merging across callers
- Gene annotation with BioMart
- OncoKB annotation for cancer genes
- Circos plot generation for somatic variants with CNAs

### 5. Quality Control

- MultiQC report generation
- Software version tracking

## Key Parameters

### Input/Output

- `--input` - Samplesheet with BAM files and optional VCF files
- `--outdir` - Output directory
- `--fasta` - Reference genome FASTA
- `--fai` - Reference genome index

### Workflow Control

- `--skip_phasing` - Skip phasing and use pre-phased data
- `--skip_somatic` - Skip somatic SV calling
- `--germline` - Enable germline SV calling

### Tool-Specific

- `--clair3_model` - Clair3 model (default: r1041_e82_400bps_sup_v500)
- `--clair3_platform` - Platform type (default: ont)
- `--somatic_callers` - Somatic SV callers (default: "severus,savana,nanomonsv")
- `--germline_callers` - Germline SV callers (default: "severus,longcallD,cutesv,sniffles")
- `--min_callers` - Minimum support for consensus SV calls (default: 2)
- `--binsize` - Bin size for CNA analysis (default: 50000)
- `--use_sv_haplotypes` - Use SV haplotypes in CNA analysis (default: true)

### Annotation

- `--gene_annotations` - Gene annotation file
- `--oncokb` - OncoKB gene list
- `--oncokb_url` - OncoKB API URL
- `--genome_build` - Genome build (default: hg38)

## Execution Profiles

- **docker** - Run with Docker containers
- **singularity** - Run with Singularity containers
- **conda** - Run with Conda environments
- **test** - Quick test with small dataset
- **test_full** - Full test dataset
- **test_sv_only** - SV calling only test
- **test_phasing** - Phasing workflow test
- **slurm** - SLURM cluster execution

## Development

### Built With

- Nextflow >= 25.04.0
- nf-core template v3.4.1
- nf-schema plugin v2.5.1
- nf-test for testing

### Author

- Asher Preska Steinberg

### Contributing

See [.github/CONTRIBUTING.md](.github/CONTRIBUTING.md) for contribution guidelines.

## CI/CD

- GitHub Actions for testing ([.github/workflows/](.github/workflows/))
- nf-test for unit testing
- Pre-commit hooks for code quality

## Citations

See [CITATIONS.md](CITATIONS.md) for references to tools used in the pipeline.

## Notes for Claude

### Working with this Repository

1. **Nextflow DSL2**: The pipeline uses Nextflow DSL2 syntax with modules and subworkflows.

2. **Channel Structure**: Key channels carry metadata tuples with structure like:
   - `[meta, bam, bai]` for BAM files
   - `[meta, vcf, tbi]` for VCF files
   - Meta maps contain: `id`, `condition` (tumor/normal), and other sample information

3. **Branching Logic**: The workflow branches based on:
   - Tumor vs Normal samples (`meta.condition`)
   - Somatic vs Germline analysis mode
   - Phasing vs non-phasing mode

4. **Testing**: Use nf-test for testing modules and subworkflows. Test files are in parallel directory structure under [tests/](tests/).

5. **Configuration**: Module-specific configurations are in [conf/modules.config](conf/modules.config). Process selectors use the pattern `withName:PROCESS_NAME`.

6. **Common Tasks**:
   - Adding a new module: Place in [modules/local/](modules/local/) or use nf-core modules
   - Modifying parameters: Update [nextflow.config](nextflow.config) and [nextflow_schema.json](nextflow_schema.json)
   - Adding tests: Create `.nf.test` files alongside modules/workflows

### Architecture Notes

- The pipeline supports multiple entry points: full analysis, phasing only, or SV calling only
- Somatic analysis requires tumor/normal pairs
- Germline analysis works with normal samples only
- CNA analysis integrates SV haplotypes when `use_sv_haplotypes` is true
- Consensus SV calling uses MINDA to merge results from multiple callers
