# shahcompbio/nanogenome: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0] - 2026-01-27

### Added

- AnnotSV integration for comprehensive structural variant annotation with gene, regulatory, and clinical databases
- New `--skip_annotsv` parameter to optionally skip AnnotSV annotation
- New `--annotsv_dir` parameter to specify pre-downloaded AnnotSV annotations directory
- New `--annotsv_genomebuild` parameter to specify genome build for AnnotSV (default: GRCh38)
- TSV2BEDPE module for converting annotated SV TSV files to BEDPE format
- SURVIVOR/BEDPETOVCF module for converting BEDPE to VCF format for AnnotSV input
- Automatic download of AnnotSV annotations if `--annotsv_dir` is not provided

### Changed

- Updated SV annotation workflow to include AnnotSV annotation step
- Updated Minda Docker image
- Modified Wakhan Dockerfile to build from Debian base image for improved compatibility
- Improved vcf2tsv.py with additional SV type handling
- Updated nf-core template to version 3.5.1

### Fixed

- Fixed output file naming in annotation workflow
- Fixed SV type extraction from MINDA output

## [1.1.0] - 2025-11-19

### Added

- Karyoplot visualization for germline structural variants using karyoploteR
- SAVANA as an alternative CNA calling tool alongside Wakhan
- Option to skip CNA analysis with `--skip_cna` parameter
- Configurable CNA tools selection with `--cna_tools` parameter (default: "wakhan,savana")
- Step change parameter for SAVANA CNA calling (`--main_cn_step_change`)
- CNA-only workflow mode for running copy number analysis independently
- Test profiles for CNA-only (`test_cna_only`) and germline-only (`test_germline_only`) workflows
- Stub tests for local modules
- Warning message when VNTR BED file is not provided
- Assertion tests for problematic chromosome filtering

### Changed

- Replaced Circos plots with karyoplot visualization for germline SV analysis
- LongcallD output is now filtered to include only structural variants (excludes small indels)
- Increased SLURM queue size from 50 to 200 for improved throughput
- Limited number of parallel gene annotation processes for resource management
- Updated queue selection logic for SLURM executor

### Fixed

- Fixed Circos plot generation when CNA data is missing
- Fixed karyoplot legend display
- Fixed phasing workflow issues
- Fixed germline-only workflow configuration
- Fixed typos in configuration

### Removed

- Removed ASCAT module (experimental feature)
- Removed Wakhan build script

## [1.0.0] - 2025-11-01

Initial release of shahcompbio/nanogenome, created with the [nf-core](https://nf-co.re/) template.

### Added

#### Core Workflows

- Phasing workflow with Clair3 variant calling and LongPhase phasing
- BAM haplotagging with WhatsHap for haplotype-aware analysis
- Somatic structural variant (SV) calling workflow with tumor-normal paired analysis
- Germline structural variant calling workflow (optional, enabled via `--germline`)
- Haplotype-resolved copy number aberration (CNA) analysis with Wakhan
- SV annotation workflow with gene and oncology knowledge base integration
- Circos plot generation for integrated visualization of SVs and CNAs

#### Somatic SV Calling

- SEVERUS integration for long-read SV detection
- SAVANA integration for deep-learning based somatic SV classification
- NanoMonSV integration for assembly-based SV calling
- MINDA consensus calling to merge results from multiple SV callers
- Configurable minimum caller support (`--min_callers`, default: 2)

#### Germline SV Calling

- SEVERUS support for germline variant detection
- Sniffles integration for germline SV calling
- CuteSV integration for lightweight SV detection
- LongcallD integration for dedicated germline analysis
- MINDA consensus calling across germline callers

#### Copy Number Analysis

- Wakhan integration for haplotype-specific CNA detection
- Automatic ploidy and purity estimation
- Haplotype-resolved copy number segments (HP1 and HP2)
- SV haplotype integration for improved phasing (`--use_sv_haplotypes`)
- Configurable bin size for CNA analysis (`--binsize`, default: 50000)
- Loss of heterozygosity (LOH) detection
- Subclonal CNA detection support

#### Annotation Features

- BioMart integration for automatic gene annotation
- OncoKB integration for cancer gene annotation (tumor suppressor genes and oncogenes)
- Automatic download of annotation databases if not provided
- Support for hg38 and hg19 genome builds
- SV type classification and strand information extraction
- Comprehensive TSV output with gene and OncoKB annotations

#### Workflow Flexibility

- Multiple execution modes:
  - Full pipeline (phasing + SV calling + CNA)
  - Phasing only mode (`--skip_somatic`)
  - SV calling only mode (`--skip_phasing`)
  - Germline analysis mode (`--germline`)
- Support for pre-phased VCF input when skipping phasing
- Configurable SV caller selection via comma-separated lists
- Regional filtering support via BED file (`--filter_bed`)
- Minimum SV size filtering (`--min_size`, default: 50)
- Breakpoint tolerance configuration (`--tolerance`, default: 100bp)

#### Platform Support

- Oxford Nanopore Technologies (ONT) platform support
- Configurable Clair3 models (`--clair3_model`)
- Platform-specific flags for LongPhase, SAVANA, and LongcallD
- VNTR/tandem repeat annotation support (`--vntr_bed`)

#### Quality Control and Reporting

- MultiQC integration for comprehensive QC reporting
- WhatsHap phasing statistics
- Software version tracking across all modules
- Execution timeline, trace, and DAG visualization
- Pipeline information reports

#### Containerization

- Docker container support
- Singularity/Apptainer container support
- Conda environment support
- Podman, Shifter, and Charliecloud support
- Wave container provisioning support

#### Testing

- nf-test based unit testing framework
- Test profile for quick validation (`-profile test`)
- Full test dataset profile (`-profile test_full`)
- SV-only test profile (`-profile test_sv_only`)
- Phasing-only test profile (`-profile test_phasing`)
- GitHub Actions CI/CD integration

#### Documentation

- Comprehensive README with usage examples
- CLAUDE.md for AI-assisted development
- Usage documentation in docs/usage.md
- Output documentation in docs/output.md
- CITATIONS.md with references to all tools

### Dependencies

- Nextflow >= 25.04.0
- nf-schema plugin v2.5.1
- nf-core template v3.4.1
- Clair3 (variant calling)
- LongPhase (phasing)
- WhatsHap (haplotagging and phasing statistics)
- SEVERUS (SV calling)
- SAVANA (somatic SV classification)
- NanoMonSV (assembly-based SV detection)
- Sniffles (germline SV calling)
- CuteSV (germline SV calling)
- LongcallD (germline SV calling)
- MINDA (consensus SV calling)
- Wakhan (CNA analysis)
- BCFtools (VCF manipulation)
- SAMtools (BAM operations)
- Tabix (VCF indexing)
- MultiQC (quality control reporting)
- BioMart (gene annotation)
- R and Python dependencies for custom modules

[unreleased]: https://github.com/shahcompbio/nanogenome/compare/1.2.0...HEAD
[1.2.0]: https://github.com/shahcompbio/nanogenome/compare/1.1.0...1.2.0
[1.1.0]: https://github.com/shahcompbio/nanogenome/compare/1.0.0...1.1.0
[1.0.0]: https://github.com/shahcompbio/nanogenome/releases/tag/1.0.0
