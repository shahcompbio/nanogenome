# T2T Support Plan

## Goal

Add T2T-CHM13 (Telomere-to-Telomere) genome support to the nanogenome pipeline, starting with gene annotation retrieval in `bin/biomart.R`.

## Background

- Ensembl BioMart does **not** support T2T-CHM13. Rapid Release (which hosted T2T) was archived in October 2025 and never had a BioMart endpoint.
- The current `biomart.R` script only supports `hg38` and `hg19` via Ensembl BioMart queries.
- The pipeline already supports a `--gene_annotations` parameter to skip BioMart and supply a pre-built annotation file.

## Plan

### 1. Update `bin/biomart.R` to handle `t2t`

- Add a `t2t` branch that downloads the CHM13 gene annotation GFF3 from NCBI/Ensembl FTP (e.g., NCBI RefSeq GFF3 for GCF_009914755.1 / T2T-CHM13v2.0).
- Parse the GFF3 into the same tab-delimited format as the BioMart output (`ensembl_gene_id`, `hgnc_symbol`, `gene_biotype`, `description`, `chromosome_name`, `start_position`, `end_position`, `strand`).
- Write out `t2t-genes.txt` matching the existing `{genome}-genes.txt` convention.
- Consider using R packages like `rtracklayer` to parse the GFF3, or alternatively use simple text parsing with base R.

### 2. Update the BIOMART Nextflow module (`modules/local/biomart/main.nf`)

- Ensure the container has any additional R packages needed for GFF3 parsing (e.g., `rtracklayer`), or use base R parsing to avoid container changes.
- Update the stub section to handle `t2t` as a valid argument.

### 3. Update pipeline configuration

- In `nextflow.config`: note `t2t` as a valid `genome_build` option.
- In `nextflow_schema.json`: add `t2t` to the allowed values / description for `genome_build`.

### 4. Validate downstream compatibility

- Check that `annotategenes` and other downstream modules that consume `*-genes.txt` work correctly with T2T chromosome names (chr1-chr22, chrX, chrY, chrM, plus T2T-specific unplaced contigs).
- Check `modules/local/svkaryoplot/main.nf` which also takes `genome_build` as input.
- Review `conf/igenomes.config` for any T2T reference genome paths that may need adding.

### 5. Testing

- Add a test profile or test case for T2T support.
- Verify the GFF3 download and parsing produces correct output.

## Key Decisions Still Needed

- **Which FTP source for T2T gene annotations?** NCBI RefSeq GFF3 is the most complete option. Specific URL to confirm.
- **R parsing approach:** `rtracklayer::import()` is cleanest but requires an extra package in the container. Base R parsing avoids container changes but is more code.
- **Chromosome naming:** T2T uses `chr1`-`chr22`, `chrX`, `chrY`, `chrM` — confirm this matches expectations in downstream modules.

## Files to Modify

- `bin/biomart.R` — main changes
- `modules/local/biomart/main.nf` — possible container/dependency update
- `nextflow.config` — document t2t option
- `nextflow_schema.json` — schema validation update
- `conf/modules.config` — if any module-specific config needed
- Downstream annotation modules — validation only

## Reference

- Tell Claude: "Read T2T_PLAN.md and implement the T2T support plan"
