# Plan: Insertion Classification Subworkflow

## Goal

Build a subworkflow to classify somatic insertions from nanomonsv and severus using `nanomonsv insert_classify`, with additional VNTR annotation. This is based on the validated approach from APS055 (SarcAtlas `260315_APS055_somatic_insertions`).

## Background

- The somatic SV pipeline already calls insertions with nanomonsv and severus, merges them via MINDA, and produces an annotated union SV table (`vcf2tsv` → `annotategenes` → final TSV).
- `nanomonsv insert_classify` classifies insertion sequences as L1, Alu, SVA, processed pseudogenes (PSD), or unclassified.
- APS055 showed that ~70% of unclassified insertions are VNTR expansions (overlap hg38 TRF BED). The remaining unclassified are mainly short insertions (≤100 bp) or partial RepeatMasker matches.
- Classification labels are 94% concordant with TLDR on overlapping calls.

## Pipeline Integration Point

The subworkflow should run **after** the ANNOTATE_SV subworkflow produces the annotated union SV table. It will:

1. Filter the union table for INS calls from nanomonsv/severus
2. Resolve symbolic `<INS>` sequences from nanomonsv result table and severus VCF
3. Reformat into nanomonsv's expected input format
4. Run `nanomonsv insert_classify`
5. Intersect unclassified insertions with VNTR BED to label VNTR expansions
6. Produce a final classified insertions table

## Proposed Modules

### 1. `PREP_INSERT_TABLE` (new local module)

**Purpose:** Extract insertions from the annotated union SV table, resolve `<INS>` symbolic alleles, and reformat for `nanomonsv insert_classify`.

**Inputs:**

- `[meta, annotated_sv_tsv]` — union annotated SV table from CSVTK_CONCAT
- `[meta, nanomonsv_result_txt]` — raw nanomonsv result table (has `Inserted_Seq`)
- `[meta, severus_vcf]` — severus somatic VCF (has ALT sequences for insertions)

**Outputs:**

- `[meta, inserts_tsv]` — reformatted insertion table ready for `insert_classify`

**Logic (Python script):**

- Filter union table for `SV_Type == "INS"` with nanomonsv or severus callers
- Extract nanomonsv SV_ID from `SV_callers` column → join to nanomonsv result table for `Inserted_Seq`
- For remaining `<INS>`, extract severus ID → look up ALT from severus VCF
- Reformat columns: `chrom1, base1, Dir_1, chrom2, base2, Dir_2, Alt_Seq, minda_ID, ...`
- Set `base2 = base1 + 1` (nanomonsv convention for insertions)

### 2. `NANOMONSV_INSERT_CLASSIFY` (new local module)

**Purpose:** Run `nanomonsv insert_classify` on the prepared insertion table.

**Container:** `biocontainers/nanomonsv:0.9.0--pyhdfd78af_0` (same as existing nanomonsv modules)

**Inputs:**

- `[meta, inserts_tsv]` — prepared insertion table
- `ref_fasta` — reference genome FASTA
- `ref_gtf` — gene annotation GTF (e.g., gencode v45)
- `line1_db` — LINE1 database BED (nanomonsv resource)

**Outputs:**

- `[meta, classified_tsv]` — classified insertion table with `Insert_Type`, `RMSK_Info`, etc.

**Command:**

```bash
nanomonsv insert_classify \
    ${inserts_tsv} \
    ${output_tsv} \
    ${ref_fasta} \
    ${ref_gtf} \
    ${line1_db}
```

**Note:** nanomonsv requires the RepeatMasker `famdb` directory mounted at `/usr/local/share/RepeatMasker/Libraries/famdb`. This needs to be handled via container bind mounts or baked into the container.

### 3. `CLASSIFY_VNTR_INSERTS` (new local module)

**Purpose:** Intersect unclassified insertions with VNTR BED and assign comprehensive categories.

**Inputs:**

- `[meta, classified_tsv]` — output from `insert_classify`
- `vntr_bed` — hg38 TRF VNTR BED file (already a pipeline parameter)

**Outputs:**

- `[meta, final_classified_tsv]` — final table with `Summary_Category` column

**Logic (Python script):**
Priority hierarchy for assignment:

1. Confirmed `Insert_Type` from nanomonsv (L1, Alu, SVA, PSD, etc.)
2. Breakpoint overlaps VNTR → "VNTR"
3. `RMSK_Info` populated → "Partial RMSK"
4. `SV_LEN` ≤ 100 bp → "Short (≤100 bp)"
5. Everything else → "Other unclassified"

## Proposed Subworkflow: `INSERT_CLASSIFY`

**Location:** `subworkflows/local/insert_classify/main.nf`

```
include { PREP_INSERT_TABLE         } from '../../../modules/local/prep_insert_table/main'
include { NANOMONSV_INSERT_CLASSIFY } from '../../../modules/local/nanomonsv/insert_classify/main'
include { CLASSIFY_VNTR_INSERTS     } from '../../../modules/local/classify_vntr_inserts/main'

workflow INSERT_CLASSIFY {
    take:
    annotated_sv_ch    // [meta, annotated_sv_tsv]
    nanomonsv_result   // [meta, nanomonsv_result_txt]
    severus_vcf        // [meta, severus_vcf]
    ref_fasta
    ref_gtf
    line1_db
    vntr_bed

    main:
    PREP_INSERT_TABLE(annotated_sv_ch, nanomonsv_result, severus_vcf)
    NANOMONSV_INSERT_CLASSIFY(PREP_INSERT_TABLE.out.inserts_tsv, ref_fasta, ref_gtf, line1_db)
    CLASSIFY_VNTR_INSERTS(NANOMONSV_INSERT_CLASSIFY.out.classified_tsv, vntr_bed)

    emit:
    classified_inserts = CLASSIFY_VNTR_INSERTS.out.final_classified_tsv
    versions = ...
}
```

## New Parameters Needed

| Parameter            | Description                                                | Default |
| -------------------- | ---------------------------------------------------------- | ------- |
| `--classify_inserts` | Enable insertion classification subworkflow                | `false` |
| `--ref_gtf`          | Gene annotation GTF for insert_classify                    | `null`  |
| `--line1_db`         | LINE1 database BED (nanomonsv resource)                    | `null`  |
| `--famdb_dir`        | RepeatMasker famdb directory (bind-mounted into container) | `null`  |
| `--vntr_bed`         | VNTR BED file (already exists)                             | `null`  |

## Channel Wiring in `nanogenome.nf`

The subworkflow needs access to:

- The annotated SV table output from `ANNOTATE_SV` (already `CSVTK_CONCAT.out.csv`)
- The raw nanomonsv result table → need to emit `NANOMONSV_GET.out.result_table` from `SV_CALLING_SOMATIC`
- The severus somatic VCF → already emitted as `SV_CALLING_SOMATIC.out.severus_vcf`

**Required change to `SV_CALLING_SOMATIC`:** Emit `nanomonsv_result` (the `.result.txt` table, not just the VCF).

## famdb Handling

`nanomonsv insert_classify` requires the RepeatMasker famdb files. Options:

1. **Bake into container** — build a custom container with famdb pre-installed (simplest)
2. **Bind mount parameter** — add a `--famdb_dir` parameter and bind mount into the container at the expected path
3. **Download at runtime** — use a setup process to download famdb (slow, not recommended)

**Recommendation:** Option 1 (custom container) for simplicity. The existing nanomonsv container doesn't include famdb so we'll need a new one (e.g., `nanomonsv_classify:0.9.0`).

## Implementation Order

1. Add `nanomonsv_result` emit to `SV_CALLING_SOMATIC`
2. Build `PREP_INSERT_TABLE` module + Python script in `bin/`
3. Build `NANOMONSV_INSERT_CLASSIFY` module (+ container with famdb)
4. Build `CLASSIFY_VNTR_INSERTS` module + Python script in `bin/`
5. Wire subworkflow `INSERT_CLASSIFY`
6. Add parameters to `nextflow.config` and `nextflow_schema.json`
7. Wire into `nanogenome.nf` (after ANNOTATE_SV, gated by `--classify_inserts`)
8. Add nf-test

## Decisions

1. **Single-end breakpoints:** Will be tackled as a follow-up after this subworkflow is complete. We'll use nanomonsv docs to guide that implementation.
2. **Scope:** Somatic insertions only for now.
3. **famdb handling:** Use the bind mount approach (same as APS055 pilot). Add a `--famdb_dir` parameter; configure `containerOptions` in `modules.config` to bind it to `/usr/local/share/RepeatMasker/Libraries/famdb`. This keeps the container lightweight and decouples the large famdb files.

## Standalone Execution Mode

The subworkflow should be runnable on its own (like TE calling) for cases where somatic SVs have already been called and annotated. This requires:

- A new boolean parameter: `--classify_inserts` (default: `false`)
- When `--classify_inserts` is set and `--skip_somatic` is true, the subworkflow reads pre-computed inputs from the samplesheet
- New samplesheet columns (optional): `annotated_sv_tsv`, `nanomonsv_result`, `severus_vcf`
- A standalone test config: `conf/test_insert_classify_only.config`

**Standalone invocation:**

```bash
nextflow run shahcompbio/nanogenome \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --classify_inserts \
    --skip_phasing \
    --skip_somatic \
    --skip_cna \
    --fasta ref.fa \
    --fai ref.fa.fai \
    --ref_gtf gencode.gtf \
    --line1_db LINE1.hg38.bed.gz \
    --famdb_dir /path/to/famdb \
    --vntr_bed human_GRCh38.trf.bed
```

When running as part of the full pipeline (somatic SV calling enabled), the subworkflow automatically receives its inputs from upstream processes and no extra samplesheet columns are needed.

## Implementation Status (as of 2025-05-19)

All modules, scripts, subworkflow, config, and workflow wiring are **complete**. Summary of what was built:

### Files Created

| File                                             | Description                                             |
| ------------------------------------------------ | ------------------------------------------------------- |
| `modules/local/prepinserttable/main.nf`          | PREPINSERTTABLE process                                 |
| `modules/local/nanomonsv/insertclassify/main.nf` | NANOMONSV_INSERTCLASSIFY process                        |
| `modules/local/vntrclassify/main.nf`             | VNTRCLASSIFY process                                    |
| `subworkflows/local/insertclassify/main.nf`      | INSERTCLASSIFY subworkflow (chains the 3 modules)       |
| `bin/prep_insert_table.py`                       | Extracts INS from union table, resolves `<INS>` alleles |
| `bin/vntr_classify.py`                           | VNTR intersection + priority classification             |
| `tests/insert_classify_only.nf.test`             | Stub nf-test for standalone mode                        |

### Files Modified

| File                                            | Changes                                                                                                                       |
| ----------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `workflows/nanogenome.nf`                       | Added INSERTCLASSIFY import + block gated by `params.classify_inserts`; updated ALL `ch_samplesheet` closures from 6→8 fields |
| `subworkflows/local/sv_calling_somatic/main.nf` | Emit `nanomonsv_result` channel                                                                                               |
| `subworkflows/local/phasing/main.nf`            | Updated closures from 6→8 fields                                                                                              |
| `nextflow.config`                               | Added params: `classify_inserts`, `ref_gtf`, `line1_db`, `famdb_dir`                                                          |
| `conf/modules.config`                           | Added containerOptions for famdb bind mount + publishDir for all 3 modules                                                    |
| `assets/schema_input.json`                      | Added optional columns: `severus_vcf`, `annotated_sv_tsv`, `nanomonsv_result_txt`                                             |

### Remaining TODO: Test Data for `tests/insert_classify_only.nf.test`

The nf-test uses `-stub` so files only need to **exist** at the URLs — content doesn't matter. Need to push dummy files to `apsteinberg/test-datasets` on the `nanogenome` branch:

1. **`samplesheet_insert_classify.csv`** — samplesheet with columns:

   ```
   sample,condition,bam,bai,snp_vcf,snp_tbi,severus_vcf,annotated_sv_tsv,nanomonsv_result_txt
   COLO829,tumor,https://raw.githubusercontent.com/apsteinberg/test-datasets/nanogenome/bams/COLO829.tumor.haplotagged.chr22.bam,https://raw.githubusercontent.com/apsteinberg/test-datasets/nanogenome/bams/COLO829.tumor.haplotagged.chr22.bam.bai,,,https://raw.githubusercontent.com/apsteinberg/test-datasets/nanogenome/insert_classify/dummy_severus.vcf,https://raw.githubusercontent.com/apsteinberg/test-datasets/nanogenome/insert_classify/dummy_annotated_sv.tsv,https://raw.githubusercontent.com/apsteinberg/test-datasets/nanogenome/insert_classify/dummy_nanomonsv.result.txt
   COLO829,normal,https://raw.githubusercontent.com/apsteinberg/test-datasets/nanogenome/bams/COLO829.normal.haplotagged.chr22.bam,https://raw.githubusercontent.com/apsteinberg/test-datasets/nanogenome/bams/COLO829.normal.haplotagged.chr22.bam.bai,,,,,,
   ```

2. **Dummy input files** (in `insert_classify/` directory):
   - `dummy_severus.vcf` — minimal VCF header + 1 dummy INS record
   - `dummy_annotated_sv.tsv` — TSV with correct column headers + 1 dummy row
   - `dummy_nanomonsv.result.txt` — nanomonsv result header + 1 dummy row

3. **Dummy reference files** (only if not already present):
   - `reference/gencode.v45.chr22.annotation.gtf` — minimal GTF (header line only is fine for stub)
   - `reference/LINE1.hg38.chr22.bed.gz` — gzipped empty/minimal BED

   **Already exist in repo:** `reference/GRCh38_chr22.fa`, `reference/GRCh38_chr22.fa.fai`, `reference/human_GRCh38_no_alt_analysis_set.trf.bed`

4. **Local clone location:** `/Users/preskaa/VSCodeProjects/test-datasets` (currently on `proteomegenerator3` branch — switch to `nanogenome`)

### Other Notes

- The `_` prefix convention on unused closure params triggers Nextflow linting warnings but is functionally correct
- Deprecated `Channel` → `channel` warnings are pre-existing throughout the codebase (not introduced by this work)
- `meta.yml` files for the new modules have placeholder content from `nf-core modules create` scaffold

## Future Work

- Classify single-end breakpoints (`sbnd.result.txt`) using nanomonsv's breakpoint classification modules
- Extend to germline insertions if needed

## Design Decisions

- **SAVANA insertions excluded from classification:** As of SAVANA v1.3.7, the `inserted_sequences.fa` output contains multiple per-read supporting sequences for each insertion variant rather than a single consensus sequence. The `nanomonsv insert_classify` tool expects a single resolved insertion sequence per variant, so SAVANA-only insertions cannot be classified with the current approach. Insertions called by both SAVANA and another caller (nanomonsv or severus) are still included, since the sequence is resolved from the other caller.
