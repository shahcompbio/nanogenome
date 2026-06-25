# Plan: Single Breakend (SBND) Classification Subworkflow

## Goal

Build a `SBND_CLASSIFY` subworkflow that annotates and classifies single-breakend SVs (`.sbnd.result.txt`) from `nanomonsv get`, triggered alongside the existing `INSERTCLASSIFY` subworkflow when `--classify_inserts` is set.

## Background

- `nanomonsv get` (with `--single_bnd`) already produces a `.sbnd.result.txt` file. The `NANOMONSV_GET` module already emits `sbnd_table`, but `SV_CALLING_SOMATIC` does not yet forward it.
- The upstream `postprocess_sbnd.sh` script wraps two steps:
  1. `annotate_contig.py` — generates a FASTA from sbnd contigs, aligns with BWA, annotates with RepeatMasker, then classifies each contig into: `Simple_repeat`, `Satellite`, `Plain_SV`, `L1_Mediated_Del`, or `Complex`
  2. `plot_contig.R` — generates per-contig visualizations (deferred to future work)
- The goal is to split step 1 into two Nextflow modules (annotation and classification) with reusable Python scripts, wired into a new subworkflow.

## Module Split Rationale

| Step                                               | Module                    | Why separate?                                                                            |
| -------------------------------------------------- | ------------------------- | ---------------------------------------------------------------------------------------- |
| BWA alignment + RepeatMasker annotation of contigs | `NANOMONSV_ANNOTATE_SBND` | Compute-heavy, needs external tools (BWA, RepeatMasker + famdb), cacheable independently |
| Classification using bwa.txt + rmsk.txt            | `NANOMONSV_CLASSIFY_SBND` | Lightweight Python, can be re-run with different logic without re-running BWA/RM         |

---

## Proposed Modules

### 1. `NANOMONSV_ANNOTATE_SBND`

**Location:** `modules/local/nanomonsv/annotateSbnd/main.nf`

**Purpose:** Generate a FASTA from sbnd contig sequences, align with BWA, annotate with RepeatMasker. Produces intermediate annotation files consumed by the classifier.

**Inputs:**

- `[meta, sbnd_result_txt]` — raw sbnd result table from `NANOMONSV_GET`
- `ref_fasta` — reference genome FASTA (must be BWA-indexed; index files staged alongside)
- `ref_fai` — reference genome index

**Outputs:**

- `[meta, bwa_txt]` — BWA alignment results (`.nanomonsv.bwa.txt`)
- `[meta, rmsk_txt]` — RepeatMasker annotation results (`.nanomonsv.rmsk.txt`)
- `[meta, sbnd_result_txt]` — pass-through for downstream module

**Script:** `bin/annotate_sbnd_contigs.py` (split from `annotate_contig.py`; does FASTA generation + BWA + RepeatMasker parsing only, writes `.bwa.txt` and `.rmsk.txt`)

**Command (in module):**

```bash
python3 annotate_sbnd_contigs.py ${sbnd_result_txt} ${ref_fasta} ${prefix}
```

**Container:** conda `environment.yml` with `nanomonsv`, `bwa`, `repeatmasker`, `pysam` — built via Wave at runtime (see Container Strategy section). Docker fallback: `quay.io/biocontainers/` image once built.

**famdb:** Bind-mount via `containerOptions` in `conf/modules.config` (same pattern as `NANOMONSV_INSERTCLASSIFY`):

```groovy
withName: 'NANOMONSV_ANNOTATE_SBND' {
    containerOptions = { params.famdb_dir ? "--bind ${params.famdb_dir}:/usr/local/share/RepeatMasker/Libraries/famdb" : "" }
    publishDir = [...]
}
```

**BWA index:** The module requires the reference to be pre-indexed with BWA. Since no `bwa/index` nf-core module is currently installed:

- Install `nf-core modules install bwa/index` and run it once in the subworkflow (result is cached by Nextflow)
- The subworkflow passes the resulting index directory to `NANOMONSV_ANNOTATE_SBND`

---

### 2. `NANOMONSV_VISUALIZE_SBND`

**Location:** `modules/local/nanomonsv/visualizeSbnd/main.nf`

**Purpose:** Generate per-contig PDF visualizations showing BWA alignment and RepeatMasker annotations. Runs by default when `--classify_inserts` is set; skipped with `--skip_sbnd_vis`.

**Inputs:**

- `[meta, bwa_txt, rmsk_txt]` — annotation outputs from `NANOMONSV_ANNOTATE_SBND`

**Outputs:**

- `[meta, vis_dir]` — directory of per-contig PDFs (`*.nanomonsv.sbnd_vis/`)

**Script:** `bin/plot_sbnd_contigs.R` — the existing `plot_contig.R` with minor refactoring (already accepts prefix + output_dir as args; no logic changes needed)

**Command (in module):**

```bash
Rscript plot_sbnd_contigs.R ${prefix} ${prefix}.nanomonsv.sbnd_vis
```

**Container:** `quay.io/biocontainers/r-ggrepel:0.9.6--r44hf9963bf_0` (includes R, tidyverse, ggrepel; no Docker Hub pull limits).

---

### 3. `NANOMONSV_CLASSIFY_SBND`

**Location:** `modules/local/nanomonsv/classifySbnd/main.nf`

**Purpose:** Classify sbnd contigs using BWA and RepeatMasker annotation results. Applies the priority classification logic from `annotate_contig.py`.

**Inputs:**

- `[meta, sbnd_result_txt, bwa_txt, rmsk_txt]` — from `NANOMONSV_ANNOTATE_SBND`

**Outputs:**

- `[meta, class_txt]` — classification table (`.class.txt`) with columns: `Contig_ID`, `Contig_Class` (`Simple_repeat` | `Satellite` | `Plain_SV` | `L1_Mediated_Del` | `Complex`), `SV_Key` (resolved canonical SV coordinates for `Plain_SV` and `L1_Mediated_Del`)
- versions emitted via topic channel (e.g. `tuple val("${task.process}"), val('nanomonsv'), eval("nanomonsv --version ..."), topic: versions, emit: versions_nanomonsv_classify_sbnd`)

**Script:** `bin/classify_sbnd_contigs.py` (split from `annotate_contig.py`; reads `.bwa.txt`, `.rmsk.txt`, `.sbnd.result.txt` and applies classification logic)

**Command (in module):**

```bash
python3 classify_sbnd_contigs.py ${sbnd_result_txt} ${bwa_txt} ${rmsk_txt} ${prefix}
```

**Container:** `quay.io/biocontainers/nanomonsv:0.9.0--pyhdfd78af_0` (pysam available; no rate-limited Docker Hub)

---

## Proposed Subworkflow: `SBND_CLASSIFY`

**Location:** `subworkflows/local/sbndclassify/main.nf`

```nextflow
include { NANOMONSV_ANNOTATE_SBND } from '../../../modules/local/nanomonsv/annotateSbnd/main'
include { NANOMONSV_CLASSIFY_SBND  } from '../../../modules/local/nanomonsv/classifySbnd/main'
include { NANOMONSV_VISUALIZE_SBND } from '../../../modules/local/nanomonsv/visualizeSbnd/main'

workflow SBND_CLASSIFY {
    take:
    sbnd_result_ch  // channel: [ val(meta), path(sbnd_result_txt) ]
    ref_fasta       // path: reference genome FASTA
    ref_fai         // path: reference genome FAI index
    bwa_index       // path: pre-built BWA index directory (--bwa_index param)

    main:
    // Versions flow automatically via topic channel — no ch_versions bookkeeping needed

    // Annotate sbnd contigs (BWA + RepeatMasker)
    NANOMONSV_ANNOTATE_SBND(
        sbnd_result_ch,
        ref_fasta,
        ref_fai,
        bwa_index
    )

    // Classify contigs
    NANOMONSV_CLASSIFY_SBND(NANOMONSV_ANNOTATE_SBND.out.annotations)

    // Visualize contigs (default on; skipped with --skip_sbnd_vis)
    if (!params.skip_sbnd_vis) {
        NANOMONSV_VISUALIZE_SBND(NANOMONSV_ANNOTATE_SBND.out.annotations)
    }

    emit:
    sbnd_classes = NANOMONSV_CLASSIFY_SBND.out.class_txt // channel: [ val(meta), path(class_txt) ]
    sbnd_vis     = params.skip_sbnd_vis ? Channel.empty() : NANOMONSV_VISUALIZE_SBND.out.vis_dir
}
```

---

## New Python Scripts (`bin/`)

### `bin/annotate_sbnd_contigs.py`

Refactored from `annotate_contig.py` — handles only annotation steps:

1. Parse `.sbnd.result.txt` → write temp FASTA of contig sequences
2. Run `RepeatMasker -species human <fasta>` → parse output → write `{prefix}.nanomonsv.rmsk.txt`
3. Run `bwa mem -h 200 <ref> <fasta>` → parse SAM → write `{prefix}.nanomonsv.bwa.txt`
4. Clean up temp files

### `bin/classify_sbnd_contigs.py`

Refactored from `annotate_contig.py` — handles only classification:

1. Load contig sequences from `.sbnd.result.txt`
2. Load `.rmsk.txt` and `.bwa.txt`
3. For each contig, apply priority classification:
   - `Simple_repeat` / `Satellite` — if ≥80% of contig covered by RepeatMasker simple/satellite
   - `Plain_SV` — if early BWA alignment is long (≥2000 bp) and MQ ≥40
   - `L1_Mediated_Del` — if L1HS/L1P1/L1PA2 covers ≥80% of early BWA alignments, followed by a long high-MQ alignment
   - `Complex` — everything else
4. Write `{prefix}.class.txt`

---

## Pipeline Integration Changes

### `subworkflows/local/sv_calling_somatic/main.nf`

- Add `ch_nanomonsv_sbnd = Channel.empty()` initialization
- Assign `ch_nanomonsv_sbnd = NANOMONSV_GET.out.sbnd_table`
- Add to `emit:` block: `nanomonsv_sbnd = ch_nanomonsv_sbnd`

### `workflows/nanogenome.nf`

- Add `include { SBND_CLASSIFY } from '../subworkflows/local/sbndclassify/main'`
- Within the `if (params.classify_inserts)` block, add validation that `--famdb_dir` is set (already required for INSERTCLASSIFY)
- Wire SBND_CLASSIFY using `SV_CALLING_SOMATIC.out.nanomonsv_sbnd` (or samplesheet fallback in standalone mode)
- No `ch_versions.mix()` needed — all three new modules emit versions via the `versions` topic channel, which flows automatically to MultiQC

### `conf/modules.config`

- Add `withName: 'NANOMONSV_ANNOTATE_SBND'` with famdb `containerOptions` and `publishDir`
- Add `withName: 'NANOMONSV_CLASSIFY_SBND'` with `publishDir`

### `nextflow.config`

- Add `skip_sbnd_vis = false` and `bwa_index = null`
- `--famdb_dir`, `--classify_inserts`, and `--fasta`/`--fai` already cover requirements

### `nextflow_schema.json`

- Add `skip_sbnd_vis` (boolean, default false) and `bwa_index` (string, nullable)

---

## Standalone Mode

When `--classify_inserts` and `--skip_somatic` are both set, the subworkflow reads pre-computed sbnd results from the samplesheet. This requires:

- A new optional samplesheet column: `sbnd_result_txt`
- Logic in `nanogenome.nf` (parallel to existing standalone classify_inserts logic) to route samplesheet values into the SBND_CLASSIFY input channel
- Update `assets/schema_input.json` to add the optional `sbnd_result_txt` column

**Standalone invocation:**

```bash
nextflow run shahcompbio/nanogenome \
    -profile singularity \
    --input samplesheet.csv \  # include sbnd_result_txt column
    --outdir results \
    --classify_inserts \
    --skip_phasing \
    --skip_somatic \
    --fasta ref.fa \
    --fai ref.fa.fai \
    --famdb_dir /path/to/famdb
```

---

## Container Strategy

### `NANOMONSV_ANNOTATE_SBND`

The annotation module needs BWA + RepeatMasker alongside nanomonsv. Since these are unlikely to be in the existing nanomonsv container, use a **conda `environment.yml`** with Wave auto-build:

```yaml
# modules/local/nanomonsv/annotateSbnd/environment.yml
channels:
  - conda-forge
  - bioconda
dependencies:
  - nanomonsv=0.9.0
  - bwa
  - repeatmasker
  - pysam
```

Wave (enabled via `wave { enabled = true }` in the Seqera/cloud profile) builds an OCI container on-the-fly from this file at run time — no manual Docker build required. For HPC/offline use, pre-build with `wave --conda-file environment.yml` and push to `quay.io`.

### `NANOMONSV_CLASSIFY_SBND`

Pure Python + pysam only — use existing nanomonsv container:

```groovy
container "... ? 'https://depot.galaxyproject.org/singularity/nanomonsv:0.9.0--pyhdfd78af_0' : 'quay.io/biocontainers/nanomonsv:0.9.0--pyhdfd78af_0'"
```

### `NANOMONSV_VISUALIZE_SBND`

R + tidyverse + ggrepel — use rocker image from quay.io:

```groovy
container "... ? 'https://depot.galaxyproject.org/singularity/r-ggrepel:0.9.6' : 'quay.io/biocontainers/r-ggrepel:0.9.6--r44hf9963bf_0'"
```

(This BioContainers image includes R, tidyverse, and ggrepel.)

### General container rule

All new modules should use `quay.io/biocontainers/` as the Docker fallback (not Docker Hub `biocontainers/`) to avoid Docker Hub pull rate limits.

---

## Implementation Order

1. Verify nanomonsv container for BWA/RepeatMasker (informs container decision)
2. Write `bin/annotate_sbnd_contigs.py` (split from `annotate_contig.py`)
3. Write `bin/classify_sbnd_contigs.py` (split from `annotate_contig.py`)
4. Build `NANOMONSV_ANNOTATE_SBND` module
5. Build `NANOMONSV_CLASSIFY_SBND` module
6. Copy/refactor `plot_contig.R` → `bin/plot_sbnd_contigs.R`
7. Build `NANOMONSV_VISUALIZE_SBND` module
8. Build `SBND_CLASSIFY` subworkflow
9. Update `SV_CALLING_SOMATIC` to emit `nanomonsv_sbnd`
10. Update `nanogenome.nf` to wire SBND_CLASSIFY (full pipeline + standalone modes)
11. Update `assets/schema_input.json` for `sbnd_result_txt` column
12. Update `conf/modules.config` (publishDir + containerOptions for all 3 modules)
13. Add `skip_sbnd_vis` to `nextflow.config` and `nextflow_schema.json`
14. Add nf-test (stub mode) for `tests/sbnd_classify_only.nf.test`

---

## Deferred / Future Work

- Germline single-breakend SVs (somatic only for now)
- Extend standalone test to use real sbnd data (post-stub)
