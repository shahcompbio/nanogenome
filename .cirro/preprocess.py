#!/usr/bin/env python3
"""
Cirro preprocess for shahcompbio/nanogenome.

Cirro model:  one Cirro sample  = one BAM (+ its BAI)
nanogenome:   one samplesheet `sample` = one patient; tumor and normal rows share it

Tumor/normal assignment, in order of precedence:
  1. Form fields `tumor_sample` / `normal_sample` (one pair per run; exact sample
     name or a unique substring of it). Optional `patient_id` names the pair.
  2. Sample metadata columns `condition` (tumor/normal) and `patient` (pair ID),
     for multi-pair runs.
  3. Germline-only mode with no metadata: every sample is treated as normal.

Input may come from one or several Cirro datasets; samples are pooled.
"""

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset


def _match(name: str, samples: list, label: str) -> str:
    name = name.strip()
    if name in samples:
        return name
    hits = [x for x in samples if name in x]
    if len(hits) == 1:
        return hits[0]
    raise ValueError(
        f"{label} '{name}' matched {len(hits)} samples; expected exactly one. "
        f"Available samples: {samples}"
    )


def assign_conditions(ds: PreprocessDataset, samples: list) -> pd.DataFrame:
    """Return DataFrame[sample, patient, condition] for the samples to run."""
    tumor = (ds.params.get("tumor_sample") or "").strip()
    normal = (ds.params.get("normal_sample") or "").strip()
    patient = (ds.params.get("patient_id") or "").strip()
    germline_only = ds.params.get("skip_somatic", False)

    # 1. Form fields
    if tumor or normal:
        rows = []
        if tumor:
            rows.append((_match(tumor, samples, "Tumor sample"), "tumor"))
        if normal:
            rows.append((_match(normal, samples, "Normal sample"), "normal"))
        if len({r[0] for r in rows}) != len(rows):
            raise ValueError("Tumor and normal sample fields matched the same sample.")
        pair_id = patient or rows[0][0]
        unused = sorted(set(samples) - {r[0] for r in rows})
        if unused:
            ds.logger.warning(f"Samples not assigned tumor/normal and skipped: {unused}")
        return pd.DataFrame(
            [{"sample": smp, "patient": pair_id, "condition": cond} for smp, cond in rows]
        )

    # 2. Sample metadata
    meta = ds.samplesheet.copy()
    meta.columns = [c.strip() for c in meta.columns]
    if "condition" in meta.columns:
        meta = meta[meta["sample"].isin(samples)].copy()
        if "patient" not in meta.columns:
            meta["patient"] = meta["sample"]
        meta["patient"] = meta["patient"].fillna(meta["sample"]).astype(str).str.strip()
        meta["condition"] = meta["condition"].astype(str).str.strip().str.lower()
        bad = meta[~meta["condition"].isin(["tumor", "normal"])]
        if len(bad):
            raise ValueError("'condition' must be 'tumor' or 'normal':\n" + bad.to_string())
        missing = sorted(set(samples) - set(meta["sample"]))
        if missing:
            raise ValueError(f"No 'condition' metadata for samples: {missing}")
        return meta[["sample", "patient", "condition"]]

    # 3. Germline-only: everything is a normal
    if germline_only:
        return pd.DataFrame({"sample": samples, "patient": samples, "condition": "normal"})

    raise ValueError(
        "Can't tell which sample is tumor and which is normal. Fill in the "
        "'Tumor sample' and 'Normal sample' form fields (sample names: "
        f"{samples}), or add 'condition'/'patient' columns to the sample metadata."
    )


def build_samplesheet(ds: PreprocessDataset) -> pd.DataFrame:
    files = ds.files.copy()
    ds.logger.info(f"Input files:\n{files.to_string()}")

    def kind(path: str):
        p = path.lower()
        if p.endswith(".bam"):
            return "bam"
        if p.endswith(".bai"):  # covers x.bam.bai and x.bai
            return "bai"
        return None

    files["kind"] = files["file"].map(kind)
    files = files.dropna(subset=["kind"])

    dup = files.groupby(["sample", "kind"]).size()
    if (dup > 1).any():
        raise ValueError(
            "Expected exactly one BAM and one BAI per Cirro sample; found duplicates:\n"
            + dup[dup > 1].to_string()
        )

    wide = files.pivot(index="sample", columns="kind", values="file").reset_index()
    for col in ("bam", "bai"):
        if col not in wide.columns:
            wide[col] = None
    missing = wide[wide["bam"].isna() | wide["bai"].isna()]
    if len(missing):
        raise ValueError(
            "Every sample needs both a .bam and a .bai. Missing for:\n" + missing.to_string()
        )

    meta = assign_conditions(ds, list(wide["sample"]))
    sheet = wide.merge(meta, on="sample", how="inner")

    # nanogenome's `sample` column is the pair ID
    sheet = sheet.rename(columns={"sample": "cirro_sample"})
    sheet["sample"] = sheet["patient"].astype(str).str.replace(r"\s+", "_", regex=True)
    # Optional columns (snp_vcf, snp_tbi, severus_vcf, ...) are left out entirely
    return sheet[["sample", "condition", "bam", "bai"]].sort_values(["sample", "condition"])


def check_pairs(ds: PreprocessDataset, sheet: pd.DataFrame, somatic: bool):
    counts = sheet.groupby(["sample", "condition"]).size().unstack(fill_value=0)
    for c in ("tumor", "normal"):
        if c not in counts.columns:
            counts[c] = 0
    if (counts > 1).any().any():
        raise ValueError("More than one tumor or normal BAM per patient:\n" + counts.to_string())
    if somatic:
        unpaired = counts[(counts["tumor"] != 1) | (counts["normal"] != 1)]
        if len(unpaired):
            raise ValueError(
                "Somatic mode needs exactly one tumor and one normal per patient. Unpaired:\n"
                + unpaired.to_string()
            )


def set_mode_params(ds: PreprocessDataset) -> bool:
    mode = ds.params.get("analysis_mode", "somatic")
    somatic = mode in ("somatic", "somatic_and_germline")
    germline = mode in ("germline", "somatic_and_germline")
    ds.add_param("skip_somatic", not somatic, overwrite=True)
    ds.add_param("germline", germline, overwrite=True)
    # Not a nanogenome param; drop so nf-schema doesn't warn about it
    ds.remove_param("analysis_mode", force=True)
    return somatic


def drop_helper_and_empty_params(ds: PreprocessDataset):
    """Remove form-only keys and blank optional fields so nextflow sees null, not ''."""
    for k in ("tumor_sample", "normal_sample", "patient_id"):
        ds.remove_param(k, force=True)
    for k, v in list(ds.params.items()):
        if isinstance(v, str) and not v.strip():
            ds.remove_param(k, force=True)


def check_optional_step_refs(ds: PreprocessDataset):
    """Fail before launch, listing every missing file, instead of one Nextflow error per run."""
    have = lambda k: bool(str(ds.params.get(k) or "").strip())
    problems = []
    if ds.params.get("classify_inserts"):
        miss = [k for k in ("ref_gtf", "line1_db", "bwa_index", "vntr_bed") if not have(k)]
        if miss:
            problems.append(f"'Classify somatic insertions' needs: {', '.join(miss)}")
    if ds.params.get("te_calling"):
        tools = str(ds.params.get("te_calling_tools") or "longcalld,tldr").split(",")
        if "tldr" in tools and not have("tldr_te_fasta"):
            problems.append("'Transposable element calling' (tldr) needs: tldr_te_fasta")
    if problems:
        raise ValueError(
            "Missing reference files for enabled optional steps. Either fill these in "
            "under 'Optional reference files' or turn the step off:\n  - "
            + "\n  - ".join(problems)
        )


def set_reference_params(ds: PreprocessDataset):
    fasta = ds.params.get("fasta")
    if not fasta:
        raise ValueError(
            "No reference FASTA received (params.fasta is empty). "
            "Enter an s3:// or https:// URL in the 'Genome FASTA' form field."
        )
    if not ds.params.get("fai"):
        ds.add_param("fai", f"{fasta}.fai", overwrite=True)


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()
    ds.logger.info(f"Params received from form: {ds.params}")

    somatic = set_mode_params(ds)
    set_reference_params(ds)
    check_optional_step_refs(ds)

    sheet = build_samplesheet(ds)
    check_pairs(ds, sheet, somatic)
    drop_helper_and_empty_params(ds)
    sheet.to_csv("samplesheet.csv", index=False)
    ds.logger.info(f"samplesheet.csv:\n{sheet.to_csv(index=False)}")

    ds.add_param("input", "samplesheet.csv", overwrite=True)
    ds.logger.info(f"Final params: {ds.params}")
