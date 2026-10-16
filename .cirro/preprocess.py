#!/usr/bin/env python3
"""
Cirro preprocess for shahcompbio/nanogenome.

Cirro model:  one Cirro sample  = one BAM (+ its BAI)
nanogenome:   one samplesheet `sample` = one patient; tumor and normal rows share it

So each Cirro sample needs two metadata columns:
  - condition : 'tumor' or 'normal'
  - patient   : pair ID shared by a tumor and its matched normal
                (falls back to the Cirro sample name, which suits germline-only runs)
"""

import pandas as pd
from cirro.helpers.preprocess_dataset import PreprocessDataset


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

    meta = ds.samplesheet.copy()
    meta.columns = [c.strip() for c in meta.columns]
    if "condition" not in meta.columns:
        raise ValueError("Sample metadata must include a 'condition' column (tumor/normal).")
    if "patient" not in meta.columns:
        meta["patient"] = meta["sample"]
    meta["patient"] = meta["patient"].fillna(meta["sample"]).astype(str).str.strip()
    meta["condition"] = meta["condition"].astype(str).str.strip().str.lower()

    bad = meta[~meta["condition"].isin(["tumor", "normal"])]
    if len(bad):
        raise ValueError("'condition' must be 'tumor' or 'normal':\n" + bad.to_string())

    sheet = wide.merge(meta[["sample", "patient", "condition"]], on="sample", how="left")
    if sheet["condition"].isna().any():
        raise ValueError(
            "No metadata for samples: " + ", ".join(sheet.loc[sheet["condition"].isna(), "sample"])
        )

    # nanogenome's `sample` column is the pair ID
    sheet = sheet.rename(columns={"sample": "cirro_sample"})
    sheet["sample"] = sheet["patient"].str.replace(r"\s+", "_", regex=True)
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


def set_reference_params(ds: PreprocessDataset):
    fasta = ds.params.get("fasta")
    if not fasta:
        raise ValueError("A reference FASTA must be selected.")
    if not ds.params.get("fai"):
        ds.add_param("fai", f"{fasta}.fai", overwrite=True)


if __name__ == "__main__":
    ds = PreprocessDataset.from_running()

    somatic = set_mode_params(ds)
    set_reference_params(ds)

    sheet = build_samplesheet(ds)
    check_pairs(ds, sheet, somatic)
    sheet.to_csv("samplesheet.csv", index=False)
    ds.logger.info(f"samplesheet.csv:\n{sheet.to_csv(index=False)}")

    ds.add_param("input", "samplesheet.csv", overwrite=True)
    ds.logger.info(f"Final params: {ds.params}")
