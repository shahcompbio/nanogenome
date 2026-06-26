#!/usr/bin/env python
"""
Prepare insertion table for nanomonsv insert_classify.

Extracts insertions from the annotated union SV table, resolves symbolic <INS>
sequences from nanomonsv result table and severus VCF, and reformats into
nanomonsv's expected input format.
"""

import argparse
import pandas as pd
import re


def parse_severus_vcf(vcf_path):
    """Parse severus VCF to extract insertion sequences keyed by ID."""
    seq_lookup = {}
    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = fields[1]
            sv_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            # Only keep actual sequences (not symbolic <INS>)
            if not alt.startswith("<"):
                seq_lookup[sv_id] = alt
    return seq_lookup


def extract_nanomonsv_id(callers):
    """Extract nanomonsv SV_ID (e.g. 'i_52') from SV_callers column."""
    for caller in str(callers).split(","):
        caller = caller.strip()
        if caller.startswith("i_"):
            return caller
    return None


def extract_severus_id(callers):
    """Extract severus ID from SV_callers column."""
    for caller in str(callers).split(","):
        caller = caller.strip()
        if caller.startswith("severus_"):
            return caller
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Prepare insertion table for nanomonsv insert_classify"
    )
    parser.add_argument(
        "--annotated_sv",
        required=True,
        help="Annotated union SV table from CSVTK_CONCAT",
    )
    parser.add_argument(
        "--nanomonsv_result",
        required=True,
        help="Raw nanomonsv result table (.result.txt)",
    )
    parser.add_argument(
        "--severus_vcf",
        required=True,
        help="Severus somatic VCF",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV for nanomonsv insert_classify",
    )
    args = parser.parse_args()

    # Load annotated SV table
    sv_df = pd.read_csv(args.annotated_sv, sep="\t")

    # Filter for insertions from nanomonsv or severus
    insert_df = sv_df[
        (sv_df["SV_Type"] == "INS")
        & (
            sv_df["SV_callers"].str.contains("i_", na=False)
            | sv_df["SV_callers"].str.contains("severus", na=False)
        )
    ].copy()

    if len(insert_df) == 0:
        # Write empty output with correct columns
        empty_df = pd.DataFrame(
            columns=[
                "Chr_1",
                "Pos_1",
                "Dir_1",
                "Chr_2",
                "Pos_2",
                "Dir_2",
                "Inserted_Seq",
                "minda_ID",
                "SV_Type",
                "SV_callers",
                "SV_LEN",
                "gene_name_1",
                "oncokb_gene1",
            ]
        )
        empty_df.to_csv(args.output, sep="\t", index=False)
        return

    # Load nanomonsv result table for insertion sequences
    nanomonsv_df = pd.read_csv(args.nanomonsv_result, sep="\t")
    nanomonsv_seq_lookup = {}
    if "Inserted_Seq" in nanomonsv_df.columns:
        nanomonsv_seq_lookup = nanomonsv_df.set_index("SV_ID")["Inserted_Seq"].to_dict()

    # Parse severus VCF for insertion sequences
    severus_seq_lookup = parse_severus_vcf(args.severus_vcf)

    # Extract caller IDs
    insert_df["nanomonsv_id"] = insert_df["SV_callers"].apply(extract_nanomonsv_id)
    insert_df["severus_id"] = insert_df["SV_callers"].apply(extract_severus_id)

    # Resolve <INS> from nanomonsv first
    mask_nanomonsv = (insert_df["Alt_Seq"] == "<INS>") & insert_df[
        "nanomonsv_id"
    ].notna()
    insert_df.loc[mask_nanomonsv, "Alt_Seq"] = insert_df.loc[
        mask_nanomonsv, "nanomonsv_id"
    ].map(nanomonsv_seq_lookup)

    # Fill remaining <INS> from severus
    mask_severus = (insert_df["Alt_Seq"] == "<INS>") & insert_df["severus_id"].notna()
    insert_df.loc[mask_severus, "Alt_Seq"] = insert_df.loc[
        mask_severus, "severus_id"
    ].map(severus_seq_lookup)

    # Reformat for nanomonsv insert_classify
    # Set Dir_1 and Dir_2 from orientation column if available
    if "orientation" in insert_df.columns:
        insert_df["Dir_1"] = insert_df["orientation"].str[0]
        insert_df["Dir_2"] = insert_df["orientation"].str[1]
    else:
        insert_df["Dir_1"] = "+"
        insert_df["Dir_2"] = "-"

    # Set Pos_2 = Pos_1 + 1 (nanomonsv convention for insertions)
    insert_df["Pos_2"] = insert_df["base1"] + 1

    # Select and rename columns to match nanomonsv expected format
    output_df = insert_df[
        [
            "chrom1",
            "base1",
            "Dir_1",
            "chrom2",
            "Pos_2",
            "Dir_2",
            "Alt_Seq",
            "minda_ID",
            "SV_Type",
            "SV_callers",
            "SV_LEN",
            "gene_name_1",
            "oncokb_gene1",
        ]
    ].copy()
    output_df.columns = [
        "Chr_1",
        "Pos_1",
        "Dir_1",
        "Chr_2",
        "Pos_2",
        "Dir_2",
        "Inserted_Seq",
        "minda_ID",
        "SV_Type",
        "SV_callers",
        "SV_LEN",
        "gene_name_1",
        "oncokb_gene1",
    ]

    # Drop rows where we couldn't resolve the insertion sequence
    output_df = output_df[
        output_df["Inserted_Seq"].notna() & (output_df["Inserted_Seq"] != "<INS>")
    ]

    output_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
