#!/usr/bin/env python
"""
Classify insertions with VNTR annotation and priority hierarchy.

Takes the output of nanomonsv insert_classify and intersects unclassified
insertions with a VNTR BED file, then assigns a Summary_Category using:
1. Confirmed insert_classify type (L1, Alu, SVA, PSD, etc.)
2. Breakpoint overlaps VNTR -> "VNTR"
3. RMSK_Info populated -> "Partial RMSK"
4. SV_LEN <= 100 bp -> "Short (<=100 bp)"
5. Everything else -> "Other unclassified"
"""
import argparse
import pandas as pd
import pybedtools


def assign_summary_category(row, vntr_ids):
    """Assign a single category using priority hierarchy."""
    ins_type = str(row.get("Insert_Type", "---"))
    if ins_type not in ("---", "nan", ""):
        return ins_type

    if row.get("minda_ID") in vntr_ids:
        return "VNTR"

    rmsk_info = str(row.get("RMSK_Info", "---"))
    if rmsk_info not in ("---", "nan", ""):
        return "Partial RMSK"

    sv_len = row.get("SV_LEN", 0)
    try:
        if float(sv_len) <= 100:
            return "Short (<=100 bp)"
    except (ValueError, TypeError):
        pass

    return "Other unclassified"


def main():
    parser = argparse.ArgumentParser(
        description="Classify insertions with VNTR annotation"
    )
    parser.add_argument(
        "--classified_tsv",
        required=True,
        help="Output of nanomonsv insert_classify",
    )
    parser.add_argument(
        "--vntr_bed",
        required=True,
        help="VNTR BED file (e.g., hg38 TRF BED)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV with Summary_Category column",
    )
    args = parser.parse_args()

    # Load classified insertions
    classified_df = pd.read_csv(args.classified_tsv, sep="\t")

    if len(classified_df) == 0:
        classified_df["Summary_Category"] = pd.Series(dtype=str)
        classified_df.to_csv(args.output, sep="\t", index=False)
        return

    # Identify unclassified insertions for VNTR intersection
    unclassified_mask = classified_df["Insert_Type"].isin(["---"]) | classified_df[
        "Insert_Type"
    ].isna()
    unclassified_df = classified_df[unclassified_mask]

    # Intersect unclassified breakpoints with VNTR BED
    vntr_ids = set()
    if len(unclassified_df) > 0:
        vntr_bed = pybedtools.BedTool(args.vntr_bed)
        ins_bed_df = unclassified_df[["Chr_1", "Pos_1", "minda_ID"]].copy()
        ins_bed_df["end"] = ins_bed_df["Pos_1"] + 1
        ins_bed_df = ins_bed_df[["Chr_1", "Pos_1", "end", "minda_ID"]]
        ins_bed_df.columns = ["chrom", "start", "end", "minda_ID"]
        ins_bed = pybedtools.BedTool.from_dataframe(ins_bed_df)

        vntr_hits = ins_bed.intersect(vntr_bed, wa=True).to_dataframe(
            names=["chrom", "start", "end", "minda_ID"]
        )
        if len(vntr_hits) > 0:
            vntr_ids = set(vntr_hits["minda_ID"].unique())

    # Assign summary categories
    classified_df["Summary_Category"] = classified_df.apply(
        lambda row: assign_summary_category(row, vntr_ids), axis=1
    )

    classified_df.to_csv(args.output, sep="\t", index=False)

    # Clean up pybedtools temp files
    pybedtools.cleanup()


if __name__ == "__main__":
    main()
