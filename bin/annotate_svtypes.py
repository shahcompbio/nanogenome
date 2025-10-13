#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np

# paths
caller_tables = sys.argv[1]
input_svtable = sys.argv[2]
out_svtable = sys.argv[3]

def get_bp_ids(call_ids):
    """
    get breakpoint ids from individual callers
    :param call_ids: bp ids from minda output
    :return: list of bp ids
    """
    bp_ids = []
    for call in call_ids:
        terms = call.split("_")
        bp_id = "_".join(terms[1:])
        bp_ids.append(bp_id)
    return bp_ids

# load in the input from the different individual callers as one unified pandas dataframe
all_callers = pd.DataFrame()
# make list from call tables
calls_list = caller_tables.split()
for i in np.arange(0, len(calls_list)):
    table_path = calls_list[i]
    caller_table = pd.read_csv(table_path, sep="\t")
    caller_table["caller"] = "caller_{i}".format(i=i)
    # uniform column names for merging
    # rename columns
    caller_table.columns = ["chrom", "pos", "ID", "SV_Type",
                           "Strands", "Strand", "BP_notation", "caller"]
    all_callers = pd.concat([all_callers, caller_table])
# fix the strands by merging in BP_notation column (some callers use this instead of Strands)
all_callers["Strands"] = all_callers.apply(
    lambda row: row["Strands"] if row["Strands"] != "." else row["BP_notation"],
    axis=1
)
# do the same by also merging in the Strand column
all_callers["Strands"] = all_callers.apply(
    lambda row: row["Strands"] if row["Strands"] != "." else row["Strand"],
    axis=1
)
# Drop the original BP_notation column if no longer needed
all_callers = all_callers.drop(columns=["BP_notation", "Strand"])
# load input table of consensus SV calls
all_svtable = pd.read_csv(input_svtable, sep="\t")
# search for BNDs (the unresolved breakpoint)
bnd_table = all_svtable[all_svtable["SV_Type"] == "BND"]
# reannotate translocations first because they're easy
tra_table = bnd_table[bnd_table["chrom1"] != bnd_table["chrom2"]]
tra_table["SV_Type"] = "TRA"
# now let's re-annotate unresolved breakpoints to their closest match
# based on strands and annotation from individual callers
unresolved_table = bnd_table[bnd_table["chrom1"] == bnd_table["chrom2"]]
sv_types = []
for i, row in unresolved_table.iterrows():
    call_ids = row["SV_callers"].split(",")
    # get breakpoint ids
    bp_ids = get_bp_ids(call_ids)
    # search for these ids in the all_callers table
    bp_df = all_callers[all_callers["ID"].isin(bp_ids)]
    # iterate through until we can annotate
    sv_type = "BND"
    j = 0
    while sv_type == "BND" and j < len(bp_df):
        bp_row = bp_df.iloc[j]
        strands = bp_row["Strands"]
        if bp_row["SV_Type"] != "BND":
            sv_type = bp_row["SV_Type"]
        elif strands == "++" or strands == "--":
            sv_type = "INV"
        elif strands == "+-":
            sv_type = "DEL"
        elif strands == "-+":
            sv_type = "DUP"
        j = j + 1
    # append sv type
    sv_types.append(sv_type)

# make new table
reannotated_table = unresolved_table
reannotated_table["SV_Type"] = sv_types

# stitch everything together into final output
resolved_table = all_svtable[all_svtable["SV_Type"] != "BND"]
# add translocations
final_table = pd.concat([resolved_table, tra_table, reannotated_table])
# get strands for final table to make compatible with recon plot package
strand_table = all_callers[all_callers["Strands"]!= "."]
strand_list = []
for _, row in final_table.iterrows():
    if row["SV_Type"] == "INS":
        strand = "INS"
    elif row["SV_Type"] == "DEL":
        strand = "+-"
    else:
        call_ids = row["SV_callers"].split(",")
        # get breakpoint ids
        bp_ids = get_bp_ids(call_ids)
        bp_df = strand_table[strand_table["ID"].isin(bp_ids)]
        strand = "."
        i = 0
        while strand == "." and i < len(bp_df):
            bp_row = bp_df.iloc[i]
            strand = bp_row["Strands"]
            i = i + 1
    strand_list.append(strand)
final_table["strands"] = strand_list
# rearrange table to make pretty
col_order = ["minda_ID",
             "chrom1",
             "base1",
             "chrom2",
             "base2",
             "strands",
             "Ref_Seq",
             "Alt_Seq",
             "SV_Type",
             "SV_LEN",
             "SV_callers",
             "gene_name_1",
             "alt_gene_name_1",
             "gene_name_2",
             "alt_gene_name_2",
             "oncokb_gene1",
             "oncokb_gene2"]
# write output
final_table.to_csv(out_svtable, sep="\t", index=False)
