#!/usr/bin/env python
import sys
import argparse
import numpy as np
import pandas as pd
import os
from tqdm import tqdm

def check_oncokb(gsvs, oncokb):
    gene1 = list(gsvs["gene_name_1"])
    gene2 = list(gsvs["gene_name_2"])
    oncodat = pd.read_csv(oncokb, sep="\t")
    oncogenes = list(oncodat["Hugo Symbol"])
    gene1_oncogenes = [False] * len(gene1)
    gene2_oncogenes = gene1_oncogenes
    for i in np.arange(0, len(gene1)):
        g1 = gene1[i]
        g2 = gene2[i]
        if g1 in oncogenes:
            gene1_oncogenes[i] = True
        if g2 in oncogenes:
            gene2_oncogenes[i] = True
    gsvs["oncogene_gene1"] = gene1_oncogenes
    gsvs["oncogene_gene2"] = gene2_oncogenes
    return gsvs


def convert_sign(sign):
    mapping = {'+': 1, '-': -1}
    return mapping.get(sign, 0)  # Returns 0 if the sign is not "+" or "-"

def _fetch_gene_names(brk, refdat, window = 0):
    """
    fetch gene names
    :param brk: breakpoint tuple
    :param refdat: reference tsv file of gene positions in reference genome can create with R txdbmaker
    :param window: if you want to expand the search window for a gene name
    :return: hugo gene symbol and other overlapping hugo gene symbols at the position
    """
    chrom, pos = brk
    ### filter by chromosome
    df = refdat[refdat["chromosome_name"] == chrom]
    ## filter by position
    if window == 0:
        df = df[df["start_position"] < pos]
        df = df[df["end_position"] > pos]
    # expand window on either side
    else:
        df = df[df["start_position"] - (window/2) < pos]
        df = df[df["end_position"] + (window/2) > pos]
    ## generate gene names ...
    genes = list(df["hgnc_symbol"])
    if len(genes) == 1:
        gene_name = str(genes[0])
        alt_gene_names = 'nan'
    elif len(genes) > 1:
        gene_name = str(genes[0])
        alt_gene_names = genes[1:]
    else:
        gene_name = 'nan'
        alt_gene_names = 'nan'
    return gene_name, alt_gene_names

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("gene_annotation", help="gene annotation tsv")
    p.add_argument("svcalls", help="sv calls tsv")
    p.add_argument("oncokb", help="oncokb tsv")
    args = p.parse_args()
    gene_annotation = args.gene_annotation
    svcalls = args.svcalls
    oncokb = args.oncokb
    ## load SV calls
    svdat = pd.read_csv(svcalls, sep="\t")
    ## load gene annotation
    refdat = pd.read_csv(gene_annotation, sep="\t")
    ## iterate through svdat ...
    gene_names = {1:[], 2:[]}
    alt_gene_names = {1:[], 2:[]}
    svs = svdat

    ## set up progress bar ...
    pbar = tqdm(total=len(svs))
    for rix, row in svs.iterrows():
        if rix % 10 == 0: pbar.update(10)
        _, chrom1 = row["chrom1"].split("chr")
        pos1 = int(row["base1"])
        _, chrom2 = row["chrom2"].split("chr")
        pos2 = int(row["base2"])
        svtype = row["SV_Type"]
        # set breakpoints
        brks = [(chrom1, pos1), (chrom2, pos2)]
        if svtype == "INS":
            gene_name = 'nan'
            window = 0 # define precision of search window ...
            ## continue to expand the search window until we find a close match
            while gene_name == 'nan' or len(gene_name) < 1:
                gene_name, alt_gene_name = _fetch_gene_names(brks[0], refdat, window)
                ## stop looking and give up after expanding the search to 10kb
                if window == 10000:
                    gene_name = 'nan'
                    alt_gene_name = 'nan'
                    break
                window = window + 1000
            gene_names[1].append(gene_name)
            gene_names[2].append(gene_name)
            alt_gene_names[1].append(alt_gene_name)
            alt_gene_names[2].append(alt_gene_name)
        else:
            for i in np.arange(1, len(brks)+1):
                gene_name = 'nan'
                window = 0  # define precision of search window ...
                ## continue to expand the search window until we find a close match
                while gene_name == 'nan' or len(gene_name) < 1:
                    gene_name, alt_gene_name = _fetch_gene_names(brks[i-1], refdat, window)
                    ## stop looking and give up after expanding the search to 10kb
                    if window == 10000:
                        gene_name = 'nan'
                        alt_gene_name = 'nan'
                        break
                    window = window + 1000
                gene_names[i].append(gene_name)
                alt_gene_names[i].append(alt_gene_name)
    for ix in [1, 2]:
        svs[f'gene_name_{ix}'] = gene_names[ix]
        svs[f'alt_gene_name_{ix}'] = alt_gene_names[ix]


    print("check oncokb ....")
    svs = check_oncokb(svs, oncokb)
    svs.to_csv("annotated_sv_calls.tsv", sep='\t', index=False)
    # close progress bar
    pbar.close()
