#!/usr/bin/env python
import os
import re
import numpy as np
import pandas as pd
import argparse

def convert_hg19(df):
    chroms1 = list(df["chrom1"])
    chroms1 = ["chr"+chr for chr in chroms1]
    df["chrom1"] = chroms1
    chroms2 = list(df["chrom2"])
    chroms2 = ["chr"+chr for chr in chroms2]
    df["chrom2"] = chroms2
    return df



def proc_info(row):
    field = row['INFO'].split(';')
    infos = {}
    for each_info in field:
        if '=' in each_info:
            key, value = each_info.split('=')
            infos[key] = value  # infos['SVTYPE'] -> 'BND'
        else:
            infos[each_info] = True  # infos['IMPRECISE'] -> True
    return infos


def proc_samples(row):
    formats = row['FORMAT'].split(':')
    tumor_gts = row['TUMOR'].split(':')
    assert 'DR' in formats
    has_normal = False
    if 'NORMAL' in row:
        has_normal = True
        normal_gts = row['NORMAL'].split(':')
        return (has_normal,
                (dict(zip(formats, tumor_gts)),
                 dict(zip(formats, normal_gts))))
    return (has_normal,
            dict(zip(formats, tumor_gts)))

def infer_svtype(row, infos):
    """
    infer svtype from strand info or SVTYPE field

    :param infos: vcf info field
    """
    # if BND, infer from STRANDS
    if row['#CHROM'] != infos['CHR2']:
        sv_type = 'TRA'
    elif infos['SVTYPE'] == 'BND':
        assert 'STRANDS' in infos, f'infos does not have STRANDS:\n{row}'
        strands = infos['STRANDS']
        if strands == "+-":
            sv_type = "DEL"
        elif strands == "-+":
            sv_type = "DUP"
        elif strands == "--" or strands == "++":
            sv_type = "INV"
        else:
            sv_type = "BND"
    # else, use SVTYPE field
    else:
        sv_type = infos['SVTYPE']
    return sv_type


def get_svlen(infos):
        return int(infos['SVLEN'])


class MAF:
    def __init__(self, tumor_id, survivor=True):
        self.chrom1s = []
        self.pos1s = []
        self.chrom2s = []
        self.pos2s = []
        self.refs = []
        self.alts = []
        self.filters = []
        self.svtypes = []
        self.svlens = []
        self.callers = []
        self.minda_ID = []
        self.orientation = []
        self.info_keys = {'SVLEN', 'SVTYPE', 'SUPP_VEC'}
        self.tumor_id = tumor_id  # tumor id only
        self.has_normal = False  # does GT include normal

    def __repr__(self):
        return f'MAF of tumor_id {self.tumor_id}'

    def proc_vcf(self, vcf_path, survivor=False):
        if vcf_path.endswith('.gz'):
            vcf_file = gzip.open(vcf_path, 'rt', encoding='utf-8')
        else:
            vcf_file = open(vcf_path, 'r')
        for line in vcf_file:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line.rstrip().split('\t')
                assert len(header) <= 11, f'ERROR: header\n{header}\ntoo long'
                for ix in range(9, len(header)):
                    if header[ix] == self.tumor_id:  # tumor
                        header[ix] = 'TUMOR'
                    ## we can fix this when we move to somatic calls
                    # else:
                    #     header[ix] = 'NORMAL'
                self.header = header
                continue
            field = line.strip().split('\t')
            row = dict(zip(self.header, field))

            self.add_data(row, survivor=survivor)

    def add_data(self, row, survivor=False):
        infos = proc_info(row)
        remove_sv = False
        svlen = get_svlen(infos)
        svtype = infer_svtype(row, infos)
        assert len(self.info_keys & set(infos.keys())) == len(self.info_keys)
        if not remove_sv:
            self.chrom1s.append(row['#CHROM'])
            self.pos1s.append(int(row['POS']))
            self.chrom2s.append(infos["CHR2"])
            self.pos2s.append(infos["END"])
            self.minda_ID.append(row['ID'])
            self.refs.append(row['REF'])
            self.alts.append(row['ALT'])
            self.filters.append(row['FILTER'])
            self.orientation.append(infos['STRANDS'])
            self.svtypes.append(svtype)
            self.svlens.append(svlen)
            self.callers.append(infos["SUPP_VEC"])
            # self.rnames.append(rnames)
## for now, I'm giving this the same format as the fusions and nanomonsv
    def to_df(self):
        self.df = pd.DataFrame({
            'minda_ID': self.minda_ID,
            'chrom1': self.chrom1s,
            'base1': self.pos1s,
            'chrom2': self.chrom2s,
            'base2': self.pos2s,
            'Ref_Seq': self.refs,
            'Alt_Seq': self.alts,
            'SV_Type': self.svtypes,
            'SV_LEN': self.svlens,
            'orientation': self.orientation,
            'SV_callers': self.callers,
        })

     #   self.df['rnames'] = self.rnames

        return self.df


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("vcf")
    p.add_argument("--ref", default="hg38")
    args = p.parse_args()
    # make maf object
    mafobj = MAF('SAMPLE', survivor=False)
    vcf = args.vcf
    reference = args.ref
    mafobj.proc_vcf(vcf, survivor=False)
    maf = mafobj.to_df()
    if reference == "hg19":
        maf = convert_hg19(maf)
    prefix = os.path.splitext(os.path.basename(vcf))[0]
    tsv = f"{prefix}.tsv"
    print(f'writing {tsv}')
    maf.to_csv(tsv, sep='\t', index=False)
     # chunking for gene annotation
    print(f'chunking table for gene annotation')
    # make each chunk 20 rows or less
    chunk_size = 20
    num_splits = int(np.ceil(len(maf) / chunk_size))
    # Split the DataFrame
    split_dfs = np.array_split(maf, num_splits)
    # Convert split arrays back into DataFrames
    split_dfs = [pd.DataFrame(split) for split in split_dfs]
    for i in np.arange(0, len(split_dfs)):
        split_df = split_dfs[i]
        split_df.to_csv(f"chunk_{i}.tsv", sep='\t', index=False)
