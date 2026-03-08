#!/usr/bin/env python
import sys

tumor_id = sys.argv[1]
normal_id = sys.argv[2]
infile = sys.argv[3]
outfile = sys.argv[4]
t_idx = n_idx = None
with open(infile) as fin, open(outfile, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        if fields[0] == "Hugo_Symbol":
            try:
                t_idx = fields.index("Tumor_Sample_Barcode")
                n_idx = fields.index("Matched_Norm_Sample_Barcode")
            except ValueError:
                pass
            fout.write(line)
        else:
            if t_idx is not None:
                fields[t_idx] = tumor_id
            if n_idx is not None:
                fields[n_idx] = normal_id
            fout.write("\t".join(fields) + "\n")
