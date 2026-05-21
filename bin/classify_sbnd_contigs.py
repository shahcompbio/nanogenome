#!/usr/bin/env python3
"""
Classify single-breakend SV contigs using BWA and RepeatMasker annotation results.
Usage: python classify_sbnd_contigs.py <sbnd_result_txt> <bwa_txt> <rmsk_txt> <output_prefix>
Outputs: {output_prefix}.class.txt
Adapted from nanomonsv/misc/subscript_sbnd/annotate_contig.py
"""
import csv, sys
from collections import namedtuple

REPEAT_RATIO_THRES = 0.8
L1_NAMES = {"L1HS", "L1P1", "L1PA2"}
MIN_ALIGN_LEN = 2000
MIN_MQ = 40
EARLY_THRES = 100

Rmsk = namedtuple("Rmsk", "QA_Start QA_End QA_Strand RName RClass RA_Start RA_End")
Bwa  = namedtuple("Bwa",  "QA_Start QA_End QA_Strand TA_Chr TA_Start TA_End MQ")

def revcomp(seq):
    comp = {"A":"T","C":"G","G":"C","T":"A","W":"W","S":"S","M":"K","K":"M","R":"Y","Y":"R","B":"V","V":"B","D":"H","H":"D","N":"N"}
    return "".join(comp.get(b, b) for b in reversed(seq))

class Contig:
    def __init__(self, cid, seq):
        parts = cid.split(",")
        self.bp_chr, self.bp_pos, self.bp_strand = parts[0], int(parts[1]), parts[2]
        self.contig, self.contig_len = seq, len(seq)
        self.contig_class = "None"
        self.rmsk_info, self.bwa_info = [], []
        self.chr1, self.pos1, self.dir1 = self.bp_chr, self.bp_pos, self.bp_strand
        self.chr2 = self.pos2 = self.dir2 = self.inseq = None

    def sv_key(self):
        if self.chr2 is None: return "---"
        if self.inseq == "": self.inseq = "---"
        if self.chr1 > self.chr2 or (self.chr1 == self.chr2 and self.pos1 > self.pos2):
            self.chr1, self.chr2 = self.chr2, self.chr1
            self.pos1, self.pos2 = self.pos2, self.pos1
            self.dir1, self.dir2 = self.dir2, self.dir1
            self.inseq = revcomp(self.inseq)
        return f"{self.chr1},{self.pos1},{self.dir1},{self.chr2},{self.pos2},{self.dir2},{self.inseq}"

    def classify(self):
        self._simple_satellite()
        if self.contig_class == "None": self._plain_sv()
        if self.contig_class == "None": self._l1_mediated_del()
        if self.contig_class == "None": self.contig_class = "Complex"

    def _simple_satellite(self):
        sizes = {"Simple_repeat": 0, "Satellite": 0, "Satellite/centr": 0}
        for r in self.rmsk_info:
            if r.RClass in sizes: sizes[r.RClass] += r.QA_End - r.QA_Start
        for cls, size in sizes.items():
            if size / self.contig_len > REPEAT_RATIO_THRES:
                self.contig_class = cls; return

    def _plain_sv(self):
        for r in self.bwa_info:
            if r.QA_Start < EARLY_THRES and r.QA_End - r.QA_Start >= MIN_ALIGN_LEN and r.MQ >= MIN_MQ:
                self.contig_class = "Plain_SV"
                self.chr2 = r.TA_Chr
                self.pos2, self.dir2 = (r.TA_Start, "-") if r.QA_Strand == "+" else (r.TA_End, "+")
                self.inseq = self.contig[: r.QA_Start - 1]; return

    def _l1_mediated_del(self):
        self.bwa_info = sorted(self.bwa_info, key=lambda x: x.QA_Start)
        l1_segs = [(r.QA_Start, r.QA_End) for r in self.rmsk_info if r.RName in L1_NAMES]
        l1_inter, qsize = 0, 0
        for i in range(min(len(self.bwa_info) - 1, 2)):
            r = self.bwa_info[i]
            qsize += r.QA_End - r.QA_Start
            for ls in l1_segs:
                if ls[1] >= r.QA_Start and ls[0] <= r.QA_End:
                    l1_inter += min(ls[1], r.QA_End) - max(ls[0], r.QA_Start)
            if float(l1_inter) / max(qsize, 1) <= REPEAT_RATIO_THRES: continue
            nxt = self.bwa_info[i + 1]
            nxt_len = nxt.QA_End - nxt.QA_Start
            if nxt_len >= MIN_ALIGN_LEN and nxt.MQ >= MIN_MQ:
                self._assign_l1(nxt); return
            if nxt_len >= MIN_ALIGN_LEN and nxt.TA_Chr == self.bp_chr:
                if (nxt.QA_Strand == self.bp_strand == "+" and -100 < nxt.TA_Start - self.bp_pos < 10000) or \
                   (nxt.QA_Strand == self.bp_strand == "-" and -10000 < nxt.TA_End - self.bp_pos < 100):
                    self._assign_l1(nxt); return

    def _assign_l1(self, r):
        self.contig_class = "L1_Mediated_Del"
        self.chr2 = r.TA_Chr
        self.pos2, self.dir2 = (r.TA_Start, "-") if r.QA_Strand == "+" else (r.TA_End, "+")
        self.inseq = self.contig[: r.QA_Start - 1]

if __name__ == "__main__":
    sbnd_result_txt, bwa_txt, rmsk_txt, output_prefix = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

    contigs = {}
    with open(sbnd_result_txt) as fh:
        for F in csv.DictReader(fh, delimiter="\t"):
            cid = f'{F["Chr_1"]},{F["Pos_1"]},{F["Dir_1"]},{F["SV_ID"]}'
            contigs[cid] = Contig(cid, F["Contig"])

    with open(rmsk_txt) as fh:
        for F in csv.DictReader(fh, delimiter="\t"):
            cid = F["Contig_ID"]
            if cid not in contigs: continue
            contigs[cid].rmsk_info.append(Rmsk(int(F["Query_Align_Start"]), int(F["Query_Align_End"]), F["Query_Align_Strand"], F["Repeat_NAME"], F["Repeat_Class"], int(F["Repeat_Align_Start"]), int(F["Repeat_Align_End"])))

    with open(bwa_txt) as fh:
        for F in csv.DictReader(fh, delimiter="\t"):
            cid = F["Contig_ID"]
            if cid not in contigs: continue
            contigs[cid].bwa_info.append(Bwa(int(F["Query_Align_Start"]), int(F["Query_Align_End"]), F["Query_Align_Strand"], F["Target_Align_Chromosome"], int(F["Target_Align_Start"]), int(F["Target_Align_End"]), int(F["Mapping_Quality"])))

    with open(output_prefix + ".class.txt", "w") as hout:
        print("Contig_ID\tContig_Class\tSV_Key", file=hout)
        for cid, c in contigs.items():
            c.classify()
            print(f"{cid}\t{c.contig_class}\t{c.sv_key()}", file=hout)
