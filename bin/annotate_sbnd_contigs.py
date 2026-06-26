#!/usr/bin/env python3
"""
Annotate single-breakend SV contig sequences with BWA mem and RepeatMasker.
Usage: python annotate_sbnd_contigs.py <sbnd_result_txt> <bwa_index_prefix> <output_prefix>
Adapted from nanomonsv/misc/subscript_sbnd/annotate_contig.py
"""

import csv, os, shutil, subprocess, sys
import pysam

# Set CSV field size limit to maximum safe value (handles OverflowError on
# platforms where C long is smaller than sys.maxsize)
maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt / 10)


def proc_rmsk(input_file, output_file):
    with open(input_file) as hin, open(output_file, "w") as hout:
        print(
            "Contig_ID\tContig_Len\tQuery_Align_Start\tQuery_Align_End\tQuery_Align_Strand\tRepeat_NAME\tRepeat_Class\tRepeat_Align_Start\tRepeat_Align_End",
            file=hout,
        )
        for line in hin:
            F = line.rstrip("\n").split()
            if len(F) <= 1 or F[0] in ["SW", "score", "There"]:
                continue
            contigs = F[4].split(",")
            contig_id, contig_len = ",".join(contigs[:-1]), contigs[-1]
            if F[8] == "C":
                F[8] = "-"
            ra_start, ra_end = (F[11], F[12]) if F[8] == "+" else (F[13], F[12])
            print(
                f"{contig_id}\t{contig_len}\t{F[5]}\t{F[6]}\t{F[8]}\t{F[9]}\t{F[10]}\t{ra_start}\t{ra_end}",
                file=hout,
            )


def proc_sam(input_sam, output_file):
    with open(output_file, "w") as hout:
        print(
            "Contig_ID\tContig_Len\tQuery_Align_Start\tQuery_Align_End\tQuery_Align_Strand\tTarget_Align_Chromosome\tTarget_Align_Start\tTarget_Align_End\tMapping_Quality",
            file=hout,
        )
        samfile = pysam.AlignmentFile(input_sam, "r")
        for read in samfile.fetch():
            if read.is_unmapped or read.is_secondary:
                continue
            query_strand = "-" if read.is_reverse else "+"
            query_length = read.infer_read_length()
            cigartuples = read.cigartuples
            left_clip = cigartuples[0][1] if cigartuples[0][0] == 5 else 0
            right_clip = cigartuples[-1][1] if cigartuples[-1][0] == 5 else 0
            if not read.is_supplementary:
                if query_strand == "+":
                    query_start, query_end = (
                        read.query_alignment_start + 1,
                        read.query_alignment_end,
                    )
                else:
                    query_start, query_end = (
                        query_length - read.query_alignment_end + 1,
                        query_length - read.query_alignment_start,
                    )
            else:
                if query_strand == "+":
                    query_start, query_end = left_clip + 1, query_length - right_clip
                else:
                    query_start, query_end = right_clip + 1, query_length - left_clip
            contigs = read.query_name.split(",")
            contig_id, contig_len = ",".join(contigs[:-1]), contigs[-1]
            print(
                f"{contig_id}\t{contig_len}\t{query_start}\t{query_end}\t{query_strand}\t{read.reference_name}\t{read.reference_start + 1}\t{read.reference_end}\t{read.mapping_quality}",
                file=hout,
            )
        samfile.close()


if __name__ == "__main__":
    sbnd_result_txt, bwa_index_prefix, output_prefix = (
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
    )
    tmp_fasta = output_prefix + ".tmp.fasta"
    tmp_rmsk_dir = output_prefix + ".tmp.rmsk"
    tmp_bwa_sam = output_prefix + ".tmp.bwa.sam"

    with open(sbnd_result_txt) as hin, open(tmp_fasta, "w") as hout:
        for F in csv.DictReader(hin, delimiter="\t"):
            print(
                f'>{F["Chr_1"]},{F["Pos_1"]},{F["Dir_1"]},{F["SV_ID"]},{len(F["Contig"])}\n{F["Contig"]}',
                file=hout,
            )

    os.makedirs(tmp_rmsk_dir, exist_ok=True)
    subprocess.check_call(
        ["RepeatMasker", "-species", "human", tmp_fasta, "-dir", tmp_rmsk_dir]
    )
    boutput_prefix = os.path.basename(output_prefix)
    proc_rmsk(
        os.path.join(tmp_rmsk_dir, boutput_prefix + ".tmp.fasta.out"),
        output_prefix + ".nanomonsv.rmsk.txt",
    )

    with open(tmp_bwa_sam, "w") as hout:
        subprocess.check_call(
            ["bwa", "mem", "-h", "200", bwa_index_prefix, tmp_fasta], stdout=hout
        )
    proc_sam(tmp_bwa_sam, output_prefix + ".nanomonsv.bwa.txt")

    shutil.rmtree(tmp_rmsk_dir)
    os.remove(tmp_fasta)
    os.remove(tmp_bwa_sam)
