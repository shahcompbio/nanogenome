#!/bin/bash
## specify params
outdir=/data1/shahs3/users/preskaa/SarcAtlas/data/nanogenome/test
pipelinedir=$HOME/nanogenome
mkdir -p ${outdir}
cd ${outdir}

nextflow run ${pipelinedir}/main.nf \
    -profile test,singularity \
    -work-dir ${outdir}/work \
    --outdir ${outdir} \
    --gene_annotations /data1/shahs3/reference/ref-sarcoma/GRCh38/v45/Homo_sapiens.GRCh38.111.annotations-genes.txt \
    --oncokb /data1/shahs3/reference/ref-sarcoma/241115_oncokb_cancerGeneList.tsv \
    -resume
