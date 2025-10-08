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
    --annotsv_annotations /data1/papaemme/isabl/ref/homo_sapiens/GRCh37d5/germline_svs/3.4.4/
    -resume
