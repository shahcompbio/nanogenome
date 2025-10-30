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
    --germline false \
    --skip_somatic \
    -resume
