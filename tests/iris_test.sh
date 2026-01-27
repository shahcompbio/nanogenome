#!/bin/bash
## specify params
outdir=/data1/shahs3/users/preskaa/SarcAtlas/data/nanogenome/sv_test
pipelinedir=$HOME/nanogenome
mkdir -p ${outdir}
cd ${outdir}
annotsv_dir=/data1/shahs3/reference/ref-sarcoma/AnnotSV_annotations/

nextflow run ${pipelinedir}/main.nf \
    -profile test_sv_only,singularity \
    -work-dir ${outdir}/work \
    --outdir ${outdir} \
    --annotsv_dir ${annotsv_dir} \
    -resume
