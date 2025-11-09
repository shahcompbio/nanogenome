#!/bin/bash
## specify params
outdir=/data1/shahs3/users/preskaa/SarcAtlas/data/nanogenome/test
pipelinedir=$HOME/nanogenome
mkdir -p ${outdir}
cd ${outdir}

nextflow run ${pipelinedir}/main.nf \
    -profile test_cna_only,singularity \
    -work-dir ${outdir}/work \
    --outdir ${outdir} \
    --ascat_allele_files /data1/shahs3/reference/ref-sarcoma/GRCh38/ascat/G1000_alleles_hg38.zip \
    --ascat_loci_files /data1/shahs3/reference/ref-sarcoma/GRCh38/ascat/G1000_loci_hg38.zip \
    -resume
