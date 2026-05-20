// Annotate single-breakend SV contig sequences with BWA mem and RepeatMasker
process NANOMONSV_ANNOTATESBND {
    tag "$meta.id"
    label 'process_high'
    stageInMode 'copy'

    // Container built with Wave from environment.yml and pushed to quay.io/shahlab_singularity.
    conda "${moduleDir}/environment.yml"
    container 'quay.io/shahlab_singularity/nanomonsv-annotatesbnd:nanomonsv-0.9.0_bwa-0.7.18_repeatmasker-4.1.7_pysam-0.22.1--e4220aa301415d2a'

    input:
    tuple val(meta), path(sbnd_result_txt)
    path ref_fasta
    path ref_fai
    path bwa_index

    output:
    tuple val(meta), path("${meta.id}.nanomonsv.bwa.txt"), path("${meta.id}.nanomonsv.rmsk.txt"), emit: annotations
    tuple val("${task.process}"), val('bwa'),          eval('bwa 2>&1 | grep -m1 Version | sed "s/Version: //"'),            topic: versions, emit: versions_bwa
    tuple val("${task.process}"), val('repeatmasker'), eval('RepeatMasker -v 2>&1 | head -1 | sed "s/RepeatMasker version //"'), topic: versions, emit: versions_repeatmasker

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args ?: ''
    def prefix           = task.ext.prefix ?: "${meta.id}"
    def bwa_index_prefix = "${bwa_index}/${ref_fasta.baseName}"
    """
    annotate_sbnd_contigs.py \\
        ${sbnd_result_txt} \\
        ${bwa_index_prefix} \\
        ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nanomonsv.bwa.txt
    touch ${prefix}.nanomonsv.rmsk.txt
    """
}
