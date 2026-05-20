// Annotate single-breakend SV contig sequences with BWA mem and RepeatMasker
process NANOMONSV_ANNOTATESBND {
    tag "$meta.id"
    label 'process_high'
    stageInMode 'copy'

    // No pre-built container exists for this tool combination (nanomonsv + bwa + repeatmasker).
    // Requires Wave (wave { enabled = true }) to build the container from environment.yml at runtime.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-nanomonsv-bwa-repeatmasker:0.9.0':
        'quay.io/biocontainers/mulled-v2-nanomonsv-bwa-repeatmasker:0.9.0' }"

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
