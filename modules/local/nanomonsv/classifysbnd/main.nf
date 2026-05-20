// Classify single-breakend SV contigs using BWA and RepeatMasker annotations
process NANOMONSV_CLASSIFYSBND {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanomonsv:0.9.0--pyhdfd78af_0':
        'quay.io/biocontainers/nanomonsv:0.9.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sbnd_result_txt), path(bwa_txt), path(rmsk_txt)

    output:
    tuple val(meta), path("${meta.id}.class.txt"), emit: class_txt
    tuple val("${task.process}"), val('nanomonsv'), eval('nanomonsv --version 2>&1 | sed "s/^nanomonsv //"'), topic: versions, emit: versions_nanomonsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    classify_sbnd_contigs.py \\
        ${sbnd_result_txt} \\
        ${bwa_txt} \\
        ${rmsk_txt} \\
        ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.class.txt
    """
}
