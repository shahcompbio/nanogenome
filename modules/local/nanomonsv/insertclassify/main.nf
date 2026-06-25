// Classify insertion sequences using nanomonsv insert_classify
process NANOMONSV_INSERTCLASSIFY {
    tag "${meta.id}"
    label 'process_medium'
    stageInMode 'copy'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/nanomonsv:0.9.0--pyhdfd78af_0'
        : 'biocontainers/nanomonsv:0.9.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(inserts_tsv)
    path ref_fasta
    path bwa_index, stageAs: 'bwa_index' // BWA index directory; copy-staged to avoid concurrent access issues
    path ref_gtf
    path line1_db
    path line1_db_tbi // tabix index staged alongside line1_db

    output:
    tuple val(meta), path("${meta.id}_somatic_inserts.classified.tsv"), emit: classified_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ln -sf bwa_index/${ref_fasta.baseName}.* .
    nanomonsv insert_classify \\
        ${inserts_tsv} \\
        ${prefix}_somatic_inserts.classified.tsv \\
        ${ref_fasta} \\
        ${ref_gtf} \\
        ${line1_db} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanomonsv: \$(echo \$(nanomonsv --version 2>&1) | sed 's/^nanomonsv //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_somatic_inserts.classified.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanomonsv: \$(echo \$(nanomonsv --version 2>&1) | sed 's/^nanomonsv //')
    END_VERSIONS
    """
}
