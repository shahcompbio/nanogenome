// Extract insertions from annotated SV table and resolve symbolic <INS> sequences
process PREPINSERTTABLE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://quay.io/shahlab_singularity/prepinserttable:pandas-2.2.0_python-3.12'
        : 'quay.io/shahlab_singularity/prepinserttable:pandas-2.2.0_python-3.12'}"

    input:
    tuple val(meta), path(annotated_sv_tsv), path(nanomonsv_result), path(severus_vcf)

    output:
    tuple val(meta), path("${meta.id}_somatic_inserts.tsv"), emit: inserts_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    prep_insert_table.py \\
        --annotated_sv ${annotated_sv_tsv} \\
        --nanomonsv_result ${nanomonsv_result} \\
        --severus_vcf ${severus_vcf} \\
        --output ${prefix}_somatic_inserts.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_somatic_inserts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
