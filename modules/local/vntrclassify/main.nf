// Classify unclassified insertions using VNTR intersection and priority hierarchy
process VNTRCLASSIFY {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'docker://quay.io/shahlab_singularity/vntrclassify:pybedtools-0.10.0_bedtools-2.31.1'
        : 'quay.io/shahlab_singularity/vntrclassify:pybedtools-0.10.0_bedtools-2.31.1'}"

    input:
    tuple val(meta), path(classified_tsv)
    path vntr_bed

    output:
    tuple val(meta), path("${meta.id}_somatic_inserts.final_classified.tsv"), emit: final_classified_tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vntr_classify.py \\
        --classified_tsv ${classified_tsv} \\
        --vntr_bed ${vntr_bed} \\
        --output ${prefix}_somatic_inserts.final_classified.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pybedtools: \$(python -c "import pybedtools; print(pybedtools.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_somatic_inserts.final_classified.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pybedtools: \$(python -c "import pybedtools; print(pybedtools.__version__)")
    END_VERSIONS
    """
}
