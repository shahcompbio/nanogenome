// annotate sv types
process SVTYPES {
    tag "$meta.id"
    label 'process_single'
    publishDir "${params.outdir}/annotated_sv", mode: 'copy', overwrite: true

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "preskaa/annotate_genes:v240817"

    input:
    tuple val(meta), path(caller_tables, stageAs: "?/*"), path(merged_table)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.tsv"), emit: annotated_sv
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def set = meta.min_callers == 1 ? "union" : "consensus"
    """
    annotate_svtypes.py \\
        \"${caller_tables}\" \\
        ${merged_table} \\
        ${prefix}.${set}.annotated_sv.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def set = meta.min_callers == 1 ? "union" : "consensus"
    """
    touch ${prefix}.${set}.annotated_sv.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}
