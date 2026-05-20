// Merge per-contig SBND visualization PDFs into a single combined PDF
process NANOMONSV_MERGESBNDPDFS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppler:25.07.0':
        'quay.io/biocontainers/poppler:25.07.0' }"

    input:
    tuple val(meta), path(vis_dir)

    output:
    tuple val(meta), path("${meta.id}.nanomonsv.sbnd.pdf"), emit: merged_pdf
    tuple val("${task.process}"), val('pdfunite'), eval('pdfunite --version 2>&1 | head -1'), topic: versions, emit: versions_pdfunite

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pdfunite \\
        ${vis_dir}/*.pdf \\
        ${prefix}.nanomonsv.sbnd.pdf \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nanomonsv.sbnd.pdf
    """
}
