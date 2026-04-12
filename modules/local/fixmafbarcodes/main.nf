// Update MAF sample barcodes (BAMs lack RG headers so variant_utils cannot extract IDs)
process FIXMAFBARCODES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "mondrianscwgs/variant_calling:v0.1.4"

    input:
    tuple val(meta), path(maf, stageAs: "?/*")

    output:
    tuple val(meta), path("*.maf"), emit: maf
    tuple val("${task.process}"), val('fixmafbarcodes'), eval("python --version"), topic: versions, emit: versions_fixmafbarcodes

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fix_maf_barcodes.py \\
        ${meta.id} \\
        ${meta.id}_NORMAL \\
        ${maf} \\
        ${prefix}.maf
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.maf
    """
}
