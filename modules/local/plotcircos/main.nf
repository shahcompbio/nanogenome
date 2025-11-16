// plot circos of SVs and CNAs
process PLOTCIRCOS {
    tag "${meta.id}"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "shahlab_singularity/r-tidyverse-circlize:4.5.1--0.4.16"

    input:
    tuple val(meta), path(sv_table), path(cn_hp1_bed), path(cn_hp2_bed)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.svg"), emit: circos
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: getFilePrefix("${sv_table}")
    """
    plotcircos.R \\
        ${sv_table} \\
        ${prefix}.circos.svg \\
        ${cn_hp1_bed} \\
        ${cn_hp2_bed} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        circlize: \$(Rscript -e "library(circlize); cat(as.character(packageVersion('circlize')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: getFilePrefix("${sv_table}")
    """
    touch ${prefix}.circos.svg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        circlize: \$(Rscript -e "library(circlize); cat(as.character(packageVersion('circlize')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
// get prefix for output circos
def getFilePrefix(filename) {
    def prefix = filename.split(/\.annotated_sv\.tsv/)[0]
    return prefix
}
