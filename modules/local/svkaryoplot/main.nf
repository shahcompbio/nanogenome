// karyoplot to show germline indel distribution
process SVKARYOPLOT {
    tag "${meta.id}"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "shahlab_singularity/r-tidyverse-karyoploter:4.5.1--1.36.0"

    input:
    tuple val(meta), path(sv_table)
    val genome_build

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.svg"), emit: karyoplot
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def set = meta.min_callers == 1 ? "union" : "consensus"
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.condition}.${set}"
    """
    svkaryoplot.R \\
        ${sv_table} \\
        ${prefix}.karyoplot.svg \\
        ${genome_build} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        karyoploteR: \$(Rscript -e "library(circlize); cat(as.character(packageVersion('karyoploteR')))")
        BSgenome: \$(Rscript -e "library(circlize); cat(as.character(packageVersion('BSgenome')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def set = meta.min_callers == 1 ? "union" : "consensus"
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.condition}.${set}"
    """
    echo ${args}

    touch ${prefix}.karyoplot.svg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        karyoploteR: \$(Rscript -e "library(circlize); cat(as.character(packageVersion('karyoploteR')))")
        BSgenome: \$(Rscript -e "library(circlize); cat(as.character(packageVersion('BSgenome')))")
        dplyr: \$(Rscript -e "library(dplyr); cat(as.character(packageVersion('dplyr')))")
    END_VERSIONS
    """
}
