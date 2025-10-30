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
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    // TODO nf-core: If the module doesn't use arguments ($args), you SHOULD remove:
    //               - The definition of args `def args = task.ext.args ?: ''` above.
    //               - The use of the variable in the script `echo $args ` below.
    """
    echo ${args}

    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotcircos: \$(plotcircos --version)
    END_VERSIONS
    """
}
// get prefix for output circos
def getFilePrefix(filename) {
    def prefix = filename.split(/\.annotated_sv\.tsv/)[0]
    return prefix
}
