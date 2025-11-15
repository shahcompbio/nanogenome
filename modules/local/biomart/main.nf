// retrieve gene annotations from biomart
process BIOMART {
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "biocontainers/bioconductor-biomart:2.62.0--r44hdfd78af_0"

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    path ("*-genes.txt"), emit: gene_annotation
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // default genome build is hg38
    def args = task.ext.args ?: 'hg38'
    """
    biomart.R ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        biomaRt: \$(Rscript -e "library(biomaRt); cat(as.character(packageVersion('biomaRt')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: 'hg38'
    """
    touch ${args}-genes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        biomaRt: \$(Rscript -e "library(biomaRt); cat(as.character(packageVersion('biomaRt')))")
    END_VERSIONS
    """
}
