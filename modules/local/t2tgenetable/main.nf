// retrieve gene annotations for T2T-CHM13v2.0 genome using BiocT2T
process T2TGENETABLE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "shahlab_singularity/bioct2t:260412--4a244e4"

    output:
    path ("*-genes.txt"), emit: gene_annotation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    t2tgenetable.R

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        BiocT2T: \$(Rscript -e "library(BiocT2T); cat(as.character(packageVersion('BiocT2T')))")
        GenomicFeatures: \$(Rscript -e "library(GenomicFeatures); cat(as.character(packageVersion('GenomicFeatures')))")
    END_VERSIONS
    """

    stub:
    """
    touch t2t-genes.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        BiocT2T: stub
        GenomicFeatures: stub
    END_VERSIONS
    """
}
