// longcalld for calling germline structural variants
process LONGCALLD {
    tag "${meta.id}"
    label 'process_high'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "biocontainers/longcalld:0.0.5--h7d57edc_0"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta1), path(fasta)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.vcf"), emit: vcf
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def contig_args = task.ext.contig_args ?: '--autosome-XY'
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    longcallD \\
        call -t16 \\
        ${fasta} \\
        ${bam} \\
        ${contig_args} \\
        ${args} \\
        > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longcallD: \$(longcallD --version)
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
        longcalld: \$(longcalld --version)
    END_VERSIONS
    """
}
