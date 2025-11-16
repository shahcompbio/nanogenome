// longcalld for calling germline structural variants
process LONGCALLD {
    tag "${meta.id}"
    label 'process_high'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "biocontainers/longcalld:0.0.6--h7d57edc_0"

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
        call -t ${task.cpus} \\
        ${fasta} \\
        ${bam} \\
        ${contig_args} \\
        ${args} \\
        > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longcallD: \$(longcallD call -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longcallD: \$(longcallD call -v)
    END_VERSIONS
    """
}
