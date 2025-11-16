// whatshap to haplotag reads
process WHATSHAP_HAPLOTAG {
    tag "${meta.id}"
    label 'process_high'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/whatshap:2.8--py39h2de1943_0'
        : 'biocontainers/whatshap:2.8--py39h2de1943_0'}"

    input:
    path ref_fasta
    path ref_fai
    tuple val(meta), path(bam), path(bai), path(vcf), path(vcf_tbi)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.haplotagged.bam"), emit: bam
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.condition}"
    """
    whatshap haplotag \\
        ${args} \\
        --reference ${ref_fasta} \\
        ${vcf} \\
        ${bam} \\
        -o ${prefix}.haplotagged.bam \\
        --ignore-read-groups \\
        --tag-supplementary \\
        --skip-missing-contigs \\
        --output-threads ${task.cpus}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.condition}"
    """
    touch ${prefix}.haplotagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        whatshap: \$(whatshap --version)
    END_VERSIONS
    """
}
