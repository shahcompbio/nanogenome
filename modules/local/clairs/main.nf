process CLAIRS {
    tag "${meta.id}"
    label 'process_high'

    container "docker://hkubal/clairs:v0.4.4"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(norm_bam), path(norm_bai)
    path ref_fasta
    path ref_fai
    val platform

    output:
    tuple val(meta), path("output/output.vcf.gz"),     emit: vcf
    tuple val(meta), path("output/output.vcf.gz.tbi"), emit: tbi
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    run_clairs \\
        --tumor_bam_fn ${tumor_bam} \\
        --normal_bam_fn ${norm_bam} \\
        --ref_fn ${ref_fasta} \\
        --threads ${task.cpus} \\
        --platform ${platform} \\
        --output_dir output \\
        --enable_indel_calling \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairs: \$(run_clairs --version 2>&1 | grep -oP 'v[0-9.]+')
    END_VERSIONS
    """
}
