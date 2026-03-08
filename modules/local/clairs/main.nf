process CLAIRS {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://docker.io/hkubal/clairs:v0.4.4'
        : 'docker.io/hkubal/clairs:v0.4.4'}"

    input:
    tuple val(meta), path(norm_bam), path(norm_bai), path(tumor_bam), path(tumor_bai)
    path ref_fasta
    path ref_fai

    output:
    tuple val(meta), path("**/snv.vcf.gz"), emit: snv_vcf
    tuple val(meta), path("**/snv.vcf.gz.tbi"), emit: snv_vcf_tbi
    tuple val(meta), path("**/indel.vcf.gz"), emit: indel_vcf
    tuple val(meta), path("**/indel.vcf.gz.tbi"), emit: indel_vcf_tbi
    tuple val("${task.process}"), val('clairs'), eval("run_clairs --version 2>&1 | sed -n 's/.*\\(v[0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_clairs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_clairs \\
        --tumor_bam_fn ${tumor_bam} \\
        --normal_bam_fn ${norm_bam} \\
        --ref_fn ${ref_fasta} \\
        --threads ${task.cpus} \\
        --output_dir ${prefix} \\
        --enable_indel_calling \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}
    mkdir -p ${prefix}
    touch ${prefix}/snv.vcf.gz
    touch ${prefix}/snv.vcf.gz.tbi
    touch ${prefix}/indel.vcf.gz
    touch ${prefix}/indel.vcf.gz.tbi
    """
}
