//  phase-correct and generates rephased snp VCFs
process WAKHAN_HAPCORRECT {
    tag "${meta1.id}"
    label 'process_high'
    stageInMode 'copy'
    publishDir "${params.outdir}/wakhan/${meta1.id}", mode: 'copy', overwrite: false
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/wakhan:94effdd"

    input:
    tuple val(meta), path(ref_fasta)
    tuple val(meta1), path(bam), path(phased_vcf)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta1), path("hapcorrect_out/*", arity: '3..*'), emit: wakhanHPOutput
    tuple val(meta1), path("**/rephased.vcf.gz"), emit: rephased_vcf, optional: true
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tabix ${phased_vcf}
    wakhan \\
        hapcorrect \\
        ${args} \\
        ${args1} \\
        ${args2} \\
        --threads ${task.cpus} \\
        --reference ${ref_fasta}  \\
        --target-bam ${bam} \\
        --normal-phased-vcf ${phased_vcf} \\
        --genome-name ${meta1.id} \\
        --out-dir-plots hapcorrect_out

    WAKHAN_VERSION=\$(python3 -c "
    import sys
    sys.path.insert(0, '/opt/wakhan/Wakhan')
    from src.__version__ import __version__
    print(__version__)
    ")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wakhan: \$WAKHAN_VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wakhan: \$(wakhan --version)
    END_VERSIONS
    """
}
