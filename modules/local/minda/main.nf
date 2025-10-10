// merge SV calls from different callers
process MINDA {
    tag "${meta.id}"
    label 'process_low'
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "preskaa/minda:v250408"

    input:
    tuple val(meta), path(savana_vcf), path(severus_vcf), path(nanomonsv_vcf), val(min_support)
    val tolerance
    val min_size

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.minda_*.vcf"), emit: ensemble_vcf
    tuple val(meta), path("${meta.id}_min_callers_${min_support}"), emit: minda_out
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def set = min_support == 1 ? "union" : "consensus"
    """
    /minda/minda.py \\
        ensemble \\
        ${args} \\
        --vcfs \\
        ${savana_vcf} \\
        ${severus_vcf} \\
        ${nanomonsv_vcf} \\
        --sample_name ${meta.id} \\
        --min_support ${min_support} \\
        --tolerance ${tolerance} \\
        --min_size ${min_size} \\
        --out_dir ${prefix}_min_callers_${min_support}
    # rename file
    cp ${prefix}_min_callers_${min_support}/${prefix}_minda_ensemble.vcf ${prefix}.minda_${set}.vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minda: v250408
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def set = min_support == 1 ? "union" : "consensus"
    """
    touch ${prefix}.minda_${set}.vcf
    mkdir -p ${prefix}_min_callers_${min_support}
    touch ${prefix}_min_callers_${min_support}/${prefix}_minda_ensemble.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minda: v250408
    END_VERSIONS
    """
}
