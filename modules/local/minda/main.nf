// merge SV calls from different callers
process MINDA {
    tag "${meta.id}"
    label 'process_low'
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "shahlab_singularity/minda:20260108b-8b6d81c"

    input:
    tuple val(meta), path(vcfs, arity: "2..*", stageAs: "?/*")
    val tolerance
    val min_size

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*_minda_*.vcf"), emit: ensemble_vcf
    tuple val(meta), path("*_min_callers_${meta.min_callers}"), emit: minda_out
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def set = meta.min_callers == 1 ? "union" : "consensus"
    def min_support = meta.min_callers
    """
    /minda/minda.py \\
        ensemble \\
        ${args} \\
        --vcfs \\
        ${vcfs} \\
        --sample_name ${meta.id} \\
        --min_support ${min_support} \\
        --tolerance ${tolerance} \\
        --min_size ${min_size} \\
        --out_dir ${prefix}_min_callers_${min_support}
    # rename file
    cp ${prefix}_min_callers_${min_support}/${meta.id}_minda_ensemble.vcf ${prefix}_minda_${set}.vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minda: v250408
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_support = meta.min_callers
    def set = meta.min_callers == 1 ? "union" : "consensus"
    """
    touch ${prefix}_minda_${set}.vcf
    mkdir -p ${prefix}_min_callers_${min_support}
    touch ${prefix}_min_callers_${min_support}/${meta.id}_minda_ensemble.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minda: v250408
    END_VERSIONS
    """
}
