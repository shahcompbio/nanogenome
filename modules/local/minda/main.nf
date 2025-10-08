// merge SV calls from different callers
process MINDA {
    tag "${meta.id}"
    label 'process_low'
    // rename vcf for union vs consensus
    publishDir "minda", mode: 'copy', overwrite: true, saveAs: { filename ->
        if (filename.endsWith("_minda_ensemble.vcf")) {
            min_support == 1 ? "${meta.id}_union.vcf" : "${meta.id}_consensus.vcf"
        }
        else {
            filename
        }
    }

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "preskaa/minda:v250408"

    input:
    tuple val(meta), path(savana_vcf), path(severus_vcf), path(nanomonsv_vcf), val(min_support)
    val tolerance
    val min_size
    path filter_bed

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("**/*minda_ensemble.vcf"), emit: ensemble_vcf
    tuple val(meta), path("${meta.id}_min_callers_${min_support}"), emit: minda_out
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    /minda/minda.py \\
        ensemble \\
        ${args} \\
        ${savana_vcf} \\
        ${severus_vcf} \\
        ${nanomonsv_vcf} \\
        --sample_name ${meta.id} \\
        --min_support ${min_support} \\
        --tolerance ${tolerance} \\
        --min_size ${min_size} \\
        --bed ${filter_bed}
        --out_dir ${prefix}_min_callers_${min_support}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minda: \$(minda --version)
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
        minda: \$(minda --version)
    END_VERSIONS
    """
}
