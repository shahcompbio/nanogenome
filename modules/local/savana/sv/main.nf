// savana for somatic sv calling
process SAVANA_SV {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/savana", mode: 'copy', overwrite: true, saveAs: { filename -> "${meta.id}.${filename}" }

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/savana:1.3.6--pyhdfd78af_0"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(norm_bam), path(norm_bai)
    path ref_fasta

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("**/classified.somatic.vcf"), emit: somatic_vcf
    tuple val(meta), path("**/read_support.tsv"), emit: read_support
    tuple val(meta), path("**/inserted_sequences.fa"), emit: inserted_seqs
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def model_args = task.ext.model_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    savana \\
        ${args} \\
        ${model_args} \\
        --tumor ${tumor_bam} \\
        --normal ${norm_bam} \\
        --reference ${ref_fasta} \\
        --outdir ${prefix} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        savana: \$(savana --version)
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
        savana: \$(savana --version)
    END_VERSIONS
    """
}
