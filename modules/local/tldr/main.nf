// tldr for transposable element detection from long-read data
process TLDR {
    tag "${meta.id}"
    label 'process_high'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "shahlab_singularity/tldr:260224--7c6dfda"

    input:
    tuple val(meta), path(bams, arity: '1..*'), path(bais, arity: '1..*')
    path te_ref_fasta
    path ref_genome
    path ref_genome_fai

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.table.txt"), emit: tsv
    tuple val("${task.process}"), val('tldr'), eval("tldr --version"), topic: versions, emit: versions_tldr

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bams_list = bams.join(',')
    """
    tldr --bams ${bams_list} \\
        --procs ${task.cpus} \\
        --elts ${te_ref_fasta} \\
        --ref ${ref_genome} \\
        --trdcol \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.table.txt
    """
}
