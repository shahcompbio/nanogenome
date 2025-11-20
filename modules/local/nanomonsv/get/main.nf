// collect nanomonsv sv results
process NANOMONSV_GET {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/nanomonsv/${meta.id}", mode: 'copy', overwrite: true

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/nanomonsv:0.8.0--pyhdfd78af_0'
        : 'biocontainers/nanomonsv:0.8.0--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(norm_bam), path(norm_bai), path(nanomonsv_parse)
    path ref_fasta
    path ref_fai

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.nanomonsv.result.txt"), emit: result_table
    tuple val(meta), path("*.nanomonsv.result.vcf"), emit: somatic_vcf
    tuple val(meta), path("*.sbnd.result.txt"), emit: sbnd_table
    tuple val(meta), path("*supporting_read.txt"), emit: read_support
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    nanomonsv get \\
        ${meta.id}.tumor \\
        ${tumor_bam} \\
        ${ref_fasta} \\
        --control_prefix ${meta.id}.normal \\
        --control_bam ${norm_bam} \\
        --processes ${task.cpus} \\
        --single_bnd \\
        --use_racon \\
        --max_memory_minimap2 8 \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanomonsv: \$(echo \$(nanomonsv --version 2>&1) | sed 's/^nanomonsv //')
        mafft: \$(echo \$(mafft --version 2>&1) | sed 's/^v//; s/ (.*//')
        racon: \$(echo \$(racon --version 2>&1) | sed 's/^v//')
        tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^tabix (htslib) //; s/ Copyright.*//')
        bgzip: \$(echo \$(bgzip --version 2>&1) | sed 's/^bgzip (htslib) //; s/ Copyright.*//')
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${meta.id}.tumor.nanomonsv.result.txt
    touch ${meta.id}.tumor.nanomonsv.result.vcf
    touch ${meta.id}.tumor.sbnd.result.txt
    touch ${meta.id}.tumor.nanomonsv.supporting_read.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanomonsv: \$(echo \$(nanomonsv --version 2>&1) | sed 's/^nanomonsv //')
        mafft: \$(echo \$(mafft --version 2>&1) | sed 's/^v//; s/ (.*//')
        racon: \$(echo \$(racon --version 2>&1) | sed 's/^v//')
        tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^tabix (htslib) //; s/ Copyright.*//')
        bgzip: \$(echo \$(bgzip --version 2>&1) | sed 's/^bgzip (htslib) //; s/ Copyright.*//')
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
