// copy number analysis with wakhan
process WAKHAN_CNA {
    tag "${meta.id}"
    label 'process_high'
    stageInMode 'copy'
    publishDir "${params.outdir}/wakhan/${meta.id}", mode: 'copy', overwrite: true, saveAs: { filename -> filename.startsWith("solution_") ? "cna_solutions/${filename}" : filename }

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/wakhan:94effdd"

    input:
    tuple val(meta), path(bam), path(bai), path(phased_vcf), path(phased_vcf_tbi), path(severus_vcf), path(wakhanHPOutput)
    path ref_fasta

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("solutions_ranks.tsv"), emit: solutions_ranks
    tuple val(meta), path("solution_*", arity: '1..*'), emit: wakhanCNAOutput
    tuple val(meta), path("coverage_plots"), emit: coverage_plots
    tuple val(meta), path("*_ploidy_purity.html"), emit: ploidy_purity_html
    tuple val(meta), path("*_optimized_peak.html"), emit: optimized_peak_html
    tuple val(meta), path("solution_1/**/*_subclonal_segments_HP_1.bed"), emit: HP1_bed
    tuple val(meta), path("solution_1/**/*_subclonal_segments_HP_2.bed"), emit: HP2_bed

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
    wakhan cna \\
        ${args} \\
        ${args1} \\
        ${args2} \\
        --threads ${task.cpus} \\
        --reference ${ref_fasta} \\
        --target-bam ${bam} \\
        --normal-phased-vcf ${phased_vcf} \\
        --genome-name ${meta.id} \\
        --breakpoints ${severus_vcf} \\
        --use-sv-haplotypes \\
        --out-dir-plots .

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
        wakhan: \$(wakhan --version)
    END_VERSIONS
    """
}
