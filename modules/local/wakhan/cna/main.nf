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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch solutions_ranks.tsv
    mkdir -p solution_1 coverage_plots
    touch ${prefix}_ploidy_purity.html
    touch ${prefix}_optimized_peak.html
    mkdir -p solution_1/subclonal_segments
    touch solution_1/subclonal_segments/${prefix}_subclonal_segments_HP_1.bed
    touch solution_1/subclonal_segments/${prefix}_subclonal_segments_HP_2.bed

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
}
