// rephasing and copy number analysis with wakhan (full workflow)
process WAKHAN_REPHASE_CNA {
    tag "${meta.id}"
    label 'process_high'
    stageInMode 'copy'
    publishDir "${params.outdir}/wakhan/${meta.id}", mode: 'copy', overwrite: true, saveAs: { filename -> filename.startsWith("solution_") ? "cna_solutions/${filename}" : filename }

    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/wakhan:260128-crankycrank-bc48900"

    input:
    tuple val(meta), path(bam), path(bai), path(phased_vcf), path(phased_vcf_tbi), path(severus_vcf)
    path ref_fasta
    path ref_fai

    output:
    // Rephasing outputs
    tuple val(meta), path("**/rephased.vcf.gz"), emit: rephased_vcf, optional: true
    // CNA outputs
    tuple val(meta), path("solutions_ranks.tsv"), emit: solutions_ranks
    tuple val(meta), path("solution_*", arity: '1..*'), emit: wakhanCNAOutput
    tuple val(meta), path("coverage_plots"), emit: coverage_plots
    tuple val(meta), path("*_ploidy_purity.html"), emit: ploidy_purity_html
    tuple val(meta), path("*_optimized_peak.html"), emit: optimized_peak_html
    tuple val(meta), path("solution_1/**/*_subclonal_segments_HP_1.bed"), emit: HP1_bed
    tuple val(meta), path("solution_1/**/*_subclonal_segments_HP_2.bed"), emit: HP2_bed

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args1 = task.ext.args1 ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    wakhan \\
        ${args} \\
        ${args1} \\
        ${args2} \\
        ${args3} \\
        --threads ${task.cpus} \\
        --reference ${ref_fasta} \\
        --target-bam ${bam} \\
        --normal-phased-vcf ${phased_vcf} \\
        --genome-name ${meta.id} \\
        --breakpoints ${severus_vcf} \\
        --out-dir-plots .

    WAKHAN_VERSION=\$(python3 -c "
    import sys
    sys.path.insert(0, '/opt/Wakhan')
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

    touch rephased.vcf.gz
    touch rephased.vcf.gz.tbi
    touch solutions_ranks.tsv
    mkdir -p solution_1/subclonal_segments coverage_plots
    touch ${prefix}_ploidy_purity.html
    touch ${prefix}_optimized_peak.html
    touch solution_1/subclonal_segments/${prefix}_subclonal_segments_HP_1.bed
    touch solution_1/subclonal_segments/${prefix}_subclonal_segments_HP_2.bed

    WAKHAN_VERSION=\$(python3 -c "
    import sys
    sys.path.insert(0, '/opt/Wakhan')
    from src.__version__ import __version__
    print(__version__)
    ")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wakhan: \$WAKHAN_VERSION
    END_VERSIONS
    """
}
