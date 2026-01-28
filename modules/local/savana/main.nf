// classify putative structural variants as somatic or germline
process SAVANA {
    tag "${meta.id}"
    label 'process_high'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "quay.io/biocontainers/savana:1.3.6--pyhdfd78af_0"

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(norm_bam), path(norm_bai), path(snp_vcf)
    path ref_fasta
    path ref_fai
    path contigs_list

    output:
    // sv calling outputs
    tuple val(meta), path("**/*.classified.somatic.vcf"), emit: somatic_vcf
    tuple val(meta), path("**/*.sv_breakpoints_read_support.tsv"), emit: read_support
    tuple val(meta), path("**/*.inserted_sequences.fa"), emit: inserted_seqs
    // cna outpus
    tuple val(meta), path("**/*_raw_read_counts.tsv"), emit: read_counts, optional: true
    tuple val(meta), path("**/*_read_counts_mnorm_log2r_segmented.tsv"), emit: log2r, optional: true
    tuple val(meta), path("**/*_fitted_purity_ploidy.tsv"), emit: purity_ploidy, optional: true
    tuple val(meta), path("**/*_segmented_absolute_copy_number.tsv"), emit: cnv, optional: true

    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def model_args = task.ext.model_args ?: ''
    // arguments for copy number analysis
    def cna_args = snp_vcf ? "--snp_vcf ${snp_vcf} --cna_threads ${task.cpus}" : ''
    def contigs_arg = contigs_list ? "--contigs ${contigs_list}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    savana \\
        -t ${tumor_bam} \\
        -n ${norm_bam} \\
        --ref ${ref_fasta} \\
        --outdir ${prefix} \\
        --threads ${task.cpus} \\
        ${args} \\
        ${model_args} \\
        ${cna_args} \\
        ${contigs_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        savana: \$(MPLBACKEND=Agg MPLCONFIGDIR=/tmp savana --version 2>&1 | tail -n 1 | sed 's/SAVANA //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.classified.somatic.vcf
    touch ${prefix}/${prefix}.sv_breakpoints_read_support.tsv
    touch ${prefix}/${prefix}.inserted_sequences.fa
    touch ${prefix}/${prefix}_raw_read_counts.tsv
    touch ${prefix}/${prefix}_read_counts_mnorm_log2r_segmented.tsv
    touch ${prefix}/${prefix}_fitted_purity_ploidy.tsv
    touch ${prefix}/${prefix}_segmented_absolute_copy_number.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        savana: \$(MPLBACKEND=Agg MPLCONFIGDIR=/tmp savana --version 2>&1 | tail -n 1 | sed 's/SAVANA //')
    END_VERSIONS
    """
}
