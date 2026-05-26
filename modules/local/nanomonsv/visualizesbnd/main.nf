// Visualize single-breakend SV contig annotations with BWA and RepeatMasker tracks
process NANOMONSV_VISUALIZESBND {
    tag "${meta.id}"
    label 'process_low'

    // Container built with Wave from environment.yml and pushed to quay.io/shahlab_singularity.
    conda "${moduleDir}/environment.yml"
    container 'quay.io/shahlab_singularity/nanomonsv-visualizesbnd:r-tidyverse-2.0.0_r-ggrepel-0.9.6--293b10b22ffcfa3f'

    input:
    tuple val(meta), path(bwa_txt), path(rmsk_txt)

    output:
    tuple val(meta), path("${meta.id}.nanomonsv.sbnd_vis"), emit: vis_dir
    tuple val("${task.process}"), val('r-ggrepel'), eval('Rscript -e "cat(as.character(packageVersion(\'ggrepel\')))"'), topic: versions, emit: versions_ggrepel

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_sbnd_contigs.R \\
        ${prefix} \\
        ${prefix}.nanomonsv.sbnd_vis \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.nanomonsv.sbnd_vis
    """
}
