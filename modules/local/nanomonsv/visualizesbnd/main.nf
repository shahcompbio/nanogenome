// Visualize single-breakend SV contig annotations with BWA and RepeatMasker tracks
process NANOMONSV_VISUALIZESBND {
    tag "$meta.id"
    label 'process_low'

    // r-ggrepel 0.9.6 has no pre-built quay.io/biocontainers image.
    // Requires Wave (wave { enabled = true }) to build the container from environment.yml at runtime.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-ggrepel:0.9.6--r44hf9963bf_0':
        'quay.io/biocontainers/r-ggrepel:0.9.6--r44hf9963bf_0' }"

    input:
    tuple val(meta), path(bwa_txt), path(rmsk_txt)

    output:
    tuple val(meta), path("${meta.id}.nanomonsv.sbnd_vis"), emit: vis_dir
    tuple val("${task.process}"), val('r-ggrepel'), eval('Rscript -e "cat(as.character(packageVersion(\'ggrepel\')))"'), topic: versions, emit: versions_ggrepel

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript plot_sbnd_contigs.R \\
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
