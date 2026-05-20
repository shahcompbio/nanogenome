// Subworkflow: Classify and visualize single-breakend SVs from nanomonsv
include { NANOMONSV_ANNOTATESBND  } from '../../../modules/local/nanomonsv/annotatesbnd/main'
include { NANOMONSV_CLASSIFYSBND  } from '../../../modules/local/nanomonsv/classifysbnd/main'
include { NANOMONSV_VISUALIZESBND } from '../../../modules/local/nanomonsv/visualizesbnd/main'
include { NANOMONSV_MERGESBNDPDFS } from '../../../modules/local/nanomonsv/mergesbndpdfs/main'

workflow SBND_CLASSIFY {
    take:
    sbnd_result_ch // channel: [ val(meta), path(sbnd_result_txt) ]
    ref_fasta // path: reference genome FASTA
    ref_fai // path: reference genome FAI index
    bwa_index // path: pre-built BWA index directory (--bwa_index param)

    main:
    // Step 1: Align contigs with BWA, annotate with RepeatMasker
    NANOMONSV_ANNOTATESBND(
        sbnd_result_ch,
        ref_fasta,
        ref_fai,
        bwa_index,
    )

    // Step 2: Classify contigs using annotation results
    // Join sbnd_result_txt back to annotations for the classifier
    classify_input_ch = sbnd_result_ch
        .map { meta, sbnd -> [meta.id, meta, sbnd] }
        .join(
            NANOMONSV_ANNOTATESBND.out.annotations.map { meta, bwa, rmsk -> [meta.id, bwa, rmsk] },
            by: 0
        )
        .map { _id, meta, sbnd, bwa, rmsk -> [meta, sbnd, bwa, rmsk] }

    NANOMONSV_CLASSIFYSBND(classify_input_ch)

    // Step 3: Visualize and merge PDFs (default on; skipped with --skip_sbnd_vis)
    if (!params.skip_sbnd_vis) {
        NANOMONSV_VISUALIZESBND(NANOMONSV_ANNOTATESBND.out.annotations)
        NANOMONSV_MERGESBNDPDFS(NANOMONSV_VISUALIZESBND.out.vis_dir)
    }

    emit:
    sbnd_classes = NANOMONSV_CLASSIFYSBND.out.class_txt // channel: [ val(meta), path(class_txt) ]
    sbnd_vis     = params.skip_sbnd_vis ? Channel.empty() : NANOMONSV_VISUALIZESBND.out.vis_dir // channel: [ val(meta), path(vis_dir) ]
    sbnd_pdf     = params.skip_sbnd_vis ? Channel.empty() : NANOMONSV_MERGESBNDPDFS.out.merged_pdf // channel: [ val(meta), path(merged_pdf) ]
}
