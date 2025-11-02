// somatic structural variant calling subworkflow
include { SEVERUS         } from '../../../modules/nf-core/severus/main'
include { SAVANA_CLASSIFY } from '../../../modules/local/savana/classify/main'
include { NANOMONSV_PARSE } from '../../../modules/local/nanomonsv/parse/main'
include { NANOMONSV_GET   } from '../../../modules/local/nanomonsv/get/main'
workflow SV_CALLING_SOMATIC {
    take:
    sv_callers   // val: list of sv callers to use
    input_sv_ch // channel: [ val(meta), tumor_bam, tumor_bai, norm_bam, norm_bai, vcf ]
    hap_bam_ch  // channel [(val), bam, bai]
    vntr_bed     // val: bed file of known VNTRs for severus
    ref_fasta    // val: reference fasta file
    ref_fai      // val: reference fasta index file

    main:

    ch_versions = Channel.empty()
    // run severus if specified
    if (sv_callers.split(',').contains('severus')) {
        SEVERUS(input_sv_ch, [[id: "ref"], vntr_bed])
        ch_versions = ch_versions.mix(SEVERUS.out.versions.first())
    }
    // run savana if specified
    if (sv_callers.split(',').contains('savana')) {
        SAVANA_CLASSIFY(
            input_sv_ch.map { meta, tumor_bam, tumor_bai, norm_bam, norm_bai, _vcf ->
                tuple(meta, tumor_bam, tumor_bai, norm_bam, norm_bai)
            },
            ref_fasta,
            ref_fai,
        )
        ch_versions = ch_versions.mix(SAVANA_CLASSIFY.out.versions.first())
    }
    // run nanomonsv if specified
    if (sv_callers.split(',').contains('nanomonsv')) {

        NANOMONSV_PARSE(hap_bam_ch)
        ch_versions = ch_versions.mix(NANOMONSV_PARSE.out.versions.first())
        // Combine all outputs into a single channel
        parse_out_ch = NANOMONSV_PARSE.out.parse_out
            .map { meta, parse_out ->
                tuple(meta.id, parse_out)
            }
            .groupTuple(by: 0)
        // now hand off to nanomonsv get
        input_get_ch = input_sv_ch
            .map { meta, tumor_bam, tumor_bai, norm_bam, norm_bai, _vcf ->
                tuple(meta.id, meta, tumor_bam, tumor_bai, norm_bam, norm_bai)
            }
            .join(parse_out_ch, by: 0)
            .map { _id, meta, tumor_bam, tumor_bai, norm_bam, norm_bai, parse_out ->
                tuple(meta, tumor_bam, tumor_bai, norm_bam, norm_bai, parse_out.flatten())
            }
        NANOMONSV_GET(input_get_ch, ref_fasta, ref_fai)
        ch_versions = ch_versions.mix(NANOMONSV_GET.out.versions.first())
    }

    emit:
    severus_vcf   = SEVERUS.out.somatic_vcf // channel: [ val(meta), [ somatic_vcf ] ]
    savana_vcf    = SAVANA_CLASSIFY.out.somatic_vcf // channel: [ val(meta), [ somatic_vcf ] ]
    nanomonsv_vcf = NANOMONSV_GET.out.somatic_vcf // channel: [ val(meta), [ somatic_vcf ] ]
    versions      = ch_versions // channel: [ versions.yml ]
}
