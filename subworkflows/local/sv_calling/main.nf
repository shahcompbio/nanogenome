// structural variant calling subworkflow
include { SEVERUS } from '../../../modules/nf-core/severus/main'

workflow SV_CALLING {
    take:
    sv_callers   // val: list of sv callers to use
    ch_hap_bam   // channel: [ val(meta), [ haplotagged_bam ] ]
    ch_hap_bai   // channel: [ val(meta), [ haplotagged_bai ] ]
    rephased_vcf // channel: [ val(meta), [ rephased_vcf ] ]
    vntr_bed     // val: bed file of known VNTRs for severus

    main:

    ch_versions = Channel.empty()
    // construct input channel for SV callers
    // branch tumor and normal
    hap_bam_ch = ch_hap_bam
        .join(ch_hap_bai, by: 0)
        .branch { meta, bam, bai ->
            tumor: meta.condition == 'tumor'
            norm: meta.condition == 'normal'
        }
    // sv caller input channel
    input_sv_ch = hap_bam_ch.tumor
        .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
        .join(hap_bam_ch.norm.map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }, by: 0)
        .join(rephased_vcf.map { meta, vcf -> tuple(meta.id, meta, vcf) }, by: 0)
        .map { id, tumor_meta, tumor_bam, tumor_bai, norm_meta, norm_bam, norm_bai, meta3, vcf ->
            tuple([id: id], tumor_bam, tumor_bai, norm_bam, norm_bai, vcf)
        }
    // run severus if specified
    if (sv_callers.split(',').contains('severus')) {
        SEVERUS(input_sv_ch, [[id: "ref"], vntr_bed])
        ch_versions = ch_versions.mix(SEVERUS.out.versions)
    }

    emit:
    severus_vcf = SEVERUS.out.somatic_vcf // channel: [ val(meta), [ somatic_vcf ] ]
    versions    = ch_versions // channel: [ versions.yml ]
}
