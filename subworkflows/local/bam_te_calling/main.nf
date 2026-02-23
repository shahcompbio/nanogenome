// find and annotate TE-mediated insertions
include { LONGCALLD         } from '../../../modules/local/longcalld/main'
include { BCFTOOLS_VIEW     } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_ANNOTATE } from '../../../modules/nf-core/bcftools/annotate/main'
workflow BAM_TE_CALLING {
    take:
    bam_ch // channel: [ val(meta), bam, bai ]
    ref_fasta // val: reference fasta file

    main:
    ch_versions = channel.empty()
    // run longcallD
    LONGCALLD(
        bam_ch,
        [[id: "ref"], ref_fasta],
    )
    ch_versions = ch_versions.mix(LONGCALLD.out.versions.first())
    // filter longcallD calls for structural variants
    BCFTOOLS_VIEW(
        LONGCALLD.out.vcf.map { meta, vcf ->
            [meta, vcf, []]
        },
        [],
        [],
        [],
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())
    // give each SV an id so we can track
    BCFTOOLS_ANNOTATE(
        BCFTOOLS_VIEW.out.vcf.map { meta, vcf ->
            [meta, vcf, [], [], []]
        },
        [],
        [],
        [],
    )
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE.out.versions.first())
    ch_longcalld_vcf = BCFTOOLS_ANNOTATE.out.vcf

    emit:
    longcalld_vcf = ch_longcalld_vcf // channel: [ val(meta), [ longcalld_vcf ] ]
}
