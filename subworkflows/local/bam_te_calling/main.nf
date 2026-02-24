// find and annotate TE-mediated insertions
include { LONGCALLD         } from '../../../modules/local/longcalld/main'
include { BCFTOOLS_VIEW     } from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_ANNOTATE } from '../../../modules/nf-core/bcftools/annotate/main'
include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT     } from '../../../modules/nf-core/samtools/sort/main'
include { WHATSHAP_STATS    } from '../../../modules/local/whatshap/stats/main'
workflow BAM_TE_CALLING {
    take:
    bam_ch // channel: [ val(meta), bam, bai ]
    ref_fasta // val: reference fasta file
    longcalld_realign // boolean: output realigned cram files with longcalld

    main:
    ch_versions = channel.empty()
    // run longcallD
    LONGCALLD(
        bam_ch,
        [[id: "ref"], ref_fasta],
    )
    ch_versions = ch_versions.mix(LONGCALLD.out.versions.first())
    // get phasing stats from longcallD
    WHATSHAP_STATS(LONGCALLD.out.vcf)
    ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions.first())
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
    // index realigned bam
    if (longcalld_realign) {
        SAMTOOLS_SORT(
            LONGCALLD.out.cram,
            [[id: "ref"], ref_fasta],
            "crai",
        )
        channel.topic("versions").view()
    }

    emit:
    longcalld_vcf  = ch_longcalld_vcf // channel: [ val(meta), [ longcalld_vcf ] ]
    whatshap_stats = WHATSHAP_STATS.out.tsv // channel: [ val(meta), [path]]
    versions       = ch_versions // channel: [ val(meta), versions ]
}
