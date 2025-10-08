// merge and annotate structural variants from different callers
include { MINDA                      } from '../../../modules/local/minda/main'
include { TABIX_TABIX                } from '../../../modules/nf-core/tabix/tabix/main.nf'
include { ANNOTSV_INSTALLANNOTATIONS } from '../../../modules/nf-core/annotsv/installannotations/main'
include { ANNOTSV_ANNOTSV            } from '../../../modules/nf-core/annotsv/annotsv/main'

workflow ANNOTATE_SV {
    take:
    sv_ch     // channel: [ val(meta), path(vcf1), path(vcf2), path(vcf3), val(min_callers) ]
    tolerance // between breakpoints to consider two SVs the same
    min_size  // minimum size of SV to consider

    main:

    ch_versions = Channel.empty()
    // run minda to combine SV
    MINDA(sv_ch, params.tolerance, params.min_size)
    ch_versions = ch_versions.mix(MINDA.out.versions.first())
    // index vcfs
    TABIX_TABIX(MINDA.out.ensemble_vcf)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())
    // make a vcf/index channel for annotsv
    vcf_index_ch = MINDA.out.ensemble_vcf
        .join(TABIX_TABIX.out.tbi, by: 0)
        .map { meta, vcf, tbi -> tuple(meta, vcf, tbi, []) }
    // run annotsv
    ANNOTSV_INSTALLANNOTATIONS()
    ANNOTSV_ANNOTSV(
        vcf_index_ch,
        ANNOTSV_INSTALLANNOTATIONS.out.annotations,
        [],
        [],
        [],
    )

    ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions.first())

    emit:
    minda_vcf     = MINDA.out.ensemble_vcf // channel: [ val(meta), [ vcf ] ]
    annotsv_table = ANNOTSV_ANNOTSV.out.tsv // channel: [ val(meta), [ vcf ] ]
    versions      = ch_versions // channel: [ versions.yml ]
}
