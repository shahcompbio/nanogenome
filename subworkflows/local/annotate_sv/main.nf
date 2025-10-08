// merge and annotate structural variants from different callers
include { MINDA                      } from '../../../modules/local/minda/main'
include { TABIX_BGZIPTABIX           } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { ANNOTSV_INSTALLANNOTATIONS } from '../../../modules/nf-core/annotsv/installannotations/main'
include { ANNOTSV_ANNOTSV            } from '../../../modules/nf-core/annotsv/annotsv/main'

workflow ANNOTATE_SV {
    take:
    sv_ch               // channel: [ val(meta), path(vcf1), path(vcf2), path(vcf3), val(min_callers) ]
    tolerance           // between breakpoints to consider two SVs the same
    min_size            // minimum size of SV to consider
    annotsv_annotations // path to annotsv annotations if already installed

    main:

    ch_versions = Channel.empty()
    // run minda to combine SV
    MINDA(sv_ch, params.tolerance, params.min_size)
    ch_versions = ch_versions.mix(MINDA.out.versions.first())
    // index vcfs
    TABIX_BGZIPTABIX(MINDA.out.ensemble_vcf)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())
    // run annotsv
    if (!annotsv_annotations) {
        ANNOTSV_INSTALLANNOTATIONS()
        ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions.first())
        annotation_ch = ANNOTSV_INSTALLANNOTATIONS.out.annotations
    }
    else {
        annotation_ch = Channel.value(annotsv_annotations)
    }


    ANNOTSV_ANNOTSV(
        TABIX_BGZIPTABIX.out.gz_tbi.map { meta, vcf, tbi ->
            tuple(meta, vcf, tbi, [])
        },
        annotation_ch,
        [],
        [],
        [],
    )
    ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions.first())

    emit:
    minda_vcf     = MINDA.out.ensemble_vcf // channel: [ val(meta), [ vcf ] ]
    annotsv_table = ANNOTSV_ANNOTSV.out.tsv // channel: [ val(meta), [ vcf ] ]
    versions      = ch_versions // channel: [ versions.yml ]
}
