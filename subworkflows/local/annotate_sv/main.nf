// merge and annotate structural variants from different callers
include { MINDA            } from '../../../modules/local/minda/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { VCF2TSV          } from '../../../modules/local/vcf2tsv/main'
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
    TABIX_BGZIPTABIX(MINDA.out.ensemble_vcf)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())
    // convert to tsv
    VCF2TSV(MINDA.out.ensemble_vcf)
    ch_versions = ch_versions.mix(VCF2TSV.out.versions.first())

    emit:
    minda_vcf = MINDA.out.ensemble_vcf // channel: [ val(meta), [ vcf ] ]
    sv_table  = VCF2TSV.out.sv_table // channel: [ val(meta), path(sv_table) ]
    versions  = ch_versions // channel: [ versions.yml ]
}
