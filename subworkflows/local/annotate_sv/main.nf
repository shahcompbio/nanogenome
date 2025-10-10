// merge and annotate structural variants from different callers
include { MINDA               } from '../../../modules/local/minda/main'
include { TABIX_BGZIPTABIX    } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'
include { ENSEMBLVEP_VEP      } from '../../../modules/nf-core/ensemblvep/vep/main'

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
    // download ensembl vep cache if needed
    ENSEMBLVEP_DOWNLOAD(
        [
            [id: "vep"],
            "GRCh38",
            "homo sapiens",
            115,
        ]
    )
    // annotate SVs
    ENSEMBLVEP_VEP(MINDA.out.ensemble_vcf.map{ meta, vcf -> tuple(meta, vcf, [])},
    "GRCh38",
    "homo_sapiens",
    115,
    ENSEMBLVEP_DOWNLOAD.out.cache.map {meta, cache -> cache},
    [[],[]],
    [])

    emit:
    minda_vcf = MINDA.out.ensemble_vcf // channel: [ val(meta), [ vcf ] ]
    versions  = ch_versions // channel: [ versions.yml ]
}
