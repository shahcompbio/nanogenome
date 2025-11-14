include { SEVERUS         } from '../../../modules/nf-core/severus/main'
include { CUTESV          } from '../../../modules/nf-core/cutesv/main'
include { SNIFFLES        } from '../../../modules/nf-core/sniffles/main'
include { LONGCALLD       } from '../../../modules/local/longcalld/main'
include { PIGZ_UNCOMPRESS } from '../../../modules/nf-core/pigz/uncompress/main'

workflow SV_CALLING_GERMLINE {
    take:
    sv_callers // val: list of sv callers to use
    input_sv_ch // channel: [ val(meta), norm_bam, norm_bai, vcf ]
    vntr_bed // val: bed file of known VNTRs for severus
    ref_fasta // val: reference fasta file

    main:

    ch_versions = Channel.empty()
    ch_severus_vcf = Channel.empty()
    ch_cutesv_vcf = Channel.empty()
    ch_sniffles_vcf = Channel.empty()
    ch_longcalld_vcf = Channel.empty()

    // run severus if specified
    if (sv_callers.split(',').contains('severus')) {
        SEVERUS(
            input_sv_ch.map { meta, norm_bam, norm_bai, vcf ->
                tuple(meta, norm_bam, norm_bai, [], [], vcf)
            },
            [[id: "ref"], vntr_bed],
        )
        ch_versions = ch_versions.mix(SEVERUS.out.versions.first())
        ch_severus_vcf = SEVERUS.out.all_vcf
    }
    // run cutesv if specified
    if (sv_callers.split(',').contains('cutesv')) {
        CUTESV(
            input_sv_ch.map { meta, norm_bam, norm_bai, _vcf ->
                tuple(meta, norm_bam, norm_bai)
            },
            [[id: "ref"], ref_fasta],
        )
        ch_versions = ch_versions.mix(CUTESV.out.versions.first())
        ch_cutesv_vcf = CUTESV.out.vcf
    }
    // run sniffles
    if (sv_callers.split(',').contains('sniffles')) {
        SNIFFLES(
            input_sv_ch.map { meta, norm_bam, norm_bai, _vcf ->
                tuple(meta, norm_bam, norm_bai)
            },
            [[id: "ref"], ref_fasta],
            [[], []],
            true,
            [],
        )
        ch_versions = ch_versions.mix(SNIFFLES.out.versions.first())
        // uncompress sniffles vcf for minda
        PIGZ_UNCOMPRESS(SNIFFLES.out.vcf)
        ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS.out.versions.first())
        ch_sniffles_vcf = PIGZ_UNCOMPRESS.out.file
    }
    // run longcallD
    if (sv_callers.split(',').contains('longcallD')) {
        longcalld_input_ch = input_sv_ch.map { meta, bam, bai, _snp_vcf ->
            tuple(meta, bam, bai)
        }
        LONGCALLD(
            longcalld_input_ch,
            [[id: "ref"], ref_fasta],
        )
        ch_versions = ch_versions.mix(LONGCALLD.out.versions.first())
        ch_longcalld_vcf = LONGCALLD.out.vcf
    }

    emit:
    severus_vcf   = ch_severus_vcf // channel: [ val(meta), [ germline_vcf ] ]
    cutesv_vcf    = ch_cutesv_vcf // channel: [ val(meta), [ germline_vcf ] ]
    sniffles_vcf  = ch_sniffles_vcf // channel: [ val(meta), [ germline_vcf ] ]
    longcalld_vcf = ch_longcalld_vcf // channel: [ val(meta), [ germline_vcf ] ]
    versions      = ch_versions // channel: [ versions.yml ]
}
