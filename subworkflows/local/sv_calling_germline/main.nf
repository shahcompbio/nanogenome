include { SEVERUS   } from '../../../modules/nf-core/severus/main'
include { CUTESV    } from '../../../modules/nf-core/cutesv/main'
include { SNIFFLES  } from '../../../modules/nf-core/sniffles/main'
include { LONGCALLD } from '../../../modules/local/longcalld/main'
include { PIGZ_UNCOMPRESS } from '../../../modules/nf-core/pigz/uncompress/main'

workflow SV_CALLING_GERMLINE {
    take:
    sv_callers     // val: list of sv callers to use
    ch_hap_bam     // channel: [ val(meta), [ haplotagged_bam ] ]
    ch_hap_bai     // channel: [ val(meta), [ haplotagged_bai ] ]
    rephased_vcf   // channel: [ val(meta), [ rephased_vcf ] ]
    vntr_bed       // val: bed file of known VNTRs for severus
    ref_fasta      // val: reference fasta file
    ch_samplesheet // channel: [ val(meta), bam]

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
    input_sv_ch = hap_bam_ch.norm
        .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
        .join(rephased_vcf.map { meta, vcf -> tuple(meta.id, meta, vcf) }, by: 0)
        .map { id, meta, norm_bam, norm_bai, _meta2, vcf ->
            tuple(meta, norm_bam, norm_bai, vcf)
        }
    // input_sv_ch.view()
    // run severus if specified
    if (sv_callers.split(',').contains('severus')) {
        SEVERUS(
            input_sv_ch.map { meta, norm_bam, norm_bai, vcf ->
                tuple(meta, norm_bam, norm_bai, [], [], vcf)
            },
            [[id: "ref"], vntr_bed],
        )
        ch_versions = ch_versions.mix(SEVERUS.out.versions.first())
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
    }
    // run longcallD
    if (sv_callers.split(',').contains('longcallD')) {
        // split ch_samplesheet into a tumor and normal channel
        ch_samplesheet
            .branch { meta, bam, bai ->
                tumor: meta.condition == 'tumor'
                norm: meta.condition == 'normal'
            }
            .set { bam_ch }
        LONGCALLD(bam_ch.norm,
        [[id: "ref"], ref_fasta])
        ch_versions = ch_versions.mix(LONGCALLD.out.versions.first())
    }

    emit:
    severus_vcf   = SEVERUS.out.all_vcf // channel: [ val(meta), [ germline_vcf ] ]
    cutesv_vcf    = CUTESV.out.vcf // channel: [ val(meta), [ germline_vcf ] ]
    sniffles_vcf  = PIGZ_UNCOMPRESS.out.file // channel: [ val(meta), [ germline_vcf ] ]
    longcalld_vcf = LONGCALLD.out.vcf // channel: [ val(meta), [ germline_vcf ] ]
    versions      = ch_versions // channel: [ versions.yml ]
}
