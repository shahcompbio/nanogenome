// import modules
include { CLAIR3            } from '../../../modules/nf-core/clair3/main'
include { LONGPHASE_PHASE   } from '../../../modules/nf-core/longphase/phase/main'
include { WAKHAN_HAPCORRECT } from '../../../modules/local/wakhan/hapcorrect/main'
include { TABIX_TABIX       } from '../../../modules/nf-core/tabix/tabix/main'
include { WHATSHAP_HAPLOTAG } from '../../../modules/local/whatshap/haplotag/main'
include { SAMTOOLS_INDEX    } from '../../../modules/nf-core/samtools/index/main'
include { WHATSHAP_STATS } from '../../../modules/local/whatshap/stats/main'
/*
 * haplotag bam files
 */
workflow HAPLOTAG {
    take:
    ch_samplesheet  // channel: [ val(meta), [ bam ] ]
    clair3_model    // val: clair3 model specification
    clair3_platform // val: clair3 platform specification
    fasta           // val: reference fasta
    fai             // val: reference fasta fai

    main:

    ch_versions = Channel.empty()
    // split ch_samplesheet into a tumor and normal channel
    ch_samplesheet
        .branch { meta, bam, bai ->
            tumor: meta.condition == 'tumor'
            norm: meta.condition == 'normal'
        }
        .set { bam_ch }
    // run clair3 for germline snps
    clair_input_ch = bam_ch.norm.map { meta, bam, bai ->
        tuple(meta, bam, bai, clair3_model, [], clair3_platform)
    }
    CLAIR3(clair_input_ch, [[id: "ref"], fasta], [[id: "ref"], fai])
    ch_versions = ch_versions.mix(CLAIR3.out.versions)
    // run longphase to phase SNPs
    longphase_input_ch = bam_ch.norm
        .join(CLAIR3.out.vcf, by: 0)
        .map { meta, bam, bai, vcf ->
            tuple(meta, bam, bai, vcf, [], [])
        }
    LONGPHASE_PHASE(longphase_input_ch, [[id: "ref"], fasta], [[id: "ref"], fai])
    ch_versions = ch_versions.mix(LONGPHASE_PHASE.out.versions)
    // phase correct tumor bam using phased SNPs
    hapcorrect_input_ch = bam_ch.tumor
        .map { meta, bam, bai -> tuple(meta.id, meta, bam) }
        .join(LONGPHASE_PHASE.out.vcf.map { meta, vcf -> tuple(meta.id, meta, vcf) }, by: 0)
        .map { id, tumor_meta, bam, norm_meta, vcf -> tuple(tumor_meta, bam, vcf) }
    WAKHAN_HAPCORRECT([[id: "ref"], fasta], hapcorrect_input_ch)
    ch_versions = ch_versions.mix(WAKHAN_HAPCORRECT.out.versions)
    // tabix rephased vcf if it exists
    rephased_vcf_ch = WAKHAN_HAPCORRECT.out.rephased_vcf
        .concat(LONGPHASE_PHASE.out.vcf)
        .first()
    rephased_vcf_ch.view()
    TABIX_TABIX(rephased_vcf_ch)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)
    // compute phasing statistics
    WHATSHAP_STATS(rephased_vcf_ch)
    ch_versions = ch_versions.mix(WHATSHAP_STATS.out.versions)
    // run whatshap haplotag to tag both tumor and normal bams
    hap_input_ch = ch_samplesheet
        .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
        .combine(
            rephased_vcf_ch.map { meta, vcf -> tuple(meta.id, meta, vcf) },
            by: 0
        )
        .combine(
            TABIX_TABIX.out.tbi.map { meta, tbi -> tuple(meta.id, meta, tbi) },
            by: 0
        )
        .map { id, bam_meta, bam, bai, vcf_meta, vcf, tbi_meta, tbi ->
            tuple(bam_meta, bam, bai, vcf, tbi)
        }
    WHATSHAP_HAPLOTAG(params.fasta, fai, hap_input_ch)
    ch_versions = ch_versions.mix(WHATSHAP_HAPLOTAG.out.versions)
    // index haplotagged bams
    SAMTOOLS_INDEX(WHATSHAP_HAPLOTAG.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    bam            = WHATSHAP_HAPLOTAG.out.bam // channel: [ val(meta), [ bam ] ]
    bai            = SAMTOOLS_INDEX.out.bai // channel: [ val(meta), [ bai ] ]
    rephased_vcf   = rephased_vcf_ch // channel: [ val(meta), [ vcf ] ]
    bam_snps       = hap_input_ch // channel: [ val(meta), bam, bai, vcf, tbi ] ]
    wakhanHPOutput = WAKHAN_HAPCORRECT.out.wakhanHPOutput // channel: [ val(meta), [ path ] ]
    whatshap_stats = WHATSHAP_STATS.out.tsv // channel: [ val(meta), [path]]
    versions       = ch_versions // channel: [ versions.yml ]
}
