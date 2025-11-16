// somatic copy number analysis
include { WAKHAN_REPHASE_CNA   } from '../../../modules/local/wakhan/rephase_cna/main'
include { SAVANA as SAVANA_CNA } from '../../../modules/local/savana/main'


workflow BAM_CNV_CALLING_SOMATIC {
    take:
    tools // val (cnv_tools)
    cna_input_ch // channel: [ val(meta), [ tumor_bam, tumor_bai, norm_bam, norm_bai, snp_vcf, snp_tbi, sv_vcf]]
    fasta // path (fasta)
    fai // path (fai)
    contigs_list // path: contigs list file for savana (optional)

    main:

    ch_versions = channel.empty()
    ch_hp1_bed = channel.empty()
    ch_hp2_bed = channel.empty()
    ch_savana_cnv = channel.empty()

    // run wakahn for cna
    if (tools.split(',').contains('wakhan')) {
        wakhan_input_ch = cna_input_ch.map { meta, tumor_bam, tumor_bai, _norm_bam, _norm_bai, snp_vcf, snp_tbi, sv_vcf ->
            tuple(meta, tumor_bam, tumor_bai, snp_vcf, snp_tbi, sv_vcf)
        }
        WAKHAN_REPHASE_CNA(wakhan_input_ch, fasta, fai)
        ch_versions = ch_versions.mix(WAKHAN_REPHASE_CNA.out.versions.first())
        ch_hp1_bed = WAKHAN_REPHASE_CNA.out.HP1_bed
        ch_hp2_bed = WAKHAN_REPHASE_CNA.out.HP2_bed
    }
    // run savana for cna
    if (tools.split(',').contains('savana')) {
        savana_input_ch = cna_input_ch.map { meta, tumor_bam, tumor_bai, norm_bam, norm_bai, snp_vcf, _snp_tbi, _sv_vcf ->
            tuple(meta, tumor_bam, tumor_bai, norm_bam, norm_bai, snp_vcf)
        }
        SAVANA_CNA(savana_input_ch, fasta, fai, contigs_list)
        ch_versions = ch_versions.mix(SAVANA_CNA.out.versions.first())
        ch_savana_cnv = SAVANA_CNA.out.cnv
    }

    emit:
    hp1_bed    = ch_hp1_bed // channel: [ val(meta), [ hp1_bed ] ]
    hp2_bed    = ch_hp2_bed // channel: [ val(meta), [ hp2_bed ] ]
    savana_cnv = ch_savana_cnv // channel: [ val(meta), [ cnv_tsv ] ]
    versions   = ch_versions // channel: [ versions.yml ]
}
