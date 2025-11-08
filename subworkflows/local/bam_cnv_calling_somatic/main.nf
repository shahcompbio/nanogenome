// somatic copy number analysis
include { WAKHAN_REPHASE_CNA      } from '../../../modules/local/wakhan/rephase_cna/main'
include { ASCAT } from '../../../modules/local/ascat/main'

workflow BAM_CNV_CALLING_SOMATIC {

    take:
    tools           // val (cnv_tools)
    wakhan_input_ch // channel: [ val(meta), [ bam, bai, snp_vcf, snp_tbi, severus_vcf]]
    ascat_input_ch // channel: [ val(meta) [input_norm, idx_norm, input_tumor, idx_tumor] ]
    genomeVersion // val (genomeVersion)
    allele_files   // path (ascat_allele_files)
    loci_files   // path (ascat_loci_files)
    fasta           // path (fasta)
    fai             // path (fai)
    gc_file // path (ascat_gc_file)
    rt_file // path (ascat_rt_file)

    main:

    ch_versions = Channel.empty()


    if (tools.split(',').contains('wakhan')) {
        WAKHAN_REPHASE_CNA(wakhan_input_ch, fasta, fai)
        ch_versions = ch_versions.mix(WAKHAN_REPHASE_CNA.out.versions.first())
    }
    if (tools.split(',').contains('ascat')) {
        ASCAT(
        ascat_input_ch,
        genomeVersion,
        allele_files,
        loci_files,
        [],
        fasta,
        gc_file,
        rt_file
        )
        ch_versions = ch_versions.mix(WAKHAN_REPHASE_CNA.out.versions.first())
    }


    emit:

    hp1_bed  =   WAKHAN_REPHASE_CNA.out.HP1_bed // channel: [ val(meta), [ hp1_bed ] ]
    hp2_bed  =   WAKHAN_REPHASE_CNA.out.HP2_bed // channel: [ val(meta), [ hp2_bed ] ]
    allelefreqs = ASCAT.out.allelefreqs // channel: [ val(meta), [ allelefreq_file ] ]
    bafs = ASCAT.out.bafs // channel: [ val(meta), [ baf_file ] ]
    logRs = ASCAT.out.logrs // channel: [ val(meta), [ logR_file ] ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
