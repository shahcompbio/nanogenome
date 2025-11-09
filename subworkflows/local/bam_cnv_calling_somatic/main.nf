// somatic copy number analysis
include { WAKHAN_REPHASE_CNA      } from '../../../modules/local/wakhan/rephase_cna/main'
include { UNZIP as UNZIP_ALLELES                    } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_LOCI                       } from '../../../modules/nf-core/unzip'
include { ASCAT } from '../../../modules/nf-core/ascat/main'

workflow BAM_CNV_CALLING_SOMATIC {

    take:
    tools           // val (cnv_tools)
    wakhan_input_ch // channel: [ val(meta), [ bam, bai, snp_vcf, snp_tbi, severus_vcf]]
    ascat_input_ch // channel: [ val(meta) [input_norm, idx_norm, input_tumor, idx_tumor] ]
    ascat_alleles   // path (ascat_allele_files)
    ascat_loci   // path (ascat_loci_files)
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
        // prepare ascat reference files
        if (!ascat_alleles) {
            allele_files = Channel.empty()
        } else if (ascat_alleles.endsWith(".zip")) {
            UNZIP_ALLELES(Channel.fromPath(file(ascat_alleles)).collect().map { it -> [[id: it[0].baseName], it] })
            allele_files = UNZIP_ALLELES.out.unzipped_archive.map { it[1] }
            ch_versions = ch_versions.mix(UNZIP_ALLELES.out.versions)
        } else {
            allele_files = Channel.fromPath(ascat_alleles).collect()
        }

        if (!ascat_loci) {
            loci_files = Channel.empty()
        } else if (ascat_loci.endsWith(".zip")) {
            UNZIP_LOCI(Channel.fromPath(file(ascat_loci)).collect().map { it -> [[id: it[0].baseName], it] })
            loci_files = UNZIP_LOCI.out.unzipped_archive.map { it[1] }
            ch_versions = ch_versions.mix(UNZIP_LOCI.out.versions)
        } else {
            loci_files = Channel.fromPath(ascat_loci).collect()
        }
        ASCAT(
        ascat_input_ch,
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
