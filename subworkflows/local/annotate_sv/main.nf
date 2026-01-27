// merge and annotate structural variants from different callers
include { MINDA                      } from '../../../modules/local/minda/main'
include { VCF2TSV                    } from '../../../modules/local/vcf2tsv/main'
include { ANNOTATEGENES              } from '../../../modules/local/annotategenes/main'
include { CSVTK_CONCAT               } from '../../../modules/nf-core/csvtk/concat/main'
include { WGET                       } from '../../../modules/nf-core/wget/main'
include { BIOMART                    } from '../../../modules/local/biomart/main'
include { ANNOTSV_INSTALLANNOTATIONS } from '../../../modules/nf-core/annotsv/installannotations/main'
include { ANNOTSV_ANNOTSV            } from '../../../modules/nf-core/annotsv/annotsv/main'
include { TSV2BEDPE                  } from '../../../modules/local/tsv2bedpe/main'
include { SURVIVOR_BEDPETOVCF        } from '../../../modules/nf-core/survivor/bedpetovcf/main'

workflow ANNOTATE_SV {
    take:
    sv_ch // channel: [ val(meta), path(vcf1), path(vcf2), path(vcf3) ]; meta [id: sample id, min_callers: min number of callers]
    tolerance // between breakpoints to consider two SVs the same
    min_size // minimum size of SV to consider
    gene_annotations // gene annotation file
    oncokb // oncokb annotation file
    oncokb_url // oncokb url to download if oncokb not provided
    skip_annotsv // boolean to run annotsv
    annotsv_dir // annotsv annotation directory

    main:

    ch_versions = Channel.empty()
    // run minda to combine SV
    MINDA(sv_ch, params.tolerance, params.min_size)
    ch_versions = ch_versions.mix(MINDA.out.versions.first())
    // convert to tsv and chunk for gene annotation
    VCF2TSV(MINDA.out.ensemble_vcf)
    ch_versions = ch_versions.mix(VCF2TSV.out.versions.first())
    // split chunks into
    // annotate chunked tsv files
    // Spread chunks into individual channel items
    VCF2TSV.out.chunks
        .flatMap { meta, chunks -> chunks.collect { chunk -> [meta, chunk] } }
        .set { tsv_chunks }
    // download oncokb if not provided
    if (!oncokb) {
        WGET([[id: 'oncokb'], oncokb_url])
        oncokb = WGET.out.outfile.map { meta, path -> path }
        ch_versions = ch_versions.mix(WGET.out.versions)
    }
    // build gene annotation table if not provided
    if (!gene_annotations) {
        BIOMART()
        gene_annotations = BIOMART.out.gene_annotation
        ch_versions = ch_versions.mix(BIOMART.out.versions)
    }
    // annotate genes
    ANNOTATEGENES(tsv_chunks, gene_annotations, oncokb)
    ch_versions = ch_versions.mix(ANNOTATEGENES.out.versions.first())
    // combine the annotated chunks back into single file per sample
    // ANNOTATEGENES.out.annotated_sv.view()
    CSVTK_CONCAT(ANNOTATEGENES.out.annotated_sv.groupTuple(), "tsv", "tsv")
    // CSVTK_CONCAT.out.csv.view()
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions.first())
    // run annotsv to annotate sv with additional annotations
    if (!skip_annotsv) {
        if (!annotsv_dir) {
            ANNOTSV_INSTALLANNOTATIONS()
            ch_annotsv_dir = ANNOTSV_INSTALLANNOTATIONS.out.annotations.map { dir -> [[id: 'annotations'], dir] }
            ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions)
        }
        else {
            ch_annotsv_dir = channel.value([[id: 'annotations'], annotsv_dir])
        }
        // run annotsv
        TSV2BEDPE(CSVTK_CONCAT.out.csv)
        // TSV2BEDPE.out.bedpe.view()
        SURVIVOR_BEDPETOVCF(TSV2BEDPE.out.bedpe)
        ch_versions = ch_versions.mix(SURVIVOR_BEDPETOVCF.out.versions.first())
        // SURVIVOR_BEDPETOVCF.out.vcf.view()
        annotsv_in_ch = SURVIVOR_BEDPETOVCF.out.vcf.map { meta, vcf -> [meta, vcf, [], []] }
        ANNOTSV_ANNOTSV(
            annotsv_in_ch,
            ch_annotsv_dir,
            [[], []],
            [[], []],
            [[], []],
        )
        ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions.first())
    }

    emit:
    minda_vcf    = MINDA.out.ensemble_vcf // channel: [ val(meta), [ vcf ] ]
    sv_table     = VCF2TSV.out.sv_table // channel: [ val(meta), path(sv_table) ]
    annotated_sv = CSVTK_CONCAT.out.csv // channel: [ val(meta), path(annotated_sv) ]
    versions     = ch_versions // channel: [ versions.yml ]
}
