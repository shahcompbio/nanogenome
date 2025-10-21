// merge and annotate structural variants from different callers
include { MINDA          } from '../../../modules/local/minda/main'
include { VCF2TSV        } from '../../../modules/local/vcf2tsv/main'
include { ANNOTATEGENES  } from '../../../modules/local/annotategenes/main'
include { CSVTK_CONCAT   } from '../../../modules/nf-core/csvtk/concat/main'
include { WGET           } from '../../../modules/nf-core/wget/main'
include { BIOMART        } from '../../../modules/local/biomart/main'
include { BCFTOOLS_QUERY } from '../../../modules/nf-core/bcftools/query/main'
include { SVTYPES        } from '../../../modules/local/svtypes/main'

workflow ANNOTATE_SV {
    take:
    sv_ch            // channel: [ val(meta), path(vcf1), path(vcf2), path(vcf3) ]; meta [id: sample id, min_callers: min number of callers]
    tolerance        // between breakpoints to consider two SVs the same
    min_size         // minimum size of SV to consider
    gene_annotations // gene annotation file
    oncokb           // oncokb annotation file
    oncokb_url       // oncokb url to download if oncokb not provided

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
    CSVTK_CONCAT(ANNOTATEGENES.out.annotated_sv.groupTuple(), "tsv", "tsv")
    // CSVTK_CONCAT.out.csv.view()
    ch_versions = ch_versions.mix(CSVTK_CONCAT.out.versions.first())
    // annotate sv types
    // collect raw calls for strand information
    sv_ch
        .flatMap { meta, vcfs ->
            vcfs.collect { vcf ->
                tuple(meta, vcf, [])
            }
        }
        .set { caller_ch }
    // take raw calls and make into tsv files
    BCFTOOLS_QUERY(caller_ch, [], [], [])
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())
    // group by number of callers used
    BCFTOOLS_QUERY.out.output
        .groupTuple()
        .join(CSVTK_CONCAT.out.csv)
        .set { calls_ch }
    // calls_ch.view()
    SVTYPES(calls_ch)
    ch_versions = ch_versions.mix(SVTYPES.out.versions.first())

    emit:
    minda_vcf    = MINDA.out.ensemble_vcf // channel: [ val(meta), [ vcf ] ]
    sv_table     = VCF2TSV.out.sv_table // channel: [ val(meta), path(sv_table) ]
    annotated_sv = SVTYPES.out.annotated_sv // channel: [ val(meta), path(annotated_sv) ]
    versions     = ch_versions // channel: [ versions.yml ]
}
