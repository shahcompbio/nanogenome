// merge and annotate structural variants from different callers
include { MINDA            } from '../../../modules/local/minda/main'
include { TABIX_BGZIPTABIX } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { VCF2TSV          } from '../../../modules/local/vcf2tsv/main'
include { CSVTK_SPLIT      } from '../../../modules/nf-core/csvtk/split/main'
include { ANNOTATEGENES    } from '../../../modules/local/annotategenes/main'
workflow ANNOTATE_SV {
    take:
    sv_ch            // channel: [ val(meta), path(vcf1), path(vcf2), path(vcf3) ]; meta [id: sample id, min_callers: min number of callers]
    tolerance        // between breakpoints to consider two SVs the same
    min_size         // minimum size of SV to consider
    gene_annotations // gene annotation file
    oncokb           // oncokb annotation file

    main:

    ch_versions = Channel.empty()
    // run minda to combine SV
    MINDA(sv_ch, params.tolerance, params.min_size)
    ch_versions = ch_versions.mix(MINDA.out.versions.first())
    // index vcfs
    TABIX_BGZIPTABIX(MINDA.out.ensemble_vcf)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())
    // convert to tsv and chunk for gene annotation
    VCF2TSV(MINDA.out.ensemble_vcf)
    ch_versions = ch_versions.mix(VCF2TSV.out.versions.first())
    // split chunks into
    // annotate chunked tsv files
    // Spread chunks into individual channel items
    VCF2TSV.out.chunks
        .flatMap { meta, chunks -> chunks.collect { chunk -> [meta, chunk] } }
        .set { tsv_chunks }
    tsv_chunks.view()
    ANNOTATEGENES(tsv_chunks, gene_annotations, oncokb)
    print("annotating svs")

    emit:
    minda_vcf = MINDA.out.ensemble_vcf // channel: [ val(meta), [ vcf ] ]
    sv_table  = VCF2TSV.out.sv_table // channel: [ val(meta), path(sv_table) ]
    versions  = ch_versions // channel: [ versions.yml ]
}
