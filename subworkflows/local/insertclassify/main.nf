// Classify somatic insertions from nanomonsv and severus
include { PREPINSERTTABLE          } from '../../../modules/local/prepinserttable/main'
include { NANOMONSV_INSERTCLASSIFY } from '../../../modules/local/nanomonsv/insertclassify/main'
include { VNTRCLASSIFY             } from '../../../modules/local/vntrclassify/main'

workflow INSERTCLASSIFY {
    take:
    annotated_sv_ch  // channel: [ val(meta), path(annotated_sv_tsv) ]
    nanomonsv_result // channel: [ val(meta), path(nanomonsv_result_txt) ]
    severus_vcf      // channel: [ val(meta), path(severus_vcf) ]
    ref_fasta        // path: reference genome FASTA
    bwa_fasta_index  // path: BWA index files (.amb, .ann, .bwt, .pac, .sa) co-located with ref_fasta
    ref_gtf          // path: gene annotation GTF
    line1_db         // path: LINE1 database BED
    vntr_bed         // path: VNTR BED file

    main:

    ch_versions = Channel.empty()

    // Combine inputs by sample ID for PREPINSERTTABLE
    prep_input_ch = annotated_sv_ch
        .map { meta, tsv -> [meta.id, meta, tsv] }
        .join(
            nanomonsv_result.map { meta, result -> [meta.id, result] },
            by: 0
        )
        .join(
            severus_vcf.map { meta, vcf -> [meta.id, vcf] },
            by: 0
        )
        .map { _id, meta, tsv, result, vcf ->
            [meta, tsv, result, vcf]
        }

    // Step 1: Prepare insertion table (resolve <INS> sequences)
    PREPINSERTTABLE(prep_input_ch)
    ch_versions = ch_versions.mix(PREPINSERTTABLE.out.versions.first())

    // Step 2: Run nanomonsv insert_classify
    NANOMONSV_INSERTCLASSIFY(
        PREPINSERTTABLE.out.inserts_tsv,
        ref_fasta,
        bwa_fasta_index,
        ref_gtf,
        line1_db,
    )
    ch_versions = ch_versions.mix(NANOMONSV_INSERTCLASSIFY.out.versions.first())

    // Step 3: VNTR classification and final categorization
    VNTRCLASSIFY(
        NANOMONSV_INSERTCLASSIFY.out.classified_tsv,
        vntr_bed,
    )
    ch_versions = ch_versions.mix(VNTRCLASSIFY.out.versions.first())

    emit:
    classified_inserts = VNTRCLASSIFY.out.final_classified_tsv // channel: [ val(meta), path(tsv) ]
    versions           = ch_versions // channel: [ path(versions.yml) ]
}
