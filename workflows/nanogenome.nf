/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { HAPLOTAG               } from '../subworkflows/local/haplotag/main'
include { SV_CALLING_SOMATIC     } from '../subworkflows/local/sv_calling_somatic/main'
include { ANNOTATE_SV            } from '../subworkflows/local/annotate_sv/main'
include { WAKHAN_REPHASE_CNA     } from '../modules/local/wakhan/rephase_cna/main'
include { SV_CALLING_GERMLINE    } from '../subworkflows/local/sv_calling_germline/main'
include { PLOTCIRCOS             } from '../modules/local/plotcircos/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nanogenome_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NANOGENOME {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // run haplotag subworkflow to haplotag bams
    if (params.skip_somatic && !params.germline) {
        println("running phasing workflow only")
    }
    HAPLOTAG(
        ch_samplesheet,
        params.clair3_model,
        params.clair3_platform,
        params.fasta,
        params.fai,
    )
    ch_versions = ch_versions.mix(HAPLOTAG.out.versions)
    // construct input for sv calling subworkflows
    haplotag_out_ch = HAPLOTAG.out.bam
        .join(HAPLOTAG.out.bai, by: 0)
        .join(HAPLOTAG.out.phased_vcf, by: 0)
    haplotag_out_ch.view()
    // add phasing stats to multiqc
    ch_multiqc_files = ch_multiqc_files.mix(HAPLOTAG.out.whatshap_stats.collect { it[1] })
    // make default channels
    sv_ch = Channel.empty()
    cna_input_ch = Channel.empty()
    // call somatic structural variants
    if (!params.skip_somatic) {
        SV_CALLING_SOMATIC(
            params.somatic_callers,
            HAPLOTAG.out.bam,
            HAPLOTAG.out.bai,
            HAPLOTAG.out.phased_vcf,
            params.vntr_bed,
            params.fasta,
            params.fai,
        )
        ch_versions = ch_versions.mix(SV_CALLING_SOMATIC.out.versions)
        // construct sv channel for annotation subworkflow
        support_ch = Channel.from(1, params.min_callers)
        sv_ch = SV_CALLING_SOMATIC.out.savana_vcf
            .join(SV_CALLING_SOMATIC.out.severus_vcf, by: 0)
            .join(SV_CALLING_SOMATIC.out.nanomonsv_vcf, by: 0)
            .combine(support_ch)
            .map { meta, vcf1, vcf2, vcf3, min_callers ->
                [
                    [
                        id: meta.id,
                        condition: "somatic",
                        min_callers: min_callers,
                    ],
                    [vcf1, vcf2, vcf3],
                ]
            }
        // branch haplotag results for cna input channel
        hap_bam_snps = HAPLOTAG.out.bam_snps.branch { meta, bam, bai, vcf, tbi ->
            tumor: meta.condition == 'tumor'
            norm: meta.condition == 'normal'
        }
        // construct cna input channel
        cna_input_ch = hap_bam_snps.tumor
            .join(
                SV_CALLING_SOMATIC.out.severus_vcf.map { meta, vcf ->
                    [meta + [condition: "tumor"], vcf]
                },
                by: 0
            )
            .map { meta, bam, bai, snp_vcf, snp_tbi, sv_vcf ->
                tuple([id: "${meta.id}", condition: "somatic"], bam, bai, snp_vcf, snp_tbi, sv_vcf)
            }
    }

    // run germline workflow
    // call germline structural variants
    if (params.germline) {
        SV_CALLING_GERMLINE(
            params.germline_callers,
            HAPLOTAG.out.bam,
            HAPLOTAG.out.bai,
            HAPLOTAG.out.phased_vcf,
            params.vntr_bed,
            params.fasta,
            ch_samplesheet,
        )
        ch_versions = ch_versions.mix(SV_CALLING_GERMLINE.out.versions)
        // construct channel of germline calls + mix with sv channel
        germline_ch = SV_CALLING_GERMLINE.out.severus_vcf
            .join(SV_CALLING_GERMLINE.out.cutesv_vcf, by: 0)
            .join(SV_CALLING_GERMLINE.out.sniffles_vcf, by: 0)
            .join(SV_CALLING_GERMLINE.out.longcalld_vcf, by: 0)
            .combine(support_ch)
            .map { meta, vcf1, vcf2, vcf3, vcf4, min_callers ->
                [
                    [
                        id: meta.id,
                        condition: "germline",
                        min_callers: min_callers,
                    ],
                    [vcf1, vcf2, vcf3, vcf4],
                ]
            }
        sv_ch = sv_ch.mix(germline_ch)
    }
    // sv_ch.view()
    // run annotation only if sv calling has been performed
    if (!params.skip_somatic || params.germline) {
        // run merge + annotate SV subworkflow
        ANNOTATE_SV(
            sv_ch,
            params.tolerance,
            params.min_size,
            params.gene_annotations,
            params.oncokb,
            params.oncokb_url,
        )
        ch_versions = ch_versions.mix(ANNOTATE_SV.out.versions)
        // run wakhan for somatic CNA
        // cna_input_ch.view()
        WAKHAN_REPHASE_CNA(cna_input_ch, params.fasta)
        ch_versions = ch_versions.mix(WAKHAN_REPHASE_CNA.out.versions.first())

        // plot results
        // ANNOTATE_SV.out.annotated_sv.view()
        ANNOTATE_SV.out.annotated_sv
            .map { meta, sv ->
                tuple("${meta.id}-${meta.condition}", meta, sv)
            }
            .branch { meta, meta1, sv ->
                somatic: meta1.condition == 'somatic'
                germline: meta1.condition == 'germline'
            }
            .set { annot_sv_ch }
        // WAKHAN_CNA.out.HP1_bed.view()
        circos_ch = annot_sv_ch.somatic
            .combine(
                WAKHAN_REPHASE_CNA.out.HP1_bed.map { meta, hp_bed ->
                    tuple("${meta.id}-${meta.condition}", meta, hp_bed)
                },
                by: 0
            )
            .combine(
                WAKHAN_REPHASE_CNA.out.HP2_bed.map { meta, hp_bed ->
                    tuple("${meta.id}-${meta.condition}", meta, hp_bed)
                },
                by: 0
            )
            .map { id, meta, sv, meta1, hp1, meta2, hp2 ->
                tuple(meta, sv, hp1, hp2)
            }
            .concat(
                annot_sv_ch.germline.map { meta, meta1, sv ->
                    [meta1, sv, [], []]
                }
            )
        // circos_ch.view()
        PLOTCIRCOS(circos_ch)
        ch_versions = ch_versions.mix(PLOTCIRCOS.out.versions.first())
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nanogenome_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
