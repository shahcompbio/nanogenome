/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PHASING                 } from '../subworkflows/local/phasing/main'
include { SV_CALLING_SOMATIC      } from '../subworkflows/local/sv_calling_somatic/main'
include { ANNOTATE_SV             } from '../subworkflows/local/annotate_sv/main'
include { SV_CALLING_GERMLINE     } from '../subworkflows/local/sv_calling_germline/main'
include { BAM_CNV_CALLING_SOMATIC } from '../subworkflows/local/bam_cnv_calling_somatic/main'
include { PLOTCIRCOS              } from '../modules/local/plotcircos/main'
include { SVKARYOPLOT             } from '../modules/local/svkaryoplot/main'
include { MULTIQC                 } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nfcore_nanogenome_pipeline'

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
    // define contigs list when null
    contigs_list = params.contigs_list ?: []
    /*
    * PHASING WORKFLOW
    */
    // run phasing subworkflow to phase variants and haplotag bams
    if (params.skip_somatic && !params.germline && params.skip_cna) {
        println("running phasing workflow only")
    }
    if (!params.skip_phasing) {
        PHASING(
            ch_samplesheet,
            params.clair3_model,
            params.clair3_platform,
            params.fasta,
            params.fai,
        )
        ch_versions = ch_versions.mix(PHASING.out.versions)
        // add phasing stats to multiqc
        ch_multiqc_files = ch_multiqc_files.mix(PHASING.out.whatshap_stats.collect { it[1] })
        // branch bam channel for downstream processing
        bam_ch = PHASING.out.bam
            .join(PHASING.out.bai, by: 0)
            .branch { meta, bam, bai ->
                tumor: meta.condition == 'tumor'
                norm: meta.condition == 'normal'
            }
        // construct snps channel for downstream processing
        snps_ch = PHASING.out.phased_vcf.map { meta, vcf -> tuple(meta, vcf, []) }
        // branch phasing results for cna input channel
        bam_snps_ch = PHASING.out.bam_snps.branch { meta, bam, bai, vcf, tbi ->
            tumor: meta.condition == 'tumor'
            norm: meta.condition == 'normal'
        }
    }
    else if (!params.skip_somatic || params.germline) {
        println("running variant calling")
        // split samplesheet into tumor/normal
        // ch_samplesheet.view()
        bam_ch = ch_samplesheet
            .map { meta, bam, bai, _snp_vcf, _snp_tbi, _sv_vcf ->
                tuple(meta, bam, bai)
            }
            .branch { meta, _bam, _bai ->
                tumor: meta.condition == 'tumor'
                norm: meta.condition == 'normal'
            }
        // construct snps channel
        snps_ch = ch_samplesheet
            .map { meta, _bam, _bai, snp_vcf, snp_tbi, _sv_vcf ->
                tuple(meta, snp_vcf, snp_tbi)
            }
            .filter { _meta, snp_vcf, _snp_tbi ->
                snp_vcf != null && snp_vcf != []
            }
        // snps_ch.view()
        // also construct bam/snps channel for wakhan
        // ch_samplesheet.view()
        bam_snps_ch = ch_samplesheet
            .map { meta, bam, bai, _snp_vcf, _snp_tbi, _sv_vcf -> tuple(meta.id, meta, bam, bai) }
            .combine(
                snps_ch.map { meta, snp_vcf, snp_tbi -> tuple(meta.id, snp_vcf, snp_tbi) },
                by: 0
            )
            .map { _id, meta, bam, bai, snp_vcf, snp_tbi ->
                tuple(meta, bam, bai, snp_vcf, snp_tbi)
            }
            .branch { meta, _bam, _bai, _snp_vcf, _snp_tbi ->
                tumor: meta.condition == 'tumor'
                norm: meta.condition == 'normal'
            }
    }
    else {
        println("run cna only")
    }

    // make default channels
    sv_ch = channel.empty()
    cna_input_ch = channel.empty()
    support_ch = channel.from(1, params.min_callers)
    /*
    * SOMATIC STRUCTURAL VARIANT CALLING WORKFLOW
    */
    // call somatic structural variants
    if (!params.skip_somatic) {
        // validate that vntr_bed is provided if severus is selected
        if (params.somatic_callers.split(',').contains('severus') && !params.vntr_bed) {
            error("ERROR: --vntr_bed is required when using SEVERUS for somatic SV calling. Please provide a VNTR BED file.")
        }
        // construct somatic sv input channel
        input_somatic_ch = bam_ch.tumor
            .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
            .join(bam_ch.norm.map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }, by: 0)
            .join(snps_ch.map { meta, vcf, _tbi -> tuple(meta.id, meta, vcf) }, by: 0)
            .map { id, tumor_meta, tumor_bam, tumor_bai, _norm_meta, norm_bam, norm_bai, _meta3, vcf ->
                tuple([id: id], tumor_bam, tumor_bai, norm_bam, norm_bai, vcf)
            }
        SV_CALLING_SOMATIC(
            params.somatic_callers,
            input_somatic_ch,
            bam_ch.tumor.mix(bam_ch.norm),
            params.vntr_bed,
            params.fasta,
            params.fai,
            contigs_list,
        )
        ch_versions = ch_versions.mix(SV_CALLING_SOMATIC.out.versions)
        // construct sv channel for annotation subworkflow
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
    }
    /*
    * GERMLINE SV CALLING WORKFLOW
    */
    // call germline structural variants
    if (params.germline) {
        // validate that vntr_bed is provided if severus is selected
        if (params.germline_callers.split(',').contains('severus') && !params.vntr_bed) {
            error("ERROR: --vntr_bed is required when using SEVERUS for germline SV calling. Please provide a VNTR BED file.")
        }
        // sv caller input channel
        // ch_samplesheet.view { v -> "samplesheet ${v}" }
        // bam_snps_ch.norm.view { v -> "normal bam_snps ${v}" }
        // bam_snps_ch.tumor.view { v -> "tumor bam_snps ${v}" }
        input_germline_ch = bam_snps_ch.norm.map { meta, bam, bai, snp_vcf, _snp_tbi ->
            tuple(meta, bam, bai, snp_vcf)
        }
        // input_germline_ch.view { v -> "input germline ${v}" }
        SV_CALLING_GERMLINE(
            params.germline_callers,
            input_germline_ch,
            params.vntr_bed,
            params.fasta,
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
    /*
    * SOMATIC CNV CALLING WORKFLOW
    */
    // run somatic cnv analysis
    hp1_bed_ch = Channel.empty()
    hp2_bed_ch = Channel.empty()
    if (!params.skip_cna) {
        // construct cna input channel
        if (!params.skip_somatic) {
            // construct somatic sv input channel
            cna_input_ch = bam_ch.tumor
                .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
                .join(bam_ch.norm.map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }, by: 0)
                .join(snps_ch.map { meta, snp_vcf, snp_tbi -> tuple(meta.id, meta, snp_vcf, snp_tbi) }, by: 0)
                .join(
                    SV_CALLING_SOMATIC.out.severus_vcf.map { meta, sv_vcf ->
                        tuple(meta.id, meta, sv_vcf)
                    },
                    by: 0
                )
                .map { id, _tumor_meta, tumor_bam, tumor_bai, _norm_meta, norm_bam, norm_bai, _meta3, snp_vcf, snp_tbi, _meta4, sv_vcf ->
                    tuple([id: id, condition: "somatic"], tumor_bam, tumor_bai, norm_bam, norm_bai, snp_vcf, snp_tbi, sv_vcf)
                }
        }
        else {
            // branch the samplesheet
            cna_ch = ch_samplesheet.branch { meta, _bam, _bai, _snp_vcf, _snp_tbi, _severus_vcf ->
                tumor: meta.condition == 'tumor'
                norm: meta.condition == 'normal'
            }
            // contruct cna input ch
            cna_input_ch = cna_ch.tumor
                .map { meta, bam, bai, _snp_vcf, _snp_tbi, sv_vcf ->
                    tuple(meta.id, bam, bai, sv_vcf)
                }
                .join(
                    cna_ch.norm.map { meta, bam, bai, snp_vcf, snp_tbi, _sv_vcf ->
                        tuple(meta.id, bam, bai, snp_vcf, snp_tbi)
                    },
                    by: 0
                )
                .map { id, tumor_bam, tumor_bai, sv_vcf, norm_bam, norm_bai, snp_vcf, snp_tbi ->
                    tuple([id: id, condition: "somatic"], tumor_bam, tumor_bai, norm_bam, norm_bai, snp_vcf, snp_tbi, sv_vcf)
                }
        }
        // run somatic cnv calling subworkflow
        BAM_CNV_CALLING_SOMATIC(
            params.cna_tools,
            cna_input_ch,
            params.fasta,
            params.fai,
            contigs_list,
        )
        ch_versions = ch_versions.mix(BAM_CNV_CALLING_SOMATIC.out.versions)
        hp1_bed_ch = BAM_CNV_CALLING_SOMATIC.out.hp1_bed
        hp2_bed_ch = BAM_CNV_CALLING_SOMATIC.out.hp2_bed
    }
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
        // circos plot for somatic SVs
        if (!params.skip_cna) {
            circos_ch = annot_sv_ch.somatic
                .combine(
                    hp1_bed_ch.map { meta, hp_bed ->
                        tuple("${meta.id}-${meta.condition}", meta, hp_bed)
                    },
                    by: 0
                )
                .combine(
                    hp2_bed_ch.map { meta, hp_bed ->
                        tuple("${meta.id}-${meta.condition}", meta, hp_bed)
                    },
                    by: 0
                )
                .map { id, meta, sv, meta1, hp1, meta2, hp2 ->
                    tuple(meta, sv, hp1, hp2)
                }
        }
        else {
            circos_ch = annot_sv_ch.somatic.map { _id, meta, sv -> [meta, sv, [], []] }
        }
        PLOTCIRCOS(circos_ch)
        ch_versions = ch_versions.mix(PLOTCIRCOS.out.versions.first())
        // karyoplot for germline SVs
        SVKARYOPLOT(
            annot_sv_ch.germline.map { _meta, meta1, sv ->
                tuple(meta1, sv)
            },
            params.genome_build,
        )
        ch_versions = ch_versions.mix(SVKARYOPLOT.out.versions.first())
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
