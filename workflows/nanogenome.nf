/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { HAPLOTAG               } from '../subworkflows/local/haplotag/main'
include { SV_CALLING_SOMATIC     } from '../subworkflows/local/sv_calling_somatic/main'
include { ANNOTATE_SV            } from '../subworkflows/local/annotate_sv/main'
include { WAKHAN_CNA             } from '../modules/local/wakhan/cna/main'
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
    HAPLOTAG(ch_samplesheet, params.clair3_model, params.clair3_platform, params.fasta, params.fai)
    ch_versions = ch_versions.mix(HAPLOTAG.out.versions)
    // add phasing stats to multiqc
    ch_multiqc_files = ch_multiqc_files.mix(HAPLOTAG.out.whatshap_stats.collect { it[1] })
    // call somatic structural variants
    SV_CALLING_SOMATIC(
        params.somatic_callers,
        HAPLOTAG.out.bam,
        HAPLOTAG.out.bai,
        HAPLOTAG.out.rephased_vcf,
        params.vntr_bed,
        params.fasta,
        params.fai,
    )
    ch_versions = ch_versions.mix(SV_CALLING_SOMATIC.out.versions)
    // run germline workflow
    // call germline structural variants
    if (params.germline == true) {
        SV_CALLING_GERMLINE(
            params.germline_callers,
            HAPLOTAG.out.bam,
            HAPLOTAG.out.bai,
            HAPLOTAG.out.rephased_vcf,
            params.vntr_bed,
            params.fasta,
            ch_samplesheet,
        )
        ch_versions = ch_versions.mix(SV_CALLING_GERMLINE.out.versions)
    }
    // merge and annotate SVs in different callers and generate both union and consensus VCFs
    support_ch = Channel.from(1, params.min_callers)
    sv_ch = SV_CALLING_SOMATIC.out.savana_vcf
        .join(SV_CALLING_SOMATIC.out.severus_vcf, by: 0)
        .join(SV_CALLING_SOMATIC.out.nanomonsv_vcf, by: 0)
        .combine(support_ch)
        .map { meta, vcf1, vcf2, vcf3, min_callers ->
            [meta + [min_callers: min_callers], vcf1, vcf2, vcf3]
        }
    sv_ch.view()
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
    // run wakhan cna
    cna_input_ch = HAPLOTAG.out.bam_snps
        .branch { meta, bam, bai, vcf, tbi ->
            tumor: meta.condition == 'tumor'
        }
        .map { meta, bam, bai, vcf, tbi ->
            tuple([id: meta.id], bam, bai, vcf, tbi)
        }
        .join(SV_CALLING_SOMATIC.out.severus_vcf, by: 0)
        .join(
            HAPLOTAG.out.wakhanHPOutput.map { meta, path ->
                tuple([id: meta.id], path)
            },
            by: 0
        )
    // cna_input_ch.view()
    WAKHAN_CNA(cna_input_ch, params.fasta)
    ch_versions = ch_versions.mix(WAKHAN_CNA.out.versions.first())
    // plot results
    circos_ch = ANNOTATE_SV.out.annotated_sv
        .map { meta, sv -> tuple([id: meta.id], meta, sv) }
        .combine(WAKHAN_CNA.out.HP1_bed, by: 0)
        .combine(WAKHAN_CNA.out.HP2_bed, by: 0)
        .map { meta, meta1, sv, hp1, hp2 ->
            tuple(meta1, sv, hp1, hp2)
        }
    // circos_ch.view()
    PLOTCIRCOS(circos_ch)
    ch_versions = ch_versions.mix(PLOTCIRCOS.out.versions.first())
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
