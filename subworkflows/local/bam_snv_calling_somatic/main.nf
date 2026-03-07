// run somatic snv & indel calling
include { CLAIRS              } from '../../../modules/local/clairs/main'
include { GUNZIP              } from '../../../modules/nf-core/gunzip/main'
include { DEEPSOMATIC         } from '../../../modules/nf-core/deepsomatic/main'
include { VCF2MAF             } from '../../../modules/nf-core/vcf2maf/main'
include { FIXMAFBARCODES      } from '../../../modules/local/fixmafbarcodes/main'
include { ENSEMBLVEP_DOWNLOAD } from '../../../modules/nf-core/ensemblvep/download/main'

workflow BAM_SNV_CALLING_SOMATIC {
    take:
    snv_input_ch // channel: [ val(meta), tumor_bam, tumor_bai, norm_bam, norm_bai ]
    fasta // /path/to/ref_fasta
    fai // /path/to ref_fai
    vep_cache // /path/to/vep_cache
    vep_assembly // parameter if downloading vep
    vep_species // specify species for vep download
    vep_cache_version // specify vep cache version
    somatic_snv_caller // snv/indel callers to use
    inhibit_vep // skip VEP annotation

    main:
    ch_versions = channel.empty()
    vcf_ch = channel.empty()
    // run somatic snv/indel calling
    if (somatic_snv_caller == "clairs") {
        CLAIRS(snv_input_ch, fasta, fai)
        // tag SNV and indel VCFs with meta.variant and combine into single channel
        snv_ch = CLAIRS.out.snv_vcf.map { meta, vcf ->
            tuple(meta + [variant: 'snv'], vcf)
        }
        indel_ch = CLAIRS.out.indel_vcf.map { meta, vcf ->
            tuple(meta + [variant: 'indel'], vcf)
        }
        vcf_ch = snv_ch.mix(indel_ch)
    }
    else if (somatic_snv_caller == "deepsomatic") {
        DEEPSOMATIC(
            snv_input_ch,
            [[id: "ref"], []],
            [[id: "ref"], fasta],
            [[id: "ref"], fai],
            [[id: "ref"], []],
        )
        vcf_ch = DEEPSOMATIC.out.vcf
    }
    else {
        error("Unsupported somatic SNV caller: ${somatic_snv_caller}. Options are: clairs, deepsomatic")
    }
    // gunzip for vcf2maf
    GUNZIP(vcf_ch)
    // download VEP cache if not provided and VEP is not inhibited
    if (!inhibit_vep && !vep_cache) {
        ENSEMBLVEP_DOWNLOAD(
            [
                [id: "vep"],
                vep_assembly,
                vep_species,
                vep_cache_version,
            ]
        )
        vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache
    }
    // vcf to maf process
    VCF2MAF(
        GUNZIP.out.gunzip,
        fasta,
        inhibit_vep ? [] : vep_cache,
    )
    ch_versions = ch_versions.mix(VCF2MAF.out.versions.first())
    FIXMAFBARCODES(VCF2MAF.out.maf)

    emit:
    maf      = FIXMAFBARCODES.out.maf // channel [ val(meta), [maf]]
    versions = ch_versions
}
