process VCF2MAF {
    tag "${meta.id}"
    label 'process_low'

    container "docker://quay.io/mondrianscwgs/variant_calling:v0.1.4"

    input:
    tuple val(meta), path(vcf), path(tumor_bam), path(norm_bam)
    path   ref_fasta
    path   vep_cache
    val    vep_fasta_suffix
    val    ncbi_build
    val    cache_version
    val    species

    output:
    tuple val(meta), path("${meta.id}.maf"), emit: maf
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    zcat ${vcf} > uncompressed.vcf

    rm -f uncompressed.vep.vcf

    # Create a vep wrapper that strips --af_gnomad (not available in this cache)
    mkdir -p vep_shim
    cat > vep_shim/vep << 'VEPEOF'
#!/bin/bash
args=()
for arg in "\$@"; do
    [ "\$arg" = "--af_gnomad" ] && continue
    args+=("\$arg")
done
exec /opt/conda/bin/vep "\${args[@]}"
VEPEOF
    chmod +x vep_shim/vep
    export PATH="\$PWD/vep_shim:\$PATH"

    vcf2maf uncompressed.vcf temp.maf \\
        ${ref_fasta} \\
        ${vep_cache} \\
        ${ncbi_build} \\
        ${cache_version} \\
        ${species} \\
        ${args}

    # Update MAF sample barcodes (BAMs lack RG headers so variant_utils cannot extract IDs)
    cat > fix_maf.py << 'PYEOF'
import sys
tumor_id  = sys.argv[1]
normal_id = sys.argv[2]
outfile   = sys.argv[3]
t_idx = n_idx = None
with open('temp.maf') as fin, open(outfile, 'w') as fout:
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue
        fields = line.rstrip('\\n').split('\\t')
        if fields[0] == 'Hugo_Symbol':
            try:
                t_idx = fields.index('Tumor_Sample_Barcode')
                n_idx = fields.index('Matched_Norm_Sample_Barcode')
            except ValueError:
                pass
            fout.write(line)
        else:
            if t_idx is not None:
                fields[t_idx] = tumor_id
            if n_idx is not None:
                fields[n_idx] = normal_id
            fout.write('\\t'.join(fields) + '\\n')
PYEOF
    python3 fix_maf.py "${meta.id}" "${meta.id}_NORMAL" "${meta.id}.maf"

    printf '"${task.process}":\\n    vcf2maf: unknown\\n' > versions.yml
    """
}
