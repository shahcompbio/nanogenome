// convert sv tsv to bedpe
process TSV2BEDPE {
    tag "${meta.id}"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "preskaa/annotate_genes:v240817"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.bedpe"), emit: bedpe
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = meta.min_callers == 1 ? 'union' : 'consensus'
    """
    #!/usr/bin/env python3
    import pandas as pd
    strat_df = pd.read_csv("${tsv}", sep="\t")
    bedpe_df = strat_df[["minda_ID", "chrom1", "base1", "chrom2", "base2", "orientation", "SV_callers", "SV_Type"]].copy()
    bedpe_df["strand1"] = bedpe_df["orientation"].str[0]
    bedpe_df["strand2"] = bedpe_df["orientation"].str[1]
    bedpe_df["end1"] = bedpe_df["base1"] + 1
    bedpe_df["end2"] = bedpe_df["base2"] + 1
    bedpe_df = bedpe_df.drop(columns=["orientation"])
    # Add score column (number of supporting callers)
    bedpe_df["score"] = bedpe_df["SV_callers"].str.count(",") + 1
    # Combine minda_ID and SV_callers into OTHER column for VCF INFO field
    bedpe_df["OTHER"] = "minda_ID=" + bedpe_df["minda_ID"].astype(str) + ";SV_callers=" + bedpe_df["SV_callers"].astype(str)
    # Reorder columns to match SURVIVOR BEDPE format (TYPE must be column 7)
    bedpe_df = bedpe_df[["chrom1", "base1", "end1", "chrom2", "base2", "end2", "SV_Type", "score", "strand1", "strand2", "OTHER"]]
    bedpe_df.to_csv("${prefix}.${suffix}.bedpe", sep="\t", index=None, header=None)

    # Write versions
    import platform
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {platform.python_version()}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = meta.min_callers == 1 ? 'union' : 'consensus'
    """
    touch ${prefix}.${suffix}.bedpe

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
