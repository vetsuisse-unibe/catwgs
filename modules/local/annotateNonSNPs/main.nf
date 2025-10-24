// Define annotateNonSNPs process
process annotateNonSNPs {
    tag "Annotating Non-SNPs with snpEff"
    scratch true
    publishDir params.finalVcfFolder, mode: 'copy', pattern: "*.ann.vcf"
    publishDir params.finalVcfFolder, mode: 'copy', pattern: "*.stats.csv"

    input:
    tuple path(nonsnp_vcf), path(nonsnp_vcf_idx)
    val snpeff_path
    val snpeff_config
    val genome_version

    output:
    path("*.rawNonSNPs.flt.ann.vcf"), emit: annotated_nonsnp_vcf
    path("*.nonSNPs.stats.csv"), emit: nonsnp_stats

    script:
    def prefix = nonsnp_vcf.baseName.replaceAll(/\.rawNonSNPs\.flt\.vcf/, '')
    """
    ${snpeff_path} -Xmx18g \\
        -csvStats ${prefix}.nonSNPs.stats.csv \\
        -c ${snpeff_config} \\
        -v ${genome_version} \\
        ${nonsnp_vcf} \\
        > ${prefix}.rawNonSNPs.flt.ann.vcf
    """
}
