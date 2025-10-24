// Define annotateSNPs process
process annotateSNPs {
    tag "Annotating SNPs with snpEff"
    scratch true
    publishDir params.finalVcfFolder, mode: 'copy', pattern: "*.ann.vcf"
    publishDir params.finalVcfFolder, mode: 'copy', pattern: "*.stats.csv"

    input:
    tuple path(snp_vcf), path(snp_vcf_idx)
    val snpeff_path
    val snpeff_config
    val genome_version

    output:
    path("*.rawSNPs.flt.ann.vcf"), emit: annotated_snp_vcf
    path("*.SNPs.stats.csv"), emit: snp_stats

    script:
    def prefix = snp_vcf.baseName.replaceAll(/\.rawSNPs\.flt\.vcf/, '')
    """
    ${snpeff_path} -Xmx18g \\
        -csvStats ${prefix}.SNPs.stats.csv \\
        -c ${snpeff_config} \\
        -v ${genome_version} \\
        ${snp_vcf} \\
        > ${prefix}.rawSNPs.flt.ann.vcf
    """
}
