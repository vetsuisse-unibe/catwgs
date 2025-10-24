// Define filterSNPs process
process filterSNPs {
    tag "Filtering SNPs from ${vcf}"
    scratch true
    publishDir params.filteredSnpVcfFolder ?: params.outdir, mode: 'copy'

    input:
    tuple path(vcf), path(vcf_index)
    path(ref)
    path(ref_fai)
    path(ref_dict)

    output:
    tuple path("${vcf.baseName.replaceAll(/\.rawSNPs/, '')}.rawSNPs.flt.vcf.gz"), path("${vcf.baseName.replaceAll(/\.rawSNPs/, '')}.rawSNPs.flt.vcf.gz.tbi"), emit: filtered_snp_vcf

    script:
    """
    mkdir -p \$SCRATCH/tmp_filter
    
    gatk --java-options "-Xmx5g -Xms5g -Djava.io.tmpdir=\$SCRATCH/tmp_filter" VariantFiltration \\
        -R ${ref} \\
        -V ${vcf} \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "QUAL < 30.0" --filter-name "QUAL30" \\
        -filter "FS > 200.0" --filter-name "FS200" \\
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \\
        -O ${vcf.baseName.replaceAll(/\.rawSNPs/, '')}.rawSNPs.flt.vcf.gz
    """
}
