// Define filterNonSNPs process
process filterNonSNPs {
    tag "Filtering non-SNPs from ${vcf}"
    scratch true
    publishDir params.filteredNonSnpVcfFolder ?: params.outdir, mode: 'copy'

    input:
    tuple path(vcf), path(vcf_index)
    path(ref)
    path(ref_fai)
    path(ref_dict)

    output:
    tuple path("${vcf.baseName.replaceAll(/\.rawNonSNPs/, '')}.rawNonSNPs.flt.vcf.gz"), path("${vcf.baseName.replaceAll(/\.rawNonSNPs/, '')}.rawNonSNPs.flt.vcf.gz.tbi"), emit: filtered_nonsnp_vcf

    script:
    """
    mkdir -p \$SCRATCH/tmp_filter
    
    gatk --java-options "-Xmx6G -Djava.io.tmpdir=\$SCRATCH/tmp_filter" VariantFiltration \\
        -R ${ref} \\
        -V ${vcf} \\
        --verbosity ERROR \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "FS > 200.0" --filter-name "FS200" \\
        -filter "ReadPosRankSum < -2.0" --filter-name "ReadPosRankSum-2" \\
        -filter "SOR > 10.0" --filter-name "SOR-10" \\
        -O ${vcf.baseName.replaceAll(/\.rawNonSNPs/, '')}.rawNonSNPs.flt.vcf.gz
    """
}
