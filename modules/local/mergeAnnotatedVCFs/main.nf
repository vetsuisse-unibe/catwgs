// Define mergeAnnotatedVCFs process
process mergeAnnotatedVCFs {
    tag "Merging annotated SNPs and Non-SNPs"
    scratch true
    publishDir params.finalVcfFolder, mode: 'copy'

    input:
    path(snp_ann_vcf)
    path(nonsnp_ann_vcf)
    path(ref)
    path(fai)
    path(dict)

    output:
    tuple path("*.var.flt.ann.vcf.gz"), path("*.var.flt.ann.vcf.gz.tbi"), emit: final_annotated_vcf

    script:
    def prefix = snp_ann_vcf.baseName.replaceAll(/\.rawSNPs\.flt\.ann\.vcf/, '')
    """
    # Sort SNP VCF
    gatk --java-options "-Xmx18G" SortVcf \\
        -I ${snp_ann_vcf} \\
        -O ${prefix}.rawSNPs.flt.ann.sorted.vcf \\
        --TMP_DIR \$SCRATCH

    # Sort NonSNP VCF
    gatk --java-options "-Xmx18G" SortVcf \\
        -I ${nonsnp_ann_vcf} \\
        -O ${prefix}.rawNonSNPs.flt.ann.sorted.vcf \\
        --TMP_DIR \$SCRATCH

    # Merge sorted files
    gatk --java-options "-Xmx18G" MergeVcfs \\
        -I ${prefix}.rawSNPs.flt.ann.sorted.vcf \\
        -I ${prefix}.rawNonSNPs.flt.ann.sorted.vcf \\
        -O ${prefix}.var.flt.ann.vcf \\
        --TMP_DIR \$SCRATCH

    # Compress and index
    bgzip ${prefix}.var.flt.ann.vcf
    tabix -p vcf ${prefix}.var.flt.ann.vcf.gz
    """
}
