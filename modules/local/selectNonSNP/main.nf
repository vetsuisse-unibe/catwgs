// Define selectNonSNP process
process selectNonSNP {
    tag "Selecting non-SNPs from ${vcf}"
    scratch true
    publishDir params.nonSnpVcfFolder ?: params.outdir, mode: 'copy'

    input:
    tuple path(vcf), path(vcf_index)
    path(ref)
    path(ref_fai)
    path(ref_dict)

    output:
    tuple path("${vcf.baseName.replaceAll(/\.vcf/, '')}.rawNonSNPs.vcf.gz"), path("${vcf.baseName.replaceAll(/\.vcf/, '')}.rawNonSNPs.vcf.gz.tbi"), emit: nonsnp_vcf

    script:
    """
    mkdir -p \$SCRATCH/tmp_select
    
    gatk --java-options "-Xmx6G -Djava.io.tmpdir=\$SCRATCH/tmp_select" SelectVariants \\
        -R ${ref} \\
        --select-type-to-exclude SNP \\
        -V ${vcf} \\
        -O ${vcf.baseName.replaceAll(/\.vcf/, '')}.rawNonSNPs.vcf.gz
    """
}
