// Define selectSNP process
process selectSNP {
    tag "Selecting SNPs from ${vcf}"
    scratch true
    publishDir params.snpVcfFolder ?: params.outdir, mode: 'copy'

    input:
    tuple path(vcf), path(vcf_index)
    path(ref)
    path(ref_fai)
    path(ref_dict)

    output:
    tuple path("${vcf.baseName.replaceAll(/\.vcf/, '')}.rawSNPs.vcf.gz"), path("${vcf.baseName.replaceAll(/\.vcf/, '')}.rawSNPs.vcf.gz.tbi"), emit: snp_vcf

    script:
    """
    mkdir -p \$SCRATCH/tmp_select
    
    gatk --java-options "-Xmx6G -Djava.io.tmpdir=\$SCRATCH/tmp_select" SelectVariants \\
        -R ${ref} \\
        --select-type-to-include SNP \\
        -V ${vcf} \\
        -O ${vcf.baseName.replaceAll(/\.vcf/, '')}.rawSNPs.vcf.gz
    """
}
