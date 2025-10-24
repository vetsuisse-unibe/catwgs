process indexgVCF {
    publishDir "${params.gVCF_folder}", mode: 'copy'
    tag "$params.tag"

    input:
              tuple val(bamID), path(vcf)
    output:
              tuple val(bamID), path("${bamID}.g.vcf.gz.tbi"), emit: indexedgVCFs
    script:
    """
    gatk --java-options "-Xmx5g" IndexFeatureFile -I ${vcf}

    """
}