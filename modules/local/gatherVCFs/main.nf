process gatherVCFs {

       tag "${bamID}"
       publishDir "${params.gVCF_folder}", mode: 'copy'

       input:
               tuple val(bamID), path(vcf)
       output:
               tuple val(bamID), path("${bamID}.g.vcf.gz"), emit: gatheredgVCFs
       script:
       //input = bam.collect{"-I ${it}"}.join(' ')
       """
       echo "${vcf.join('\n')}" | sort -V > ${bamID}.vcf.list
       gatk GatherVcfs -I ${bamID}.vcf.list -O ${bamID}.g.vcf.gz
       """
}