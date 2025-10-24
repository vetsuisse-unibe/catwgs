#!/usr/bin/env nextflow 

nextflow.enable.dsl=2

process haplotypeCaller {

       tag "${bamID}-${interval.baseName}"

       input:
               tuple val(bamID), path(bam), path(bai), path(interval) 
               tuple path(ref), path(fai), path(dict)
       
       output:
               tuple  val(bamID), path("${interval.baseName}_${bamID}.g.vcf.gz"), emit: gvcfHaplotypeCaller
               tuple  val(bamID), path(interval), path("${interval.baseName}_${bamID}.g.vcf.gz"), emit: gvcfGenotypeGVCFs
       
       script:
       def tmp_dir = "\$TMPDIR/gatk_tmp_${bamID}_${interval.baseName}"
       """
       # Create and use a unique temporary directory
       mkdir -p ${tmp_dir}
       
       # Run HaplotypeCaller with more memory and proper tmp dir
       gatk --java-options "-Xmx${task.memory.toGiga()}g -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}" \\
            HaplotypeCaller \\
            --tmp-dir ${tmp_dir} \\
            -R ${ref} \\
            -I ${bam} \\
            -pairHMM AVX_LOGLESS_CACHING \\
            -L ${interval} \\
            --native-pair-hmm-threads ${task.cpus} \\
            -O ${interval.baseName}_${bamID}.g.vcf.gz \\
            -ERC GVCF \\
            -stand-call-conf 10

       # Clean up tmp directory
       rm -rf ${tmp_dir}
       """

}