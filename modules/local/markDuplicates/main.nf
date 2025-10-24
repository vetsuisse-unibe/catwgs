#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process markDuplicates {
    tag "${sampleID}"
    publishDir "${params.dedup}", mode: 'symlink', overwrite: true
    scratch true  // Use scratch space for the entire process

    input:
        tuple val(sampleID), path(bam), path(bai)
        val(dedup_dir)

    output:
        tuple val(sampleID), path(outfile_bam), path(outfile_bai), emit: duplicateMarkedBam
        path(outfile_metrics), emit: markDuplicatesQC

    script:
        outfile_bam = sampleID + ".dedup.bam"
        outfile_bai = sampleID + ".dedup.bai"
        outfile_metrics = sampleID + "_duplicate_metrics.txt"
        def avail_mem = task.memory ? task.memory.toGiga() - 2 : 8
        def java_mem = avail_mem > 2 ? "${avail_mem}G" : "2G"
        tmp_dir = "\$PWD/tmp_markdup_${sampleID}"
        
        """
        # Create and use a unique temporary directory
        mkdir -p ${tmp_dir}

        # Print current working directory for debugging
        echo "Process working directory: \$PWD"
        echo "GATK temp directory: ${tmp_dir}"

        gatk --java-options "-Xmx${java_mem} -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}" MarkDuplicates \
            -I ${bam} \
            -O ${outfile_bam} \
            -M ${outfile_metrics} \
            --CREATE_INDEX true \
            --TMP_DIR ${tmp_dir} \
            --ASSUME_SORTED true \
            --MAX_RECORDS_IN_RAM 2000000

        # Clean up temporary directory
        rm -rf ${tmp_dir}
        """
}