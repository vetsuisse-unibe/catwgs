#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process runBWA {
    tag "${fastQbasename}|${sampleID}|${libraryID}|${rgID}"
    scratch true  // Use scratch space for the entire process
    

    input:
      tuple val(fastQbasename), path(fastqR1), path(fastqR2),val(sampleID), val(libraryID), val(rgID), val(platform), val(model), val(run_date),val(centre),val(pu)
      tuple val(assembly), val(REF)
    output:
      tuple val(fastQbasename), val(sampleID), path('*aligned.bam'), path('*aligned.bam.bai'), emit: mappedBams

    script:
      this_chunk = fastqR1.getName().split("-")[0].substring(0,4)
      outfileSam = fastQbasename + ".aligned.sam"
      outfile = fastQbasename + ".aligned.bam"
      outfile_index = outfile + ".bai"
      // Use a work directory-based temporary directory
      tmp_dir = "\$PWD/tmp_bwa_${fastQbasename}"

      // Calculate memory per thread for samtools sort
      def total_mem = task.memory.toGiga()
      def mem_per_thread = (total_mem / task.cpus).intValue()
      def sort_mem = mem_per_thread > 1 ? "${mem_per_thread - 1}G" : "1G"
      
      """
      # Create temporary directory
      mkdir -p ${tmp_dir}
      
      # Print current working directory for debugging
      echo "Process working directory: \$PWD"
      echo "BWA temp directory: ${tmp_dir}"
      
      echo "Received REF path: ${REF}"
      bwa-mem2 mem -K 100000000 -Y -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform}\\tSM:${sampleID}\\tLB:${libraryID}\\tDS:${assembly}\\tCN:UNIBE" -t ${task.cpus} ${REF} $fastqR1 $fastqR2 >$outfileSam
      samtools sort -m ${sort_mem} -@ ${task.cpus} -T "${tmp_dir}/${fastQbasename}_tmp" -o $outfile $outfileSam 
      samtools index $outfile
      
      # Clean up temporary directory
      rm -rf ${tmp_dir}
      """
}
