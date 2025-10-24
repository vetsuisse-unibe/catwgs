nextflow.enable.dsl=2

process runFastp {
    tag "${fastQbasename}"
    scratch true  // Use scratch space for the entire process

    input:
      tuple  val(fastQbasename),path(fastqR1), path(fastqR2)

    output:
      tuple val(fastQbasename), path("*${left}"), path("*${right}"), emit: trimmedReads
      tuple val(fastQbasename), path(json), path(html), emit: trimmingReports

    script:
      left = file(fastqR1).getBaseName() + ".cleaned.fastq.gz"
      right = file(fastqR2).getBaseName() + ".cleaned.fastq.gz"
      json = fastQbasename + ".fastp.json"
      html = fastQbasename + ".fastp.html"
      // Use a work directory-based temporary directory instead of params.scratch_dir
      tmp_dir = "\$PWD/tmp_fastp_${fastQbasename}"
      """
      # Create temporary directory
      mkdir -p ${tmp_dir}
      
      # Print current working directory for debugging
      echo "Process working directory: \$PWD"
      echo "Fastp temp directory: ${tmp_dir}"
      
      # Run fastp with temp directory
      TMPDIR=${tmp_dir} fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right -w ${task.cpus} -j $json -h $html
      
      # Clean up temporary directory
      rm -rf ${tmp_dir}
      """
}
