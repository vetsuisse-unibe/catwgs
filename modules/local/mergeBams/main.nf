#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process mergeBams{
      debug true
  tag "${sampleID}"
  scratch true  // Use scratch space for the entire process

  input:
  tuple val(fastQbasename),val(sampleID), path(bams),path(bais) 
  output:
  tuple val(sampleID), path("${outfile}"),path("${outfile_index}"),  emit: inputMarkDuplicates
   
  script:
  inputfile = bams.join(" I=")
  outfile = sampleID + ".merged.bam"
  outfile_index = outfile + ".bai"
  // Use a work directory-based temporary directory
  tmp_dir = "\$PWD/tmp_merge_${sampleID}"
  
  """
  # Create temporary directory
  mkdir -p ${tmp_dir}
  
  # Print current working directory for debugging
  echo "Process working directory: \$PWD"
  echo "Merge BAMs temp directory: ${tmp_dir}"
  
  java -jar \$EBROOTPICARD/picard.jar MergeSamFiles \
    I=$inputfile \
    O=$outfile \
    SORT_ORDER=coordinate \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=${tmp_dir}
    
  samtools index $outfile 
  
  # Clean up temporary directory
  rm -rf ${tmp_dir}
  """
}
