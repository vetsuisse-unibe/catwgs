// Define genotypeGVCFs process
process genotypeGVCFs {
    tag "${interval}"
    scratch true  // Use scratch space for the entire process
    publishDir params.vcfFolder, mode: 'copy'

    // Retry on transient failures, including external termination (null exit
    // status) from SLURM node/staging issues, so one bad node doesn't abort
    // the full scatter.
    errorStrategy { (task.exitStatus == null || task.exitStatus in [143,137,104,134,139,1,255]) ? 'retry' : 'finish' }
    maxRetries 3
    maxErrors -1

    input:
    tuple val(interval), path(dbfile)
    path(ref)
    path(fai)
    path(dict)

    output:
    tuple path(vcfFile), path("${vcfFile}.tbi"), emit: vcf

    script:
    fileIntervalString = "${interval}".replaceAll(':','_')
    vcfFile = "cohort_1543_${fileIntervalString}.vcf.gz"
    """
    # Create temporary directory
    mkdir -p \$SCRATCH/tmp_vj
    
    # Print current working directory for debugging
    echo "Process working directory: \$PWD"
    echo "GenotypeGVCFs temp directory: \$SCRATCH/tmp_vj"
    
    gatk --java-options "-Xmx5g -Xms5g" GenotypeGVCFs \\
        --tmp-dir \$SCRATCH/tmp_vj \\
        -R ${ref} \\
        -O ${vcfFile} \\
        -V gendb://${dbfile}
        
    # Index the VCF file
    gatk --java-options "-Xmx5g" IndexFeatureFile -I ${vcfFile}
    """
}
