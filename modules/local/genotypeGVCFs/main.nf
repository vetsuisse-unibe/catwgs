// Define genotypeGVCFs process
process genotypeGVCFs {
    tag "${interval}"
    scratch true  // Use scratch space for the entire process
    publishDir params.vcfFolder, mode: 'copy'

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
