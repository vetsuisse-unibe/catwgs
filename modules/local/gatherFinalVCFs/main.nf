// Define gatherFinalVCFs process
process gatherFinalVCFs {
    tag "Gathering final VCFs"
    scratch true  // Use scratch space for the entire process
    publishDir params.finalVcfFolder, mode: 'copy'

    input:
    path(vcf_list)
    path(ref)
    path(fai)
    path(dict)

    output:
    tuple path("cohort_final.vcf.gz"), path("cohort_final.vcf.gz.tbi"), emit: final_vcf

    script:
    """
    # Create temporary directory
    mkdir -p \$SCRATCH/tmp_gather
    
    # Print current working directory for debugging
    echo "Process working directory: \$PWD"
    echo "GatherVcfs temp directory: \$SCRATCH/tmp_gather"
    
    # Create a properly sorted VCF list based on chromosome order in the reference
    # Extract chromosome order from the reference dictionary
    grep -oP '(?<=SN:).*?(?=\\t)' ${dict} > chrom_order.txt
    
    # Create a sorted list based on chromosome and position
    cat ${vcf_list} | while read vcf; do
        # Extract chromosome/scaffold name and position from filename
        # Match both chr* and NW_* (scaffolds)
        chrom=\$(echo \$vcf | grep -oP '(chr[^_]*|NW_[0-9]+\\.[0-9]+)')
        pos=\$(echo \$vcf | grep -oP '\\d+-\\d+')
        
        # Find the index of this chromosome in the reference order
        idx=\$(grep -n "^\$chrom\$" chrom_order.txt | cut -d':' -f1)
        if [ -z "\$idx" ]; then
            # If not found, assign a large number to put it at the end
            idx=999999
        fi
        
        # Output with sortable prefix
        printf "%06d\\t%s\\t%s\\n" "\$idx" "\$pos" "\$vcf"
    done | sort -k1,1n -k2,2V | cut -f3 > sorted_vcf.list
    
    # Print the sorted list for debugging
    echo "Sorted VCF list:"
    cat sorted_vcf.list
    
    # Gather all VCFs into a single file
    gatk --java-options "-Xmx6G -Djava.io.tmpdir=\$SCRATCH/tmp_gather" GatherVcfs \\
        -R ${ref} \\
        -I sorted_vcf.list \\
        -O cohort_final.vcf.gz
    
    # Index the final VCF
    gatk --java-options "-Xmx5g" IndexFeatureFile \\
        -I cohort_final.vcf.gz
    
    # Clean up temp directory
    rm -rf \$SCRATCH/tmp_gather
    """
}
