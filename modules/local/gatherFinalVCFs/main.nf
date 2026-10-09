// Define gatherFinalVCFs process
process gatherFinalVCFs {
    tag "Gathering final VCFs"
    scratch true  // Use scratch space for the entire process
    publishDir params.finalVcfFolder, mode: 'copy'

    input:
    path(vcfs)
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
    # Chunks may be plain .vcf or bgzipped .vcf.gz
    shopt -s nullglob
    chunks=( *.vcf *.vcf.gz )
    echo "Number of VCF files: \${#chunks[@]}"

    # Extract chromosome order from the reference dictionary (preserves dict order)
    grep -oP '(?<=SN:)[^\\t]*' ${dict} > chrom_order.txt

    # Build the sorted gather list from each chunk's ACTUAL first-record contig and
    # position, read from the VCF content itself. This is robust to every contig
    # naming style (chr*, NW_* scaffolds, MT) instead of regex-parsing the filename,
    # which silently dropped MT (no 'chr'/'NW_' prefix) and any non-matching contig.
    # Empty chunks (header only) carry no variants and are skipped.
    : > sorted_keys.tsv
    for vcf in "\${chunks[@]}"; do
        first=\$(gzip -cdf "\$vcf" 2>/dev/null | grep -m1 -v '^#' || true)
        if [ -z "\$first" ]; then
            echo "Skipping empty chunk (no records): \$vcf" >&2
            continue
        fi
        chrom=\$(printf '%s' "\$first" | cut -f1)
        pos=\$(printf '%s' "\$first" | cut -f2)

        # Index of this contig in the reference dictionary order
        idx=\$(grep -nxF "\$chrom" chrom_order.txt | head -1 | cut -d':' -f1 || true)
        if [ -z "\$idx" ]; then
            # Contig not in dictionary: place at the end (should not happen)
            idx=999999
        fi

        printf "%06d\\t%012d\\t%s\\n" "\$idx" "\$pos" "\$vcf" >> sorted_keys.tsv
    done
    sort -k1,1n -k2,2n sorted_keys.tsv | cut -f3 > sorted_vcf.list

    if [ ! -s sorted_vcf.list ]; then
        echo "ERROR: no non-empty VCF chunks to gather" >&2
        exit 1
    fi

    # Two chunks starting at the same contig/position means the input folder holds
    # duplicate chunks (e.g. old .vcf next to new .vcf.gz); GatherVcfs would fail or
    # duplicate records.
    dups=\$(cut -f1,2 sorted_keys.tsv | sort | uniq -d)
    if [ -n "\$dups" ]; then
        echo "ERROR: multiple VCF chunks start at the same position:" >&2
        grep -F -f <(printf '%s\\n' "\$dups") sorted_keys.tsv >&2
        exit 1
    fi

    # Print the sorted list for debugging
    echo "Sorted VCF list (first 10):"
    head sorted_vcf.list
    echo "Total VCFs to merge: \$(wc -l < sorted_vcf.list)"

    # Gather all VCFs into a single file
    gatk --java-options "-Xmx30G -Djava.io.tmpdir=\$SCRATCH/tmp_gather" GatherVcfs \\
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
