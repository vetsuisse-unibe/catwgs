#!/bin/bash
# Script to create head.txt for variant extraction
# This creates a header file with sample names and annotation fields

set -e

VCF_FILE="final_vcf/cohort_124.var.flt.ann.vcf.gz"
OUTPUT_FILE="head.txt"

echo "Creating header file: ${OUTPUT_FILE}"

# Load required modules
module load BCFtools

# Extract header line and process it
bcftools view --header ${VCF_FILE} | grep "^#CHROM" | \
  awk 'BEGIN {OFS="\t"} {
    # Print CHROM, POS, REF, ALT, FILTER
    printf "%s\t%s\t%s\t%s\t%s", $1, $2, $4, $5, $7;
    
    # Print all sample names (starting from column 10)
    for (i=10; i<=NF; i++) {
      printf "\t%s", $i;
    }
    
    # Add annotation fields at the end
    printf "\tALLELE\tEFFECT\tIMPACT\tGENE\tGENEID\tFEATURE\tFEATUREID\tBIOTYPE\tRANK\tHGVS_C\tHGVS_P\tCDNA_POS\tCDNA_LEN\tCDS_POS\tCDS_LEN\tAA_POS\tAA_LEN\tDISTANCE\tERRORS\n";
  }' > ${OUTPUT_FILE}

echo "Header file created successfully!"
echo "Number of fields: $(head -1 ${OUTPUT_FILE} | awk '{print NF}')"
echo "First few fields: $(head -1 ${OUTPUT_FILE} | cut -f1-10)"
