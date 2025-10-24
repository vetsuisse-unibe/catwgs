#!/bin/bash
# Slurm options
#SBATCH --mail-user=vidhya.jagannathan@unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --job-name="eVCFs"
#SBATCH --chdir=.
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --partition=pibu_el8
#SBATCH --array=1-2478   
#SBATCH -o logs/region_%A_%a.out
#SBATCH -e logs/region_%A_%a.err

module load HTSlib/1.12-GCC-10.3.0
module load zstd/1.4.9-GCCcore-10.3.0

# Get the region for this array task
region=$(sed -n "${SLURM_ARRAY_TASK_ID}p" regions.wholeGenome.txt)
outfile="extract_results/$(echo $region | tr ':' '_' | tr '-' '_').zst"

# Create output directory if it doesn't exist
mkdir -p extract_results logs

# Process the region
tabix /data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/final_vcf/cohort_124.var.flt.ann.vcf.gz ${region} |
    /data/users/vjaganna/software/snpEff/scripts/vcfEffOnePerLine.pl |
    java -Xmx16g -jar /data/users/vjaganna/software/snpEff/SnpSift.jar extractFields - \
    CHROM POS REF ALT FILTER "GEN[*].GT" "ANN[*].ALLELE" "ANN[*].EFFECT" \
    "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" \
    "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \
    "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" \
    "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" \
    "ANN[*].ERRORS" | sed 's/|/\//g' | zstd > ${outfile}