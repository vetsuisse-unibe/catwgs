#!/bin/bash
#SBATCH --job-name=catwgs_v2_pipeline
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=172:00:00
#SBATCH --output=logs/nextflow_v2_%j.out
#SBATCH --error=logs/nextflow_v2_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=$USER@unibe.ch

# Create required directories if they don't exist
mkdir -p logs reports

# Load required modules
module load Java/17.0.6

# Set Nextflow configuration
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_HOME=$HOME/.nextflow

# Change to pipeline directory
cd /data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2

# Run the pipeline
echo "Starting Nextflow CATWGS v2 pipeline at $(date)"
echo "Sample sheet: $(realpath assets/sampleSheet_Sep2025.txt)"
echo "Number of samples: $(($(wc -l < assets/sampleSheet_Sep2025.txt) - 1))"

# Entry point options (see workflows/main_v2.nf):
# --entry_point markdup           : FASTQ -> dedup BAMs
# --entry_point start             : FASTQ -> dedup BAMs -> per-sample gVCFs
# --entry_point haplotypecaller   : existing dedup BAMs -> per-sample gVCFs
# --entry_point cohortmap         : all gVCFs -> joint genotyping -> filtered, annotated VCF
# --entry_point variantprocessing : existing regional VCFs -> filtered, annotated VCF

/data/users/vjaganna/software/nextflow run main.nf \
    -profile unibe \
    --entry_point start \
    -with-report reports/nextflow_v2_report_$(date +%Y%m%d_%H%M%S).html \
    -with-timeline reports/nextflow_v2_timeline_$(date +%Y%m%d_%H%M%S).html \
    -with-trace reports/nextflow_v2_trace_$(date +%Y%m%d_%H%M%S).txt \
    -resume

echo "Pipeline completed at $(date)"
