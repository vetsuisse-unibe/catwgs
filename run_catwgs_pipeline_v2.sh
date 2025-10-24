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

# Entry point options for main_v2.nf:
# --entry_point fastp         : Start from FastP (default - full pipeline)
# --entry_point markduplicates: Start from mark duplicates (requires --input_bams)
# --entry_point gathervcfs    : Start from gather VCFs (requires --input_gvcfs)

/data/users/vjaganna/software/nextflow run workflows/main_v2.nf \
    -profile unibe \
    --entry_point fastp \
    -with-report reports/nextflow_v2_report_$(date +%Y%m%d_%H%M%S).html \
    -with-timeline reports/nextflow_v2_timeline_$(date +%Y%m%d_%H%M%S).html \
    -with-trace reports/nextflow_v2_trace_$(date +%Y%m%d_%H%M%S).txt \
    -resume

echo "Pipeline completed at $(date)"
