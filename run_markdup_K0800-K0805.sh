#!/bin/bash
#SBATCH --job-name=markdup_K0800-K0805
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/markdup_K0800-K0805_%j.out
#SBATCH --error=logs/markdup_K0800-K0805_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=$USER@unibe.ch

# Create logs directory if it doesn't exist
mkdir -p logs

# Load required modules
module load Java/17.0.6

# Set Nextflow configuration
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_HOME=$HOME/.nextflow

# Change to pipeline directory
cd /data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2

# Run the pipeline
echo "Starting MarkDuplicates pipeline for K0800, K0801, K0805 at $(date)"
echo "Sample sheet: $(realpath assets/sampleSheet_March2025_K0800-K0805.txt)"
echo "Number of samples: 3"

# Clean Nextflow cache to force fresh run
rm -rf .nextflow .nextflow.*

/data/users/vjaganna/software/nextflow run main.nf \
    -profile unibe \
    --entry_point markdup \
    --samples assets/sampleSheet_March2025_K0800-K0805.txt \
    -with-report reports/markdup_K0800-K0805_report_$(date +%Y%m%d_%H%M%S).html \
    -with-timeline reports/markdup_K0800-K0805_timeline_$(date +%Y%m%d_%H%M%S).html \
    -with-trace reports/markdup_K0800-K0805_trace_$(date +%Y%m%d_%H%M%S).txt

echo "Pipeline completed at $(date)"
