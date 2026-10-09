#!/bin/bash
#SBATCH --job-name=catwgs_pipeline
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=14-00:00:00
#SBATCH --output=logs/nextflow_%j.out
#SBATCH --error=logs/nextflow_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Vidhya.Jagannathan@unibe.ch

# Create logs directory if it doesn't exist
mkdir -p logs

# Load required modules
module load Java/17.0.6

# Set Nextflow configuration
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_HOME=$HOME/.nextflow

# Change to pipeline directory
cd /data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2

# Run the pipeline (variant processing: re-gather all per-region VCFs incl. MT +
# NW_ scaffolds, then select -> filter -> annotate -> merge over the full cohort)
echo "Starting Nextflow CATWGS pipeline at $(date)"
echo "Entry point: variantprocessing (gather from per-region VCFs in vcfFolder)"
echo "vcf folder: /data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/vcf/"
echo "Number of per-region VCFs: $(ls /data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/vcf/*.vcf 2>/dev/null | wc -l)"

/data/users/vjaganna/software/nextflow run main.nf \
    -resume \
    -profile unibe \
    --entry_point variantprocessing \
    --regional_vcf_glob '*.vcf' \
    -with-report reports/nextflow_report_$(date +%Y%m%d_%H%M%S).html \
    -with-timeline reports/nextflow_timeline_$(date +%Y%m%d_%H%M%S).html \
    -with-trace reports/nextflow_trace_$(date +%Y%m%d_%H%M%S).txt

echo "Pipeline completed at $(date)"
