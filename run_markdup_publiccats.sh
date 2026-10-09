#!/bin/bash
#SBATCH --job-name=markdup_publiccats
#SBATCH --partition=pibu_el8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4-00:00:00
#SBATCH --output=logs/markdup_publiccats_%j.out
#SBATCH --error=logs/markdup_publiccats_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=Vidhya.Jagannathan@unibe.ch

mkdir -p logs reports

module load Java/17.0.6
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_HOME=$HOME/.nextflow

cd /data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2

DEDUP=/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dedup_bams/publiccats
mkdir -p "$DEDUP"

echo "Starting markDuplicates-only pipeline (fastp -> BWA -> markDuplicates) at $(date)"
echo "Samples: assets/sampleSheet_publiccats.txt (Boo=SRR25054056, Onyx=SRR25054055)"
echo "Dedup BAM output: $DEDUP"

/data/users/vjaganna/software/nextflow run main.nf \
    -profile unibe \
    --entry_point markdup \
    --samples assets/sampleSheet_publiccats.txt \
    --dedup "$DEDUP" \
    -with-report reports/markdup_publiccats_report_$(date +%Y%m%d_%H%M%S).html \
    -with-timeline reports/markdup_publiccats_timeline_$(date +%Y%m%d_%H%M%S).html \
    -with-trace reports/markdup_publiccats_trace_$(date +%Y%m%d_%H%M%S).txt \
    -resume

echo "Pipeline completed at $(date)"
