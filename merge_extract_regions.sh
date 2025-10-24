#!/bin/bash
# Slurm options
#SBATCH --mail-user=vidhya.jagannathan@unibe.ch
#SBATCH --mail-type=fail,end
#SBATCH --job-name="mergeTxtVCFs"
#SBATCH --chdir=.
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=pibu_el8
#SBATCH -o logs/merge_%j.out
#SBATCH -e logs/merge_%j.err

module load HTSlib/1.12-GCC-10.3.0
module load zstd/1.4.9-GCCcore-10.3.0


# Combine all results in the original order
cat head.txt > final_output.txt

# Combine all results in order, skipping header lines
while read region; do
    filename="extract_results/$(echo $region | tr ':' '_' | tr '-' '_').zst"
    if [ -f "$filename" ]; then
        zstd -dc "$filename" | grep -v "^CHROM" >> final_output.txt
    else
        echo "Warning: Missing file for region $region" >&2
    fi
done < /data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/regions.wholeGenome.txt

# # Compress final output
# zstd final_output.txt -o cohort.1582.flt.var.ann.zst
# rm final_output.txt header.txt