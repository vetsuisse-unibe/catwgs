#!/bin/bash
# Master submission script for variant extraction pipeline
# Submits extraction array job followed by merge job

echo "=========================================="
echo "Submitting Variant Extraction Pipeline"
echo "Cohort: cohort_136 (136 samples)"
echo "Regions: 2478 genomic regions"
echo "Timestamp: $(date)"
echo "=========================================="

# Step 1: Create header file first
echo ""
echo "Step 1: Creating header file..."
bash create_header.sh

if [ ! -f "head.txt" ]; then
    echo "ERROR: head.txt was not created successfully!"
    exit 1
fi

echo "  - head.txt created successfully"
echo "  - Number of fields: $(head -1 head.txt | awk '{print NF}')"

# Step 2: Submit extraction array job (2478 regions in parallel)
# Optional arg $1 = an upstream job ID (e.g. VCF indexing) to wait on (afterok).
echo ""
echo "Step 2: Submitting extraction array job..."
DEP_ARG=""
if [ -n "${1:-}" ]; then
    DEP_ARG="--dependency=afterok:$1"
    echo "  - extraction will wait for upstream job $1 (afterok)"
fi
JOB1=$(sbatch --parsable $DEP_ARG extract_regions.sh)
echo "  - extract_regions.sh submitted: Job ID $JOB1 (2478 array tasks)"

# Step 3: Submit merge job (depends on all extraction tasks completing)
echo ""
echo "Step 3: Submitting merge job (depends on extraction)..."
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 merge_extract_regions.sh)
echo "  - merge_extract_regions.sh submitted: Job ID $JOB2 (depends on $JOB1)"

echo ""
echo "=========================================="
echo "All jobs submitted successfully!"
echo "=========================================="
echo ""
echo "Job Dependency Chain:"
echo "  $JOB1 (extract_regions - 2478 parallel tasks) -> $JOB2 (merge_extract_regions)"
echo ""
echo "Monitor jobs with:"
echo "  - All jobs: squeue -u $USER"
echo "  - Array job: squeue -j $JOB1"
echo "  - Merge job: squeue -j $JOB2"
echo ""
echo "Check progress:"
echo "  - Completed regions: ls extract_results/*.zst | wc -l"
echo "  - Expected regions: 2478"
echo ""
echo "Output will be in: final_output.txt"
echo ""
