#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process createCohortMap {
    tag "Creating cohort map for samples"
    publishDir "${params.cohortMapFolder}/assets", mode: 'copy'

    // Always re-scan the gVCF folder so newly added gVCFs are picked up under
    // -resume. The only input is the folder path string, so without this the task
    // would cache on the path and reuse a stale map. The map is cheap to rebuild;
    // downstream genomicsDB stays cache-gated on the map's contents.
    cache false

    input:
    val gvcf_folder

    output:
    path "cohort_*.sample_map", emit: cohortMapFile

    script:
    """
    # Find all gVCF files and create the cohort map (maxdepth 1 to avoid subdirectories)
    SAMPLE_COUNT=\$(find ${gvcf_folder} -maxdepth 1 -name "*.g.vcf.gz" | wc -l)
    COHORT_FILE="cohort_\${SAMPLE_COUNT}.sample_map"

    # sort for deterministic output and write in one shot (truncate, never append)
    find ${gvcf_folder} -maxdepth 1 -name "*.g.vcf.gz" | sort | while read vcf; do
        sample=\$(basename \$vcf .g.vcf.gz)
        echo -e "\${sample}\t\${vcf}"
    done > \${COHORT_FILE}
    """
}
