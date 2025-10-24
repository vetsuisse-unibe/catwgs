#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process createCohortMap {
    tag "Creating cohort map for samples"
    publishDir "${params.cohortMapFolder}", mode: 'copy'

    input:
    val gvcf_folder

    output:
    path("cohort.sample_map"), emit: cohortMapFile

    script:
    """
    #!/usr/bin/env python3
    
    import os 
    import glob

    # Create the output file
    with open("cohort.sample_map", "w") as outfile:
        # Find all gVCF files in the specified directory
        for gvcf in glob.glob("${gvcf_folder}/*.g.vcf.gz"):
            # Extract the sample name from the filename
            sample = os.path.basename(gvcf).split('.')[0]
            # Write the sample and file path to the cohort map
            outfile.write(f"{sample}\\t{gvcf}\\n")
    """
}