process createCohortMapFile {
    tag "Gathering gVCFs from folder"
    
    input:
    path gvcf_folder

    output:
    path 'all_gvcfs.txt', emit: gvcf_list

    script:
    """
    #!/usr/bin/env python3
    import os
    import glob

    with open('all_gvcfs.txt', 'w') as outfile:
        for gzVCF in glob.glob('${gvcf_folder}/*.gz'):
            file = os.path.basename(gzVCF)
            tmp = file.split(".")
            outfile.write(f"{tmp[0]}\\t{gzVCF}\\n")
    """
}