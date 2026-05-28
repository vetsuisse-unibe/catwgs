// Define genomicsDB process
process genomicsDB {
    tag "${interval}"
    scratch true  // Use scratch space for the entire process
    
    // Retry on known transient codes AND on external termination (null exit
    // status), which is how SLURM node/staging failures surface. Without the
    // null case a single bad node aborts the whole multi-thousand-task run.
    errorStrategy { (task.exitStatus == null || task.exitStatus in [143,137,104,134,139,1,255]) ? 'retry' : 'finish' }
    maxRetries 3

    // Allow many transient retries across the full scatter before giving up.
    maxErrors -1
    time { 24.hour * task.attempt }

    input:
    tuple path(samples), val(interval)

    output:
    tuple val(interval), path(dbfile), emit: genomicsdb_out

    script:
    fileIntervalString = "${interval}".replaceAll(':','_')
    dbfile = "chunkDB_${fileIntervalString}"
    
    // Calculate the number of samples to determine appropriate memory settings
    """
    # Create temporary directory with unique name to avoid conflicts
    TEMP_DIR="\$SCRATCH/tmp_genomicsdb_${fileIntervalString}_\$RANDOM"
    mkdir -p \$TEMP_DIR
    
    # Print current working directory and environment for debugging
    echo "Process working directory: \$PWD"
    echo "GenomicsDB temp directory: \$TEMP_DIR"
    echo "Available memory: \$(free -h)"
    echo "Disk space: \$(df -h \$TEMP_DIR)"
    
    # Count number of samples for logging
    SAMPLE_COUNT=\$(wc -l < ${samples})
    echo "Processing \$SAMPLE_COUNT samples for interval ${interval}"
    
    # Clean up any existing workspace with the same name
    if [ -d "${dbfile}" ]; then
        echo "Removing existing workspace directory: ${dbfile}"
        rm -rf "${dbfile}"
    fi
    
    # Set more conservative Java options with better garbage collection
    # Use smaller batch size to reduce memory pressure
    echo "Starting GenomicsDBImport at \$(date)"
    gatk --java-options "-Xmx20g -Xms8g -XX:+UseG1GC -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+DisableExplicitGC -Djava.io.tmpdir=\$TEMP_DIR" \\
        GenomicsDBImport \\
        --tmp-dir \$TEMP_DIR \\
        --L ${interval} \\
        --sample-name-map ${samples} \\
        --batch-size 50 \\
        --genomicsdb-workspace-path ${dbfile} \\
        --genomicsdb-shared-posixfs-optimizations
    
    # Check the exit status
    GATK_EXIT=\$?
    echo "GenomicsDBImport completed with exit code \$GATK_EXIT at \$(date)"
    
    # Clean up temp directory
    echo "Cleaning up temporary directory: \$TEMP_DIR"
    rm -rf \$TEMP_DIR
    
    # Verify the output exists
    if [ -d "${dbfile}" ]; then
        echo "Successfully created GenomicsDB workspace: ${dbfile}"
        ls -la "${dbfile}"
    else
        echo "ERROR: Failed to create GenomicsDB workspace: ${dbfile}"
        exit 1
    fi
    
    exit \$GATK_EXIT
    """
}