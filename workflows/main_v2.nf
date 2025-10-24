#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Add entry point parameters
params.tag = 'default' // Initialize the tag parameter with a default value
params.entry_point = 'fastp'  // Default start point
params.input_bams = null      // For starting from markDuplicates
params.input_gvcfs = null     // For starting from gatherVCFs

// Include processes
include { runFastp } from '../modules/local/fastp/main'
include { runBWA } from '../modules/local/bwa-mem/main'
include { mergeBams } from '../modules/local/mergeBams/main'
include { markDuplicates } from '../modules/local/markDuplicates/main'
include { haplotypeCaller } from '../modules/local/hc/main'
include { gatherVCFs } from '../modules/local/gatherVCFs/main'
include { indexgVCF } from '../modules/local/indexgVCF/main'
include { createCohortMap } from '../modules/local/createCohortMap/main'
include { genomicsDB } from '../modules/local/genomicsDB/main'
include { genotypeGVCFs } from '../modules/local/genotypeGVCFs/main'
include { gatherFinalVCFs } from '../modules/local/gatherFinalVCFs/main'
include { selectSNP } from '../modules/local/selectSNP/main'
include { selectNonSNP } from '../modules/local/selectNonSNP/main'
include { filterSNPs } from '../modules/local/filterSNPs/main'
include { filterNonSNPs } from '../modules/local/filterNonSNPs/main'
include { annotateSNPs } from '../modules/local/annotateSNPs/main'
include { annotateNonSNPs } from '../modules/local/annotateNonSNPs/main'
include { mergeAnnotatedVCFs } from '../modules/local/mergeAnnotatedVCFs/main'

workflow {
    // Entry-point specific validation
    if (params.entry_point == 'fastp') {
        if (!params.samples) {
            error "Samples file not specified with params.samples!"
        }
        if (!params.intervals_folder) {
            error "Intervals folder not specified with params.intervals_folder!"
        }
        if (!file(params.ref).exists()) {
            error "Reference file not found: ${params.ref}"
        }
        if (!file(params.fai).exists()) {
            error "Reference index file not found: ${params.fai}"
        }
        if (!file(params.dict).exists()) {
            error "Reference dictionary file not found: ${params.dict}"
        }
    } else if (params.entry_point == 'markduplicates') {
        if (!params.input_bams) {
            error "When starting from markduplicates, you must provide input BAMs with params.input_bams"
        }
        if (!params.intervals_folder) {
            error "Intervals folder not specified with params.intervals_folder!"
        }
        if (!file(params.ref).exists()) {
            error "Reference file not found: ${params.ref}"
        }
        if (!file(params.fai).exists()) {
            error "Reference index file not found: ${params.fai}"
        }
        if (!file(params.dict).exists()) {
            error "Reference dictionary file not found: ${params.dict}"
        }
    } else if (params.entry_point == 'gathervcfs') {
        if (!params.input_gvcfs) {
            error "When starting from gathervcfs, you must provide input gVCFs with params.input_gvcfs"
        }
    } else {
        error "Invalid entry_point: ${params.entry_point}. Valid options are: fastp, markduplicates, gathervcfs"
    }
    
    // Common validation for all entry points
    if (!params.regionsFile) {
        error "Regions file not specified with params.regionsFile!"
    }
    if (!params.gVCF_folder) {
        error "gVCF folder not specified with params.gVCF_folder!"
    }
    
    // Define interval and regions channels
    if (params.entry_point in ['fastp', 'markduplicates']) {
        // Create interval channel for HaplotypeCaller
        def correctIntervalPath = "/data/references/Felis_catus/NCBI/F.catus_Fca126_mat1.0/Sequence/interval-files-folder-48"
        log.info "Intervals folder path: ${correctIntervalPath}"
        def intervalsFolder = file(correctIntervalPath)
        log.info "Intervals folder exists: ${intervalsFolder.exists()}"
        log.info "Intervals folder absolute path: ${intervalsFolder.toAbsolutePath()}"
        
        interval_ch = Channel
            .fromPath("${correctIntervalPath}/*-scattered.interval_list")
            .ifEmpty { error "No interval files found in ${correctIntervalPath}" }
    }
    
    // Create regions channel for all entry points
    log.info "Regions file path: ${params.regionsFile}"
    def regionsFile = file(params.regionsFile)
    log.info "Regions file exists: ${regionsFile.exists()}"
    log.info "Regions file absolute path: ${regionsFile.toAbsolutePath()}"
    
    regions_ch = Channel
        .fromPath(params.regionsFile)
        .splitCsv(header: false)
        .ifEmpty { error "No regions found in ${params.regionsFile}" }

    // Create main input channel based on entry point
    if (params.entry_point == 'fastp') {
        // Debug sample sheet path
        log.info "Sample sheet path: ${params.samples}"
        def sampleFile = file(params.samples)
        log.info "Sample file exists: ${sampleFile.exists()}"
        log.info "Sample file absolute path: ${sampleFile.toAbsolutePath()}"
        
        Channel
            .fromPath(params.samples)
            .ifEmpty { error "Cannot find samples file: ${params.samples}" }
            .splitCsv(header: true, sep: ';', strip: true)
            .map { row -> 
                def r1 = file(row.R1)
                def r2 = file(row.R2)
                if (!r1.exists()) error "R1 file not found: ${row.R1}"
                if (!r2.exists()) error "R2 file not found: ${row.R2}"
                [
                    fastQbasename: row.fastQbasename,
                    sampleID: row.sampleID,
                    libraryID: row.libraryID,
                    rgID: row.rgID,
                    platform: row.platform,
                    model: row.model,
                    centre: row.center,
                    date: row.run_date,
                    pu: row.pu,
                    R1: r1,
                    R2: r2
                ]
            }
            .set { base_channel }

        // Create fastqc input channel
        fastqc_input = base_channel
            .map { row -> tuple(
                row.fastQbasename,
                row.R1,
                row.R2
            )}

        // Create BWA input channel
        bwa_mem_input = base_channel
            .map { row ->
                tuple(
                    row.fastQbasename,
                    row.R1,
                    row.R2,
                    row.sampleID,
                    row.libraryID,
                    row.rgID,
                    row.platform,
                    row.model,
                    row.date,
                    row.centre,
                    row.pu
                )
            }

        // Run FastP
        trimmedReads = runFastp(fastqc_input).trimmedReads

        // Join trimmed reads with metadata for BWA
        combined_for_bwa_mem = trimmedReads.join(
            bwa_mem_input.map { row ->
                tuple(
                    row[0], // fastQbasename as the join key
                    row[3], // sampleID
                    row[4], // libraryID
                    row[5], // rgID
                    row[6], // platform
                    row[7], // model
                    row[8], // date
                    row[9], // centre
                    row[10] // pu
                )
            }
        ).map { fastQbasename, r1, r2, sampleID, libraryID, rgID, platform, model, date, centre, pu ->
            tuple(fastQbasename, r1, r2, sampleID, libraryID, rgID, platform, model, date, centre, pu)
        }

        // Run BWA
        mappedBams = runBWA(
            combined_for_bwa_mem, 
            tuple(params.assembly, params.ref)
        ).mappedBams

        // Group BAMs by sample ID
        mappedBams_by_sampleID = mappedBams
            .map { fastQbasename, sampleId, bam, bai -> 
                tuple(fastQbasename, sampleId, bam, bai) 
            }
            .groupTuple(by: 1)

        // Split into single and multiple lanes
        channels = mappedBams_by_sampleID
            .branch {
                singleLane: it[0].size() == 1
                multipleLanes: it[0].size() > 1
            }

        singleLaneSamples = channels.singleLane
            .map { fastQbasenames, sampleId, bams, bais -> 
                tuple(sampleId, bams[0], bais[0])
            }

        // Merge multiple lane BAMs
        mergedBams = mergeBams(channels.multipleLanes).inputMarkDuplicates

        // Combine single and merged BAMs
        allBams = singleLaneSamples
            .mix(mergedBams)
            .ifEmpty { error "No BAM files to process" }
            // Verify that each tuple has at least two elements
            .map { it -> 
                if (it != null && it.size() < 2) {
                    log.warn "Warning: Tuple for sample ${it ? it[0] : 'unknown'} has fewer than 2 elements in allBams"
                }
                return it
            }

    } else if (params.entry_point == 'markduplicates') {
        // Create channel from input BAMs
        allBams = Channel
            .fromPath(params.input_bams)
            .ifEmpty { error "Cannot find input BAMs file: ${params.input_bams}" }
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def bam = file(row.bam)
                def bai = file("${row.bam}.bai")
                if (!bam.exists()) error "BAM file not found: ${row.bam}"
                if (!bai.exists()) error "BAI file not found: ${row.bam}.bai"
                tuple(row.sample_id, bam, bai)
            }
    } else if (params.entry_point == 'gathervcfs') {
        // Create channel from input gVCFs
        gatheredgVCFs = Channel
            .fromPath(params.input_gvcfs)
            .ifEmpty { error "Cannot find input gVCFs file: ${params.input_gvcfs}" }
            .map { gvcf ->
                def sample_id = gvcf.simpleName
                // Ensure tuple has at least two elements
                tuple(sample_id, gvcf)
            }
            // Verify that each tuple has at least two elements
            .map { it -> 
                if (it != null && it.size() < 2) {
                    log.warn "Warning: Tuple for sample ${it ? it[0] : 'unknown'} has fewer than 2 elements"
                }
                return it
            }
    }

    // Continue workflow based on entry point
    if (params.entry_point in ['fastp', 'markduplicates']) {
        // Mark duplicates
        duplicateMarkedBam = markDuplicates(
            allBams, 
            params.dedup
        ).duplicateMarkedBam

        // Combine with intervals for HaplotypeCaller
        bamHaplotypeCaller = duplicateMarkedBam
            .combine(interval_ch)
            .map { sampleId, bam, bai, interval -> 
                tuple(sampleId, bam, bai, interval)
            }

        // Run HaplotypeCaller
        gvcfHaplotypeCaller = haplotypeCaller(
            bamHaplotypeCaller,
            tuple(file(params.ref), file(params.fai), file(params.dict))
        ).gvcfHaplotypeCaller

        // Group gVCFs by sample
        gvcfHaplotypeCaller = gvcfHaplotypeCaller
            .groupTuple()
            .map { sampleId, gvcfs -> 
                // Ensure tuple has at least two elements with proper structure
                tuple(sampleId, gvcfs) 
            }
            // Verify that each tuple has at least two elements
            .map { it -> 
                if (it != null && it.size() < 2) {
                    log.warn "Warning: Tuple for sample ${it ? it[0] : 'unknown'} has fewer than 2 elements in gvcfHaplotypeCaller"
                }
                return it
            }

        // Gather gVCFs
        gatheredgVCFs = gatherVCFs(gvcfHaplotypeCaller)
            // Verify that each tuple has at least two elements
            .map { it -> 
                if (it != null && it.size() < 2) {
                    log.warn "Warning: Tuple for sample ${it ? it[0] : 'unknown'} has fewer than 2 elements in gatheredgVCFs"
                }
                return it
            }

        // Index gVCFs
        indexedgVCFs = indexgVCF(gatheredgVCFs)
            // Verify that each tuple has at least two elements
            .map { it -> 
                if (it != null && it.size() < 2) {
                    log.warn "Warning: Tuple for sample ${it ? it[0] : 'unknown'} has fewer than 2 elements"
                }
                return it

            }

        // Wait for indexgVCF to complete, then create cohort map
        indexedgVCFs
            .collect()
            .map { it -> params.gVCF_folder }
            .set { gvcf_folder_ch_gather }

        createCohortMap(gvcf_folder_ch_gather)
            .set { cohortMapFile }

        // Combine cohort map with regions and run final steps
        cohortMap_intervals = cohortMapFile
            .combine(regions_ch)
            .map { mapFile, region -> tuple(mapFile, region) }

        // Run genomicsDB
        genomicsDB(cohortMap_intervals)

        // Run final genotyping
        genotypeGVCFs(
            genomicsDB.out,
            params.ref,
            params.fai,
            params.dict
        )
        
        // Create a list of all VCF files for gathering
        genotypeGVCFs.out.vcf
            .map { vcf, idx -> vcf }
            .collect()
            .map { vcfs ->
                def vcfListFile = file("${params.vcfFolder}/cohort_vcf.list")
                vcfListFile.text = vcfs.join('\n')
                return vcfListFile
            }
            .set { vcf_list_ch }
            
        // Gather all VCFs into a single final VCF and index it
        gatherFinalVCFs(
            vcf_list_ch,
            file(params.ref),
            file(params.fai),
            file(params.dict)
        )
    } else if (params.entry_point == 'gathervcfs') {
        // Index gVCFs
        indexedgVCFs = indexgVCF(gatheredgVCFs)
            // Verify that each tuple has at least two elements
            .map { it -> 
                if (it != null && it.size() < 2) {
                    log.warn "Warning: Tuple for sample ${it ? it[0] : 'unknown'} has fewer than 2 elements"
                }
                return it

            }

        // Wait for indexgVCF to complete, then create cohort map
        indexedgVCFs
            .collect()
            .map { it -> params.gVCF_folder }
            .set { gvcf_folder_ch_gather }

        createCohortMap(gvcf_folder_ch_gather)
            .set { cohortMapFile }

        // Combine cohort map with regions and run final steps
        cohortMap_intervals = cohortMapFile
            .combine(regions_ch)
            .map { mapFile, region -> tuple(mapFile, region) }

        // Run genomicsDB
        genomicsDB(cohortMap_intervals)

        // Run final genotyping
        genotypeGVCFs(
            genomicsDB.out,
            params.ref,
            params.fai,
            params.dict
        )
        
        // Create a list of all VCF files for gathering
        genotypeGVCFs.out.vcf
            .map { vcf, idx -> vcf }
            .collect()
            .map { vcfs ->
                def vcfListFile = file("${params.vcfFolder}/cohort_vcf.list")
                vcfListFile.text = vcfs.join('\n')
                return vcfListFile
            }
            .set { vcf_list_ch }
            
        // Gather all VCFs into a single final VCF and index it
        gatherFinalVCFs(
            vcf_list_ch,
            file(params.ref),
            file(params.fai),
            file(params.dict)
        )
        
        // Select SNPs and non-SNPs from the final VCF
        selectSNP(
            gatherFinalVCFs.out.final_vcf,
            file(params.ref),
            file(params.fai),
            file(params.dict)
        )
        
        selectNonSNP(
            gatherFinalVCFs.out.final_vcf,
            file(params.ref),
            file(params.fai),
            file(params.dict)
        )
        
        // Filter SNPs and non-SNPs
        filterSNPs(
            selectSNP.out.snp_vcf,
            file(params.ref),
            file(params.fai),
            file(params.dict)
        )
        
        filterNonSNPs(
            selectNonSNP.out.nonsnp_vcf,
            file(params.ref),
            file(params.fai),
            file(params.dict)
        )
        
        // Annotate filtered SNPs and non-SNPs with snpEff
        annotateSNPs(
            filterSNPs.out.filtered_snp_vcf,
            params.snpeff_path,
            params.snpeff_config,
            params.genome_version
        )
        
        annotateNonSNPs(
            filterNonSNPs.out.filtered_nonsnp_vcf,
            params.snpeff_path,
            params.snpeff_config,
            params.genome_version
        )
        
        // Merge annotated SNPs and non-SNPs into final VCF
        mergeAnnotatedVCFs(
            annotateSNPs.out.annotated_snp_vcf,
            annotateNonSNPs.out.annotated_nonsnp_vcf,
            file(params.ref),
            file(params.fai),
            file(params.dict)
        )
    }
}
