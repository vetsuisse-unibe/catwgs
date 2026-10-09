#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { runFastp }           from '../modules/local/fastp/main'
include { runBWA }             from '../modules/local/bwa-mem/main'
include { mergeBams }          from '../modules/local/mergeBams/main'
include { markDuplicates }     from '../modules/local/markDuplicates/main'
include { haplotypeCaller }    from '../modules/local/hc/main'
include { gatherVCFs }         from '../modules/local/gatherVCFs/main'
include { indexgVCF }          from '../modules/local/indexgVCF/main'
include { createCohortMap }    from '../modules/local/createCohortMap/main'
include { genomicsDB }         from '../modules/local/genomicsDB/main'
include { genotypeGVCFs }      from '../modules/local/genotypeGVCFs/main'
include { gatherFinalVCFs }    from '../modules/local/gatherFinalVCFs/main'
include { selectSNP }          from '../modules/local/selectSNP/main'
include { selectNonSNP }       from '../modules/local/selectNonSNP/main'
include { filterSNPs }         from '../modules/local/filterSNPs/main'
include { filterNonSNPs }      from '../modules/local/filterNonSNPs/main'
include { annotateSNPs }       from '../modules/local/annotateSNPs/main'
include { annotateNonSNPs }    from '../modules/local/annotateNonSNPs/main'
include { mergeAnnotatedVCFs } from '../modules/local/mergeAnnotatedVCFs/main'

// Entry points:
//   markdup           : FASTQ -> fastp -> BWA -> merge lanes -> MarkDuplicates (stops at dedup BAMs)
//   start             : as markdup, then HaplotypeCaller -> per-sample gVCFs (stops at gVCFs)
//   haplotypecaller   : existing dedup BAMs in params.dedup -> per-sample gVCFs
//   cohortmap         : all gVCFs in params.gVCF_folder -> joint genotyping -> filtered, annotated cohort VCF
//   variantprocessing : existing regional VCFs in params.vcfFolder -> filtered, annotated cohort VCF
params.entry_point = 'start'

// FASTQ sample sheet -> deduplicated BAMs
workflow ALIGN_AND_MARKDUP {
    take:
        samples_file

    main:
        base_channel = Channel
            .fromPath(samples_file)
            .ifEmpty { error "Cannot find samples file: ${samples_file}" }
            .splitCsv(header: true, sep: ';', strip: true)
            .map { row ->
                [fastQbasename: row.fastQbasename, sampleID: row.sampleID, libraryID: row.libraryID, rgID: row.rgID,
                 platform: row.platform, model: row.model, center: row.center, date: row.run_date, PU: row.pu,
                 R1: file(row.R1, checkIfExists: true), R2: file(row.R2, checkIfExists: true)]
            }

        fastp_input = base_channel.map { row -> [row.fastQbasename, row.R1, row.R2] }

        bwa_meta = base_channel.map { row ->
            [row.fastQbasename, row.sampleID, row.libraryID, row.rgID,
             row.platform, row.model, row.date, row.center, row.PU]
        }

        trimmedReads = runFastp(fastp_input).trimmedReads

        mappedBams = runBWA(
            trimmedReads.join(bwa_meta, by: 0),
            Channel.value(tuple(params.assembly, params.ref))
        ).mappedBams

        lanes = mappedBams
            .groupTuple(by: 1)
            .branch {
                singleLane: it[0].size() == 1
                multipleLanes: it[0].size() > 1
            }

        singleLaneBams = lanes.singleLane.map { fqs, sampleID, bams, bais -> tuple(sampleID, bams[0], bais[0]) }
        mergedBams = mergeBams(lanes.multipleLanes).inputMarkDuplicates

        dedup = markDuplicates(singleLaneBams.mix(mergedBams), params.dedup)

    emit:
        bams = dedup.duplicateMarkedBam
}

// Dedup BAMs -> per-sample gathered and indexed gVCFs
workflow CALL_GVCFS {
    take:
        bams

    main:
        intervals = Channel
            .fromPath("${params.intervals_folder}/*-scattered.interval_list")
            .ifEmpty { error "No interval files found in ${params.intervals_folder}" }

        gvcfs = haplotypeCaller(
            bams.combine(intervals),
            Channel.value(tuple(file(params.ref), file(params.fai), file(params.dict)))
        ).gvcfHaplotypeCaller

        gathered = gatherVCFs(gvcfs.groupTuple(by: 0)).gatheredgVCFs
        indexed = indexgVCF(gathered).indexedgVCFs

    emit:
        gvcfs = gathered
        indexed = indexed
}

// Regional VCFs -> single cohort VCF -> select, filter, annotate, merge
workflow VARIANT_PROCESSING {
    take:
        regional_vcfs   // value channel: list of VCF paths

    main:
        ref_ch  = Channel.value(file(params.ref))
        fai_ch  = Channel.value(file(params.fai))
        dict_ch = Channel.value(file(params.dict))

        final_vcf = gatherFinalVCFs(regional_vcfs, ref_ch, fai_ch, dict_ch).final_vcf

        snp_raw    = selectSNP(final_vcf, ref_ch, fai_ch, dict_ch).snp_vcf
        nonsnp_raw = selectNonSNP(final_vcf, ref_ch, fai_ch, dict_ch).nonsnp_vcf

        snp_filtered    = filterSNPs(snp_raw, ref_ch, fai_ch, dict_ch).filtered_snp_vcf
        nonsnp_filtered = filterNonSNPs(nonsnp_raw, ref_ch, fai_ch, dict_ch).filtered_nonsnp_vcf

        snp_annotated = annotateSNPs(
            snp_filtered, Channel.value(params.snpeff_path), Channel.value(params.snpeff_config), Channel.value(params.genome_version)
        ).annotated_snp_vcf
        nonsnp_annotated = annotateNonSNPs(
            nonsnp_filtered, Channel.value(params.snpeff_path), Channel.value(params.snpeff_config), Channel.value(params.genome_version)
        ).annotated_nonsnp_vcf

        merged = mergeAnnotatedVCFs(snp_annotated, nonsnp_annotated, ref_ch, fai_ch, dict_ch)

    emit:
        annotated_vcf = merged.final_annotated_vcf
}

workflow CATWGS {
    def valid_entry_points = ['markdup', 'start', 'haplotypecaller', 'cohortmap', 'variantprocessing']
    if (!(params.entry_point in valid_entry_points)) {
        error "Invalid entry_point '${params.entry_point}'. Valid options: ${valid_entry_points.join(', ')}"
    }
    log.info "CATWGS entry point: ${params.entry_point}"

    if (params.entry_point in ['markdup', 'start']) {
        if (!params.samples) error "--samples is required for entry_point '${params.entry_point}'"
        dedup_bams = ALIGN_AND_MARKDUP(params.samples).bams

        if (params.entry_point == 'markdup') {
            dedup_bams.view { sampleID, bam, bai -> "Completed: ${sampleID} -> ${bam}" }
        } else {
            CALL_GVCFS(dedup_bams)
        }

    } else if (params.entry_point == 'haplotypecaller') {
        // Unique sample IDs (column 2) from the sample sheet; BAMs must already exist in params.dedup
        def sampleIDs = file(params.samples, checkIfExists: true)
            .readLines()
            .drop(1)
            .findAll { it.trim() }
            .collect { it.split(';')[1].trim() }
            .unique()
        log.info "Found ${sampleIDs.size()} unique samples in ${params.samples}"

        def bam_tuples = sampleIDs.collect { sampleID ->
            tuple(sampleID,
                  file("${params.dedup}/${sampleID}.dedup.bam", checkIfExists: true),
                  file("${params.dedup}/${sampleID}.dedup.bai", checkIfExists: true))
        }
        CALL_GVCFS(Channel.fromList(bam_tuples))

    } else if (params.entry_point == 'cohortmap') {
        regions = Channel
            .fromPath(params.regionsFile, checkIfExists: true)
            .splitCsv(header: false)

        cohortMapFile = createCohortMap(Channel.value(params.gVCF_folder)).cohortMapFile
        db = genomicsDB(cohortMapFile.combine(regions)).genomicsdb_out
        genotyped = genotypeGVCFs(db, Channel.value(file(params.ref)), Channel.value(file(params.fai)), Channel.value(file(params.dict))).vcf

        VARIANT_PROCESSING(genotyped.map { vcf, tbi -> vcf }.collect())

    } else if (params.entry_point == 'variantprocessing') {
        // Mixing plain and bgzipped chunks would gather two (possibly different) cohorts
        def vcfDirFiles = file(params.vcfFolder, checkIfExists: true).listFiles()*.name
        def n_plain = vcfDirFiles.count { it.endsWith('.vcf') }
        def n_gz    = vcfDirFiles.count { it.endsWith('.vcf.gz') }
        if (!params.regional_vcf_glob && n_plain && n_gz) {
            error "${params.vcfFolder} contains ${n_plain} .vcf and ${n_gz} .vcf.gz chunks. " +
                  "Pick one set with --regional_vcf_glob '*.vcf' or '*.vcf.gz'."
        }
        def vcf_glob = params.regional_vcf_glob ?: (n_gz ? '*.vcf.gz' : '*.vcf')

        log.info "Starting variant processing from ${params.vcfFolder}/${vcf_glob}"
        regional_vcfs = Channel
            .fromPath("${params.vcfFolder}/${vcf_glob}")
            .ifEmpty { error "No files matching ${vcf_glob} in ${params.vcfFolder}" }
            .collect()

        VARIANT_PROCESSING(regional_vcfs)
    }
}

// Lets this file be run directly; ignored when CATWGS is included from ../main.nf
workflow {
    CATWGS()
}
