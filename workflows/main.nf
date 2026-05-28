#!/usr/bin/env nextflow

// Assuming the nf-core/modules repository is accessible
include { runFastp } from '/data/projects/p926_Lynx_whole_genome_analysis/dsl2/modules/local/fastp/main'
include { runBWA } from '/data/projects/p926_Lynx_whole_genome_analysis/dsl2/modules/local/bwa-mem/main'
include { mergeBams } from '/data/projects/p926_Lynx_whole_genome_analysis/dsl2/modules/local/mergeBams/main'
include { markDuplicates } from '/data/projects/p926_Lynx_whole_genome_analysis/dsl2/modules/local/markDuplicates/main'
include { haplotypeCaller } from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/hc/main'
include { gatherVCFs } from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/gatherVCFs/main'  
include { indexgVCF } from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/indexgVCF/main'
include {createCohortMapFile} from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/createCohortMapFile/main'
// include {createCohortMap} from '/data/projects/p926_Lynx_whole_genome_analysis/dsl2/modules/local/createCohortMap/main' // Commented out - using local definition
include {genomicsDB} from '/data/projects/p926_Lynx_whole_genome_analysis/dsl2/modules/local/genomicsDB/main'
include {genotypeGVCFs} from '/data/projects/p926_Lynx_whole_genome_analysis/dsl2/modules/local/genotypeGVCFs/main'
include {gatherFinalVCFs} from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/gatherFinalVCFs/main'
include {selectSNP} from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/selectSNP/main'
include {selectNonSNP} from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/selectNonSNP/main'
include {filterSNPs} from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/filterSNPs/main'
include {filterNonSNPs} from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/filterNonSNPs/main'
include {annotateSNPs} from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/annotateSNPs/main'
include {annotateNonSNPs} from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/annotateNonSNPs/main'
include {mergeAnnotatedVCFs} from '/data/projects/p531_Felis_Catus__whole_genome_Analysis/nextFlow/dsl2/modules/local/mergeAnnotatedVCFs/main'

params.entry_point = 'start' // Options: start, haplotypecaller, cohortmap, variantprocessing

// Define the createCohortMap process once
process createCohortMap {
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

workflow LYNXWGS { 
    if (params.entry_point == 'start') {
        Channel
            .fromPath(params.samples)
            .splitCsv(header: true, sep: ';', strip: true)
            .map { row ->
                [fastQbasename: row.fastQbasename, sampleID: row.sampleID, libraryID: row.libraryID, rgID: row.rgID, 
                 platform: row.platform, model: row.model, center: row.center, date: row.run_date, PU: row.pu,
                 R1: file(row.R1), R2: file(row.R2)]
            }
            .set { base_channel }

        // channel for interval list
        interval_ch=Channel.fromPath(params.intervals_folder + '*-scattered.interval_list')

        // channel for fastqc 
        base_channel
        .map { row ->
            [row.fastQbasename,row.R1, row.R2]
        }
        .set { fastqc_input }

        // channel for bwa-mem
        base_channel
        .map { row ->
            [row.fastQbasename, row.sampleID, row.libraryID, row.rgID, 
             row.platform, row.model, row.date, row.center, row.PU]
        }
        .set { bwa_mem_input }

        // Read intervals for genomicsDB
        Channel
        .fromPath(params.regionsFile)
        .splitCsv(header: false)
        .set { regions_ch }

        //run fastp
        trimmedReads=runFastp(fastqc_input).trimmedReads
        
        //combine the trimmedReads (runFastp output) with bwa_mem_input channel
        trimmedReads
        .join(bwa_mem_input, by: 0)
        .set { combined_for_bwa_mem }
        
        //run bwa-mem
        mappedBams=runBWA(combined_for_bwa_mem,params.assembly,params.ref).mappedBams
        
        //split the mappedBams by sampleID and the number of mappedBams per sample
        mappedBams
        .groupTuple(by: 1)
        .set { mappedBams_by_sampleID }
        
        // Use the branch operator to split the data into two channels based on the lane count
        mappedBams_by_sampleID.branch {
            singleLane: it[0].size() == 1
            multipleLanes: it[0].size() > 1
        }
        .set { channels }
        
        // Assigning to new channels
        channels.singleLane.set { singleLaneSamples }
        channels.multipleLanes.set { multipleLaneSamples }
        
        //merge the bams from multiple lanes
        mergedBams=mergeBams(multipleLaneSamples).inputMarkDuplicates
        singleLaneSamples = singleLaneSamples.map { item ->
            return [item[1], item[2][0], item[3][0]]
        }
        
        // Merge single lane and multiple lane samples using mix and ifEmpty
        allBams = singleLaneSamples
            .mix(mergedBams)
            .ifEmpty { mergedBams }

        duplicateMarkedBam=markDuplicates(allBams).duplicateMarkedBam
        
        // Create interval channel for haplotype calling
        interval_ch = Channel.fromPath(params.intervals_folder + '/*-scattered.interval_list')
        
        // Combine deduped BAMs with intervals for scattered haplotype calling
        bam_intervals = duplicateMarkedBam.combine(interval_ch)
        
        // Run HaplotypeCaller to generate gVCFs
        gvcfs = haplotypeCaller(bam_intervals, Channel.value(tuple(file(params.ref), file(params.fai), file(params.dict)))).gvcfHaplotypeCaller
        
        // Group gVCFs by sample for gathering
        gvcfs_by_sample = gvcfs.groupTuple(by: 0)
        
        // Gather scattered gVCFs per sample
        gatheredgVCFs = gatherVCFs(gvcfs_by_sample).gatheredgVCFs
        
        // Index the gathered gVCFs
        indexedgVCFs = indexgVCF(gatheredgVCFs).indexedgVCFs
        
        // Read intervals for genomicsDB
        regions_ch = Channel
            .fromPath(params.regionsFile)
            .splitCsv(header: false)
        
    } else if (params.entry_point == 'haplotypecaller') {
        // Read sample sheet and create BAM tuples directly
        def allSamples = file(params.samples)
            .readLines()
            .drop(1)  // Skip header
            .collect { line -> 
                def fields = line.split(';')
                fields[1]  // sampleID is column 2 (index 1)
            }
            .unique()
        
        log.info "Found ${allSamples.size()} unique samples from sample sheet"
        log.info "Samples: ${allSamples}"
        
        bam_tuples = allSamples.collect { sample ->
                def bam = file("${params.dedup}/${sample}.dedup.bam")
                def bai = file("${params.dedup}/${sample}.dedup.bai")
                if (!bam.exists()) error "BAM file not found: ${bam}"
                if (!bai.exists()) error "BAI file not found: ${bai}"
                tuple(sample, bam, bai)
            }
        
        log.info "Created ${bam_tuples.size()} BAM tuples for processing"
        
        // Create channels
        duplicateMarkedBam = Channel.fromList(bam_tuples)
            .view { "BAM channel emitting: $it" }
        interval_ch = Channel.fromPath(params.intervals_folder + '/*-scattered.interval_list')
            .ifEmpty { error "No interval files found in ${params.intervals_folder}" }

        // Combine deduped BAMs with intervals for scattered haplotype calling
        bam_intervals = duplicateMarkedBam.combine(interval_ch)
            .view { "Combined channel emitting: $it" }
        
        // Run HaplotypeCaller to generate gVCFs
        gvcfs = haplotypeCaller(bam_intervals, Channel.value(tuple(file(params.ref), file(params.fai), file(params.dict)))).gvcfHaplotypeCaller
        
        // Group gVCFs by sample for gathering
        gvcfs_by_sample = gvcfs.groupTuple(by: 0)
        
        // Gather scattered gVCFs per sample
        gatheredgVCFs = gatherVCFs(gvcfs_by_sample).gatheredgVCFs
        
        // Index the gathered gVCFs
        indexedgVCFs = indexgVCF(gatheredgVCFs).indexedgVCFs

        // Read intervals for genomicsDB
        regions_ch = Channel
            .fromPath(params.regionsFile)
            .splitCsv(header: false)
    } else if (params.entry_point == 'cohortmap') {
        // Read intervals for genomicsDB
        regions_ch = Channel
            .fromPath(params.regionsFile)
            .splitCsv(header: false)

        // Run createCohortMap with the gVCF folder
        cohortMapFile = createCohortMap(Channel.from(params.gVCF_folder))
    }

    // Joint genotyping + variant processing only run for the 'cohortmap' and
    // 'variantprocessing' entry points. 'start' and 'haplotypecaller' stop after
    // per-sample gVCFs are generated and indexed, so they no longer trigger a
    // redundant joint-genotyping/annotation pass over a stale cohort.
    if (params.entry_point == 'cohortmap' || params.entry_point == 'variantprocessing') {

        if (params.entry_point == 'variantprocessing') {
            // Collect all VCF files from vcfFolder
            all_vcfs = Channel.fromPath("${params.vcfFolder}/*.vcf").collect()
            log.info "Starting variant processing from existing VCFs in ${params.vcfFolder}"
        } else {
            // cohortmap: cohortMapFile and regions_ch were created in the branch above
            cohortMap_intervals = cohortMapFile.combine(regions_ch)

            // run genomicsDB and final genotyping
            genomicsDB_out = genomicsDB(cohortMap_intervals)
            genotyped_vcfs = genotypeGVCFs(genomicsDB_out, Channel.value(params.ref), Channel.value(params.fai), Channel.value(params.dict))

            // Gather all regional VCFs into a single cohort VCF
            all_vcfs = genotyped_vcfs.vcf_out.map { it[1] }.collect()
        }

        // === VARIANT PROCESSING PIPELINE ===
        final_vcf = gatherFinalVCFs(
            all_vcfs,
            Channel.value(file(params.ref)),
            Channel.value(file(params.fai)),
            Channel.value(file(params.dict))
        )

        // Reference channels for variant processing
        ref_ch = Channel.value(file(params.ref))
        fai_ch = Channel.value(file(params.fai))
        dict_ch = Channel.value(file(params.dict))

        // Select SNPs and Non-SNPs in parallel
        snp_raw = selectSNP(final_vcf.final_vcf, ref_ch, fai_ch, dict_ch)
        nonsnp_raw = selectNonSNP(final_vcf.final_vcf, ref_ch, fai_ch, dict_ch)

        // Filter SNPs and Non-SNPs
        snp_filtered = filterSNPs(snp_raw.snp_vcf, ref_ch, fai_ch, dict_ch)
        nonsnp_filtered = filterNonSNPs(nonsnp_raw.nonsnp_vcf, ref_ch, fai_ch, dict_ch)

        // Annotate with snpEff
        snp_annotated = annotateSNPs(
            snp_filtered.filtered_snp_vcf,
            Channel.value(params.snpeff_path),
            Channel.value(params.snpeff_config),
            Channel.value(params.genome_version)
        )
        nonsnp_annotated = annotateNonSNPs(
            nonsnp_filtered.filtered_nonsnp_vcf,
            Channel.value(params.snpeff_path),
            Channel.value(params.snpeff_config),
            Channel.value(params.genome_version)
        )

        // Merge annotated VCFs
        mergeAnnotatedVCFs(
            snp_annotated.annotated_snp_vcf,
            nonsnp_annotated.annotated_nonsnp_vcf,
            ref_ch,
            fai_ch,
            dict_ch
        )
    }
}