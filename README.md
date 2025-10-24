# CATWGS: Felis Catus Whole Genome Sequencing Pipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17434829.svg)](https://doi.org/10.5281/zenodo.17434829)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)

A Nextflow DSL2 pipeline for processing *Felis catus* (domestic cat) whole genome sequencing data using BWA-MEM2 alignment and GATK best practices for variant calling.

## Overview

This pipeline implements a comprehensive whole genome sequencing analysis workflow including:
- **Quality Control & Trimming**: FastP for adapter trimming and quality filtering
- **Alignment**: BWA-MEM2 for efficient mapping to the reference genome
- **BAM Processing**: Merging multi-lane samples and marking duplicates
- **Variant Calling**: GATK HaplotypeCaller for gVCF generation
- **Joint Genotyping**: GenomicsDBImport and GenotypeGVCFs for cohort analysis
- **Variant Filtering**: SNP and INDEL filtering using GATK best practices
- **Annotation**: snpEff for functional variant annotation

## Pipeline Workflow

```
FASTQ Files
    ↓
FastP (Quality Control & Trimming)
    ↓
BWA-MEM2 (Alignment)
    ↓
Merge BAMs (for multi-lane samples)
    ↓
Mark Duplicates (Picard)
    ↓
HaplotypeCaller (per-sample gVCF)
    ↓
Gather VCFs (merge intervals per sample)
    ↓
Index gVCFs
    ↓
Create Cohort Map
    ↓
GenomicsDBImport (joint database)
    ↓
GenotypeGVCFs (cohort calling)
    ↓
Gather Final VCFs (merge regions)
    ↓
Select SNPs / Select INDELs
    ↓
Filter Variants
    ↓
Annotate with snpEff
    ↓
Final Annotated VCF
    ↓
[Optional: Variant Extraction Pipeline]
    ↓
Tabular Output (TSV format)
```

**Note**: The variant extraction pipeline is a post-processing step that converts the final VCF into tabular format for downstream analysis.

## Reference Genome

- **Species**: *Felis catus*
- **Assembly**: F.catus_Fca126_mat1.0 (NCBI)
- **Location**: `/data/references/Felis_catus/NCBI/F.catus_Fca126_mat1.0/`

## Requirements

### Software Dependencies
- **Nextflow**: >=21.10.3
- **Java**: 17.0.6
- **FastP**: 0.23.4
- **BWA-MEM2**: 2.2.1
- **SAMtools**: 1.13
- **Picard**: 2.25.5
- **GATK**: 4.2.6.1
- **snpEff**: (custom installation)

### System Requirements
- **Executor**: SLURM
- **Partition**: pibu_el8
- **Max Memory**: 128 GB
- **Max CPUs**: 16
- **Max Time**: 72 hours

## Input Data

### Sample Sheet Format

The pipeline requires a semicolon-separated sample sheet with the following headers:

```
fastQbasename;sampleID;libraryID;rgID;platform;model;center;run_date;pu;R1;R2
```

**Example**:
```
K0806_L001;K0806;K0806_LIB;HJYMJDSXF.1;ILLUMINA;Novogene;Unibe;2025-09-25;HJYMJDSXF.1.CCGCGGTT+CTAGCGCT;/path/to/R1.fastq.gz;/path/to/R2.fastq.gz
```

**Field Descriptions**:
- `fastQbasename`: Unique identifier for FASTQ pair (sample_lane)
- `sampleID`: Sample identifier (for grouping multi-lane data)
- `libraryID`: Library preparation identifier
- `rgID`: Read group ID (flowcell.lane)
- `platform`: Sequencing platform (ILLUMINA)
- `model`: Sequencing provider/model
- `center`: Sequencing center
- `run_date`: Sequencing run date
- `pu`: Platform unit (flowcell.lane.barcode)
- `R1`: Full path to forward reads
- `R2`: Full path to reverse reads

### Multi-Lane Samples

The pipeline automatically handles:
- Multiple lanes per sample (L1, L2, etc.)
- Multiple flowcells per lane
- Proper BAM merging before duplicate marking

## Usage

### Running the Full Pipeline

```bash
sbatch run_catwgs_pipeline_v2.sh
```

This script submits the pipeline as a SLURM job with the following specifications:
- **Job name**: catwgs_v2_pipeline
- **Partition**: pibu_el8
- **CPUs**: 4
- **Memory**: 16 GB
- **Time limit**: 172 hours (7 days)

### Pipeline Entry Points

The pipeline supports three entry points for resuming from different stages:

#### 1. Full Pipeline (Default)
Start from FastP and run the complete workflow:

```bash
nextflow run workflows/main_v2.nf \
    -profile unibe \
    --entry_point fastp \
    -resume
```

#### 2. From Mark Duplicates
Resume from mark duplicates step (requires pre-aligned BAMs):

```bash
sbatch run_catwgs_pipeline_v2_markdup.sh
```

Or manually:
```bash
nextflow run workflows/main_v2.nf \
    -profile unibe \
    --entry_point markduplicates \
    --input_bams bam_list.txt \
    -resume
```

#### 3. From Gather VCFs
Resume from gVCF gathering step (requires per-sample gVCFs):

```bash
sbatch run_catwgs_pipeline_v2_gathervcfs.sh
```

Or manually:
```bash
nextflow run workflows/main_v2.nf \
    -profile unibe \
    --entry_point gathervcfs \
    --input_gvcfs "path/to/gvcfs/*.g.vcf.gz" \
    -resume
```

## Configuration

The main configuration file is `nextflow.config`. Key parameters:

### Input/Output Paths
```groovy
params {
    samples         = "assets/sampleSheet_Sep2025.txt"
    ref             = "/data/references/.../F.catus_Fca126_mat1.0_genomic.fa"
    intervals_folder = "/data/references/.../interval-files-folder-48/"
    regionsFile     = "regions.wholeGenome.txt"
    dedup           = "/data/projects/.../dedup_bams/Fca126"
    gVCF_folder     = "/data/projects/.../gVCFs/"
    vcfFolder       = "vcf/"
    finalVcfFolder  = "final_vcf/"
}
```

### Resource Configuration

Process-specific resources are defined per tool:

| Process | CPUs | Memory | Time | Retries |
|---------|------|--------|------|---------|
| FastP | 4 | 8 GB | 4 h | 3 |
| BWA-MEM2 | 8 | 32 GB | 8 h | 3 |
| Mark Duplicates | 1 | 30 GB | 12 h | 3 |
| HaplotypeCaller | 1 | 30 GB | 12 h | 3 |
| GenomicsDB | 2 | 50 GB | 24 h | 3 |
| GenotypeGVCFs | 1 | 30 GB | 12 h | 3 |

Resources automatically scale with retry attempts.

## Output Structure

```
results/
├── fastp/                    # Quality control reports
├── bwa/                      # Aligned BAM files
├── merged_bams/              # Merged multi-lane BAMs
├── dedup_bams/               # Duplicate-marked BAMs
├── gVCFs/                    # Per-sample gVCF files
├── cohortMap/                # Cohort map for GenomicsDB
├── genomicsDB/               # GenomicsDB workspaces
├── vcf/                      # Per-region VCF files
├── final_vcf/                # Final merged VCF
│   ├── cohort.vcf.gz
│   ├── cohort_SNPs.vcf.gz
│   ├── cohort_SNPs_filtered.vcf.gz
│   ├── cohort_SNPs_filtered_annotated.vcf.gz
│   └── cohort_final_annotated.vcf.gz
└── reports/                  # Nextflow execution reports
    ├── nextflow_v2_report.html
    ├── nextflow_v2_timeline.html
    └── nextflow_v2_trace.txt
```

## Execution Reports

The pipeline generates comprehensive execution reports:

- **HTML Report**: Resource usage, task statistics, execution time
- **Timeline**: Visual timeline of task execution
- **Trace**: Detailed trace file with all task metrics

Reports are saved to `reports/` with timestamps.

## Resume Capability

The pipeline supports Nextflow's `-resume` feature:

```bash
nextflow run workflows/main_v2.nf -profile unibe -resume
```

This will:
- Skip completed tasks
- Use cached results where possible
- Resume from failed/interrupted tasks
- Enable lenient caching for HPC environments

## Monitoring

Check pipeline progress:

```bash
# View SLURM job status
squeue -u $USER

# Check Nextflow logs
tail -f logs/nextflow_v2_*.out

# Monitor resource usage
cat .nextflow.log
```

## Troubleshooting

### Common Issues

1. **Out of Memory Errors**
   - Resources automatically increase with retries
   - Check error logs: `logs/nextflow_v2_*.err`

2. **File Not Found Errors**
   - Verify sample sheet paths are absolute
   - Ensure reference files exist
   - Check FASTQ file accessibility

3. **SLURM Job Failures**
   - Review job output: `work/<hash>/.command.out`
   - Check job error: `work/<hash>/.command.err`
   - Verify partition availability: `sinfo`

4. **Pipeline Crashes**
   - Use `-resume` to continue from last successful step
   - Check `.nextflow.log` for detailed errors
   - Ensure sufficient disk space in work directory

### Getting Help

Check execution reports and logs:
```bash
# View last run report
ls -lt reports/ | head -n 4

# Check Nextflow cache
nextflow log last -f hash,name,status,exit,container
```

## Project Structure

```
dsl2/
├── workflows/
│   └── main_v2.nf              # Main workflow file
├── modules/
│   └── local/                  # Process modules
│       ├── fastp/
│       ├── bwa-mem/
│       ├── markDuplicates/
│       ├── haplotypeCaller/
│       └── ...
├── assets/
│   └── sampleSheet_Sep2025.txt # Sample sheet
├── conf/                       # Additional configs
├── bin/                        # Helper scripts
├── nextflow.config             # Main configuration
├── run_catwgs_pipeline_v2.sh   # SLURM submission script
└── README.md                   # This file
```

## Performance Optimization

The pipeline includes several optimizations:

- **Parallel Execution**: Multiple samples/intervals processed simultaneously
- **Resource Scaling**: Automatic retry with increased resources
- **Scratch Space**: Uses SLURM temporary storage
- **Symlinks**: Reduces disk I/O and storage requirements
- **Lenient Cache**: Better resume capability on HPC
- **Queue Management**: Rate limiting prevents scheduler overload

## Variant Calling Strategy

### HaplotypeCaller
- Runs per-sample on 48 genomic intervals
- Generates gVCF files for joint genotyping
- Uses `-ERC GVCF` mode for all sites

### GenomicsDB Import
- Creates genomic databases per region
- Enables efficient cohort-wide joint calling
- Handles large cohorts (100+ samples)

### Filtering Thresholds

**SNPs**:
- QD < 2.0
- FS > 60.0
- MQ < 40.0
- MQRankSum < -12.5
- ReadPosRankSum < -8.0

**INDELs**:
- QD < 2.0
- FS > 200.0
- ReadPosRankSum < -20.0

## Variant Extraction Pipeline

After the main Nextflow pipeline completes, you can extract variant information into a tabular format for downstream analysis using the extraction pipeline.

### Overview

The extraction pipeline processes the final annotated VCF file and extracts:
- Variant coordinates (CHROM, POS, REF, ALT)
- Genotypes for all samples
- Detailed snpEff annotations (gene names, effects, impacts, HGVS notation)

The pipeline uses parallel SLURM array jobs to process 2478 genomic regions simultaneously, then merges results into a single tabular file.

### Usage

Run the complete extraction pipeline:

```bash
bash submit_extract_pipeline.sh
```

This master script automatically:
1. Creates a header file with sample names and annotation fields
2. Submits a SLURM array job to extract 2478 regions in parallel
3. Submits a merge job (with dependency) to combine all results

### Pipeline Steps

#### Step 1: Create Header File

```bash
bash create_header.sh
```

**Output**: `head.txt` containing column headers:
- Basic fields: CHROM, POS, REF, ALT, FILTER
- Sample IDs (genotype columns for all 124 samples)
- Annotation fields: ALLELE, EFFECT, IMPACT, GENE, GENEID, FEATURE, etc.

#### Step 2: Extract Regions (Parallel)

```bash
sbatch extract_regions.sh
```

**SLURM Configuration**:
- **Array size**: 2478 tasks (one per genomic region)
- **CPUs**: 16 per task
- **Memory**: 20 GB per task
- **Time limit**: 36 hours
- **Partition**: pibu_el8

**Process per region**:
1. Uses `tabix` to extract region from VCF
2. Runs `vcfEffOnePerLine.pl` to split multi-annotation lines
3. Uses SnpSift `extractFields` to extract specific fields
4. Compresses output with `zstd`

**Output**: `extract_results/<region>.zst` (compressed tabular data per region)

#### Step 3: Merge Results

```bash
sbatch merge_extract_regions.sh
```

**SLURM Configuration**:
- **Dependency**: Waits for all extraction tasks to complete
- **CPUs**: 2
- **Memory**: 16 GB
- **Time limit**: 48 hours

**Process**:
1. Combines header with all region results in genomic order
2. Decompresses zstd files and filters headers
3. Creates final output: `final_output.txt`

### Input Files

The extraction pipeline requires:
- `final_vcf/cohort_124.var.flt.ann.vcf.gz` - Final annotated VCF from main pipeline
- `final_vcf/cohort_124.var.flt.ann.vcf.gz.tbi` - Tabix index
- `regions.wholeGenome.txt` - List of 2478 genomic regions (CHROM:START-END)

### Output Files

```
extract_results/
├── SUPER_1_1_1000000.zst       # Compressed per-region extracts (2478 files)
├── SUPER_1_1000001_2000000.zst
├── ...
└── chrZ_73680001_73817304.zst

final_output.txt                 # Final merged tabular output
head.txt                         # Header file
```

### Extracted Fields

The final output contains these columns:

**Variant Information**:
- `CHROM` - Chromosome
- `POS` - Position
- `REF` - Reference allele
- `ALT` - Alternate allele
- `FILTER` - Filter status (PASS/filtered)

**Genotypes** (124 sample columns):
- `K0001`, `K0002`, ..., `K0124` - Genotype for each sample (e.g., 0/0, 0/1, 1/1)

**snpEff Annotations**:
- `ALLELE` - Variant allele
- `EFFECT` - Variant effect (e.g., missense_variant, synonymous_variant)
- `IMPACT` - Impact level (HIGH, MODERATE, LOW, MODIFIER)
- `GENE` - Gene symbol
- `GENEID` - Gene identifier
- `FEATURE` - Feature type (e.g., transcript)
- `FEATUREID` - Feature identifier (e.g., transcript ID)
- `BIOTYPE` - Biotype (e.g., protein_coding)
- `RANK` - Exon/intron rank
- `HGVS_C` - HGVS coding sequence notation
- `HGVS_P` - HGVS protein notation
- `CDNA_POS/LEN` - Position and length in cDNA
- `CDS_POS/LEN` - Position and length in coding sequence
- `AA_POS/LEN` - Position and length in amino acid sequence
- `DISTANCE` - Distance to nearest feature
- `ERRORS` - Annotation errors/warnings

### Monitoring the Extraction Pipeline

Check job status:

```bash
# View all jobs
squeue -u $USER

# Check specific jobs (using Job IDs from submission output)
squeue -j <JOB_ID>

# Monitor extraction progress
ls extract_results/*.zst | wc -l
# Expected: 2478 files

# Check for failed array tasks
grep -i error logs/region_*.err

# View merge job progress
tail -f logs/merge_*.out
```

### Troubleshooting Extraction Pipeline

**Issue: Missing regions in output**
```bash
# Check which regions failed
while read region; do
    filename="extract_results/$(echo $region | tr ':' '_' | tr '-' '_').zst"
    [ ! -f "$filename" ] && echo "Missing: $region"
done < regions.wholeGenome.txt
```

**Issue: Array job failures**
```bash
# Review failed task logs
grep -i "error\|failed" logs/region_*.err

# Resubmit specific array tasks
sbatch --array=123,456,789 extract_regions.sh
```

**Issue: Insufficient disk space**
```bash
# Check disk usage
du -sh extract_results/
df -h $PWD

# Clean up if needed (after successful merge)
rm extract_results/*.zst
```

### Performance Notes

- **Parallel Efficiency**: 2478 regions processed simultaneously
- **Compression**: zstd compression reduces storage by ~80%
- **Total Runtime**: ~8-12 hours for extraction + 2-4 hours for merge
- **Output Size**: ~50-100 GB uncompressed (depending on variant count)

### Downstream Analysis

The final `final_output.txt` can be used for:
- Importing into R/Python for statistical analysis
- Loading into databases (PostgreSQL, MySQL)
- Filtering variants by effect/impact
- Identifying variants in specific genes
- Comparing genotypes across samples
- Population genetics analyses

Example loading in R:
```R
# Load data
variants <- read.table("final_output.txt", header=TRUE, sep="\t", 
                       quote="", comment.char="")

# Filter for high impact variants
high_impact <- variants[variants$IMPACT == "HIGH", ]

# Get missense variants in specific gene
gene_vars <- variants[variants$GENE == "GENE_NAME" & 
                      variants$EFFECT == "missense_variant", ]
```

## Citation

If you use this pipeline, please cite:

- **BWA-MEM2**: Vasimuddin et al., 2019
- **GATK**: McKenna et al., 2010
- **Nextflow**: Di Tommaso et al., 2017
- **snpEff**: Cingolani et al., 2012

## Authors

- **Vidhya Jagannathan**
- Institute of Genetics, University of Bern

## License

This pipeline is part of project p531_Felis_Catus__whole_genome_Analysis.

## Version History

- **v2.0** (2025): DSL2 implementation with modular architecture
- **v1.0** (Previous): Initial DSL1 implementation (legacy)

## Contact

For questions or issues:
- Email: vjaganna@unibe.ch
- Check SLURM job outputs in `logs/` directory
- Review Nextflow logs in `.nextflow.log`

---

**Last Updated**: October 2025  
**Pipeline Name**: CATWGS v2  
**Reference Assembly**: F.catus_Fca126_mat1.0
