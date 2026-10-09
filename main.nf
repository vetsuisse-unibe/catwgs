#!/usr/bin/env nextflow

/*
========================================================================================
    Felis Catus Whole Genome Analysis Pipeline
========================================================================================
    Github : https://github.com/vetsuisse-unibe/catwgs
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

log.info """\
         F E L I S   C A T U S   W G S   P I P E L I N E
         ===================================
         entry_point    : ${params.entry_point}
         samples        : ${params.samples}
         assembly       : ${params.assembly}
         ref            : ${params.ref}
         outdir         : ${params.outdir}
         """
         .stripIndent()

include { CATWGS } from './workflows/main_v2'

workflow {
    CATWGS ()
}
