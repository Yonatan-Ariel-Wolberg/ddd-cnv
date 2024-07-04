#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PARAMETERS
params.xhmm_vcf = ''
params.clamms_vcf = 'CLAMMS.cnv.filtered.vcf'
params.canoes_vcf = ''
params.survivor = '~/SURVIVOR/Debug/SURVIVOR'

// CONVERTS CLAMMS OUTPUT TO VCF
process bedtovcf_CLAMMS {
    // Tags process as 'bedtovcf_CLAMMS'
    tag { 'bedtovcf_CLAMMS' }

    // Labels process as 'survivor'
    label 'survivor'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_SURVIVOR", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path("CLAMMS.cnv.filtered.bed") from filtered_cnvs

    // Specifies output
    output:
    path("CLAMMS.cnv.filtered.bed"), emit: clamms_vcf

    // Specifies script
    script:
    """
    ${params.survivor} bedtovcf ${outdir}/out_CLAMMS/CLAMMS.cnv.filtered.bed ${ref} CLAMMS.cnv.filtered.vcf  
    """
}

// MERGE SAMPLES USING SURVIVOR
process mergeSamples_SURVIVOR {
    // Tags process as 'merge_samples_survivor'
    tag { 'merge_samples_survivor' }

    // Labels process as 'survivor'
    label 'survivor'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_SURVIVOR", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(canoes_vcf)
    path(clamms_vcf)
    path(xhmm_vcf)

    // Specifies output
    output:
    

    // Specifies script
    script:
    """
    SURVIVOR bedtovcf ${outdir}/out_CLAMMS/CLAMMS.cnv.filtered.bed ${ref} CLAMMS.cnv.filtered.vcf
    SURVIVOR merge
    """
}
