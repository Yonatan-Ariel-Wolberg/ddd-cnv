#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PROCESS THAT CONVERTS CLAMMS OUTPUT TO VCF FORMAT
process clamms_to_vcf {
    // Tags process as 'clamms_to_vcf'
    tag { 'clamms_to_vcf' }

    // Labels process as 'clamms'
    label 'clamms'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path("CLAMMS.cnv.filtered.bed")

    // Specifies output
    output:
    file("CLAMMS_CNVs_filtered.vcf"), emit: clamms_vcf

    // Specifies script
    script:
    """
    python ddd-cnv/bin/clamms_bed_to_vcf.py -i ${outdir}/out_CLAMMS/CLAMMS.cnv.filtered.bed -o ${outdir}/out_CLAMMS/CLAMMS_CNVs_filtered.vcf
    """
}