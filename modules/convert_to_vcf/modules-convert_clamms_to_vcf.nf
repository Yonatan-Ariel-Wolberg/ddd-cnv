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
    path("samples.cnv.filtered.bed")

    // Specifies output
    output:
    file("${sample_name}_clamms.vcf"), emit: clamms_vcf

    // Specifies script
    script:
    """
    python ddd-cnv/bin/convert_clamms_to_vcf.py ${outdir}/out_CLAMMS/samples.cnv.filtered.bed ${outdir}/out_CLAMMS
    """
}