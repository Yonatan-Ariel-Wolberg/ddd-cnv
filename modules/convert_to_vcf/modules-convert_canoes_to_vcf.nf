#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// PROCESS THAT CONVERTS CANOES OUTPUT FILE TO VCF FORMAT
process canoes_to_vcf {
    // Tags process as 'canoes_to_vcf'
    tag { 'canoes_to_vcf' }

    // Labels process as 'canoes'
    label 'canoes'

    // Specifies destination directory where output file is to be published
    publishDir "${outidr}/out_CANOES", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path("Sample_CNVs_filtered.csv")

    // Specifies output
    output:
    file("${sample_name}_canoes.vcf"), emit: canoes_vcf

    // Specifies script
    script:
    """
    python ddd-cnv/bin/convert_canoes_to_vcf.py ${outdir}/out_CANOES/Sample_CNVs_filtered.csv ${outdir}/out_CANOES
    """
}