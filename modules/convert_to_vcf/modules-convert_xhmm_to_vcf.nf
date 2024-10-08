#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters
params.outdir = 'out_xhmm'
outdir = file(params.outdir, type: 'dir')

// Convert XHMM XCNV format to VCF format
process XhmmToVcf {

    // Tags the process as 'xhmm_to_vcf'
    tag { 'xhmm_to_vcf' }

    // Labels the process as 'xhmm'
    label 'xhmm'

    // Specifies destination directory where output files are published
    publishDir "${params.outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path filtered_cnvs
    
    // Specifies output
    output:
    path "XHMM_CNVs_filtered.vcf", emit: filtered_vcfs
    
    // Specifies script
    script:
    """
    # Converts XHMM XCNV format to VCF format
    python /home/ywolberg/DDD_Africa/ddd-cnv/ddd-cnv/bin/xhmm_xcnv_to_vcf2.py -i /dataG/ddd/analysis/CNV_calling/TRUTHSET.2023.08.22/XHMM_filtered.xcnv -d ${params.outdir}/nextflow -c XHMM_CNVs_filtered.vcf 
    """
}

// Index VCF
process IndexVcf {

    // Tags the process as 'index_vcf'
    tag { 'index_vcf' }

    // Labels the process as 'xhmm'    
    label 'xhmm'

    // Specifies destination directory where output files are published
    publishDir "${params.outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path("XHMM_CNVs_filtered.vcf") from filtered_vcfs

    // Specifies output
    output:
    path("XHMM_CNVs_filtered.vcf.gz"), emit: bgzip_vcfs
    path("sorted_XHMM_CNVs_filtered.vcf.gz"), emit: sorted_bgzipped_vcfs
    path("sorted_XHMM_CNVs_filtered.vcf.gz.tbi"), emit: vcf_index
    path("sorted_XHMM_CNVs_filtered.vcf"), emit: unzipped_sorted_vcf
  
    // Specifies script
    script:
    """
    # Compress the VCF file with bgzip
    bgzip -c ${params.outdir}/nextflow/XHMM_CNVs_filtered.vcf > out_xhmm/XHMM_CNVs_filtered.vcf.gz
    
    # Load the samtools/1.20
    module load samtools/1.20

    # Sort the VCF file with BCFtools
    bcftools sort -o ${params.outdir}/nextflow/sorted_XHMM_CNVs_filtered.vcf.gz -O z ${params.outdir}/nextflow/XHMM_CNVs_filtered.vcf.gz 

    # Index the VCF file with Tabix
    tabix -p vcf ${params.outdir}/nextflow/sorted_XHMM_CNVs_filtered.vcf.gz

    # Unzip the sorted VCF file and redirects output to VCF file
    zcat ${params.outdir}/nextflow/sorted_XHMM_CNVs_filtered.vcf.gz | less -S > ${params.outdir}/nextflow/sorted_XHMM_CNVs_filtered.vcf
    """
}

workflow {
// Input channel
filtered_cnvs = Channel.fromPath('/dataG/ddd/analysis/CNV_calling/TRUTHSET.2023.08.22/XHMM_filtered.xcnv')

XhmmToVcf(filtered_cnvs)
IndexVcf(filtered_vcfs)
}