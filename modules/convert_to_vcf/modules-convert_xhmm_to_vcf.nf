#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Parameters

// Convert XHMM XCNV format to VCF format
process xhmm_to_vcf {
    // Tags the process as 'xhmm_to_vcf'
    tag { 'xhmm_to_vcf' }

    // Labels the process as 'xhmm'
    label 'xhmm'

    // Specifies destination directory where output files are published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite=true

    // Specifies input
    input:
    path(filtered_cnvs)
    
    // Specifies output
    output:
    path(*), emit: filtered_vcfs
    
    // Specifies script
    script:
    """
    # Converts XHMM XCNV format to VCF format
    python /home/ywolberg/DDD_Africa/ddd-cnv/ddd-cnv/bin/xhmm_xcnv_to_vcf2.py -i /dataG/ddd/analysis/CNV_calling/TRUTHSET.2023.08.22/XHMM_filtered.xcnv -d out_xhmm -c XHMM_CNVs_filtered.vcf 
    """
}

// Index VCF
process index_vcf {
    // Tags the process as 'index_vcf'
    tag { 'index_vcf' }

    // Labels the process as 'xhmm'    
    label 'xhmm'

    // Specifies destination directory where output files are published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite=true

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
    bgzip -c out_xhmm/XHMM_CNVs_filtered.vcf > out_xhmm/XHMM_CNVs_filtered.vcf.gz
    
    # Sort the VCF file with BCFtools
    bcftools sort -o sorted_XHMM_CNVs_filtered.vcf.gz -O z out_xhmm/XHMM_CNVs_filtered.vcf.gz 

    # Index the VCF file with Tabix
    tabix -p vcf sorted_XHMM_CNVs_filtered.vcf.gz

    # Unzip the sorted VCF file and redirects output to VCF file
    zcat sorted_XHMM_CNVs_filtered.vcf.gz | less -S > sorted_XHMM_CNVs_filtered.vcf
    """
}
