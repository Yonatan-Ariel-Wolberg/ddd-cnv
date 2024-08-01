#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// INPUTS
params.csv_input = '/dataG/ddd/analysis/CNV_calling/TRUTHSET.2023.08.22/CANOES_CNVs_filtered.csv'
params.sample_list = '/home/ywolberg/Yoni_DDD/EGA_sample_ids.txt'
params.canoes_to_vcf = '/home/ywolberg/Yoni_DDD/ddd-cnv/bin/canoes_csv_to_vcf.py'
params.outdir = '/home/ywolberg/Yoni_DDD'




// PROCESS THAT CONVERTS CANOES OUTPUT FILE TO VCF FORMAT
process CANOESToVCF {
    // Tags process as 'canoes_to_vcf'
    tag { 'canoes_to_vcf' }

    // Labels process as 'canoes'
    label 'canoes'

    // Specifies destination directory where output file is to be published
    publishDir "${params.outdir}/out_CANOES", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(params.csv_input)

    // Specifies output
    output:
    file("*.vcf") into canoes_vcf

    // Specifies script
    script:
    """
    python ${params.canoes_to_vcf} \
        -i ${params.csv_input} \
        -d ${params.outdir}/out_CANOES \
        -c ${params.canoes_vcf} \
        -s ${params.sample_list} 
    """
}

process IndexVCF {
    // Tags process as 'IndexVCF'
    tag { 'IndexVCF' }

    // Labels process as 'canoes'
    label 'canoes'

    // Specifies destination directory where output file is to be published
    publishDir "${params.outdir}/out_CANOES", mode: 'copy', overwrite: true

    // Loads the 'samtools/1.20' module
    module 'samtools/1.20'

    // Specifies input
    input:
    file("CANOES_CNVs_filtered.vcf") from canoes_vcf

    // Specifies output
    output:
    tuple file("CANOES_CNVs_filtered.vcf.gz"), file("sorted_CANOES_CNVs_filtered.vcf.gz"), file("sorted_CANOES_CNVs_filtered.vcf.gz.tbi"), file("sorted_CANOES_CNVs_filtered.vcf") into indexed_vcf
    
    // Specifies script
    script:
    """
    # Compresses the VCF file
    bgzip -c CANOES_CNVs_filtered.vcf > CANOES_CNVs_filtered.vcf.gz

    # Sorts the VCF file
    bcftools sort \
        -o sorted_CANOES_CNVs_filtered.vcf.gz \
        -O z \
        CANOES_CNVs_filtered.vcf.gz

    # Indexes the VCF file
    tabix -p vcf CANOES_CNVs_filtered.vcf.gz

    # Decompresses and writes out sorted .vcf.gz file to a .vcf file
    bcftools view sorted_CANOES_CNVs_filtered.vcf.gz | less -S > CANOES_CNVs_filtered.vcf
    """
}

workflow CANOESIndex {
    canoes_vcf = CANOESToVCF()
    indexed_vcf = IndexVCF(canoes_vcf)
}