#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// DEFINE PARAMETERS
ref = file(params.ref, type: 'file')
probes = file(params.ref, type: 'file')
outdir = file(params.outdir, type: 'dir')

// PROCESS THAT GENERATES READ COUNTS
process genReadCounts {
    // Tags the process as 'read_counts'
    tag { "read_counts" }

    // Labels the process as 'canoes'
    label 'canoes'

    // Specifies destination directory where output files are published
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true

    // Loads the 'bedtools/2.30.0' module
    module 'bedtools/2.30.0'

    // Specifies input
    input:
    path(bam)

    // Specifies output
    tuple val("canoes_in"), path("canoes.reads.txt"), path("sample_list"), emit: canoes_reads

    // Specifies script
    script:
    """
    # Sorts list of BAM files
    sort bam_list_unsorted.txt > bam_list_sorted.txt

    # Counts reads from BAM files list, intersecting them with regions defined in a BED file, then filters reads with mapping quality (q) >= 20
    bedtools mutlicov -bams 'cat bam_list_sorted.txt' -bed ${probes} -q 20 > canoes.reads.txt

    # Reads each line from 'bam_list_sorted.txt'
    while read line
    do

        # Extracts base filename (without directory path) of each BAM file and removes '.bam' extension
        echo `basename \$line` | sed 's/.bam//'
    
    # Writes modified names to a file name 'sample_list'
    done < bam_list_sorted.txt > sample_list
    """
}

// PROCESS THAT CALCULATES GC CONTENT
process calcGC_CANOES {

    // Tags process as 'calc_gc'
    tag { " calc_gc " }

    // Labels process as 'canoes'
    lable 'canoes'

    // Specifies destination directory where output files are published
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true

    // Specifies output
    output:
    path("gc.txt"), emit: gc_content

    // Specifies script
    script:
    """
    # Specifies Java with maximum heap size of 2000 MB, sets the temporary directory for Java operations
    # Executes GATK toolkit
    java -Xmx2000m -Djava.io.tmpdir=TEMP -jar /home/phelelani/applications/gatk-2.1-9/GenomeAnalysisTK.jar \
        -T GCContentByInterval \ # Specifies GATK tool for calculating GC content
        -L ${probes} \ # Specifies BED file containing intervals for GC content calculation
        -R ${ref} \ # Specifies reference genome used for alignment
        -o gc.txt # Specifies output
    """
}

// CALL CNVS WITH CANOES
process runCANOES {

    // Tags process as 'run_canoes'
    tag { "run_canoes" }

    // Labels process as 'canoes'
    label 'canoes'

    // Specifies destination directory where output files are published
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true

    // Specifies input
    tuple val(canoes_in), path(canoes_reads), path(sample_list)
    path(gc_content)

    // Specifies output
    path("*"), emit: canoes_out
    path("Sample_CNVs.csv"), emit: cnvs

    // Specifies script
    script:
    """
    # Removes occurrences of 'chr' globally, then deletes files starting with 'X' and 'Y'
    sed 's/chr//g; /^X/d; /^Y/d' ${canoes_reads} > canoes.reads_new.txt

    # Removes occurrences of 'chr' globally, then deletes files starting with 'X' and 'Y'
    sed 's/chr//g; /^X/d; /^Y/d' ${gc_content} > gc_new.txt

    # Runs CANOES
    run_canoes.R gc_new.txt canoes.reads_new.txt ${sample_list}
    """
}

// FILTER CNV CALLS
process filterCANOESCNVs {
    
    // Tags process as 'filter_cnvs'
    tag { "filter_cnvs" }

    // Labels process as 'canoes'

    // Specifies destination directory where output files are published
    publishDir "${outdir}/out_CANOES", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(cnvs)

    // Specifies output
    output:
    path("Sample_CNVs_filtered.csv"), emit: filtered_cnvs

    // Specifies script
    script:
    """
    # Selects rows where column 10 is greater than 80 and not equal to 'NA' and where column 4 is greater than or equal to 100
    awk '{ if( (\$10>=80) && (\$10!="NA") && (\$4>=100) ) { print }}' Sample_CNVs.csv > Sample_CNVs_filtered.csv
    """
}   