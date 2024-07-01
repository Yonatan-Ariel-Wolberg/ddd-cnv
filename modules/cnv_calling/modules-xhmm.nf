#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// INPUT REQUIREMENTS
// Reference Genome
ref         = file(params.ref, type: 'file')

// Probes
probes      = file(params.probes, type: 'file')

// XHMM Config
xhmm_conf   = file(params.xhmm_conf, type: 'file')

// Output Directory
outdir = file(params.outdir, type: 'dir')


// SPLIT BAM FILES INTO GROUPS
process groupBAMs {
    // Tags process as 'group_bams'
    tag { 'group_bams' }

    // Specifies input
    input:
    path(bams)

    // Specifies output
    output:
    path("bam_group_*"), emit: bam_groups

    // Specifies script:
    script:
    """
    # Sorts list of bam files
    sort bam_list_unsorted.txt > bam_list_sorted.txt

    # Split bam_sorted_list.txt into files (each of 5 lines) with nmeric suffixes
    split -l 5 bam_list_sorted.txt --numeric-suffixes --additional-suffix=.list bam_group_
    """
}

// RUN GATK FOR DEPTH OF COVERAGE (FOR SAMPLES IN EACH GROUP)
process gatkDOC {
    // Tags process as 'calc_doc'
    tag { 'calc_doc' }

    // Specifies that 10 tasks can be run concurrently
    maxForks 10

    // Assigns 11 GB of memory to the process
    memory '11 GB'

    // Loads gatk/4.2.5.0 module
    module 'gatk/4.2.5.0'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    tuple val(group), path(list)

    // Specifies output
    output:
    tuple val(group), path("${group}.DATA.sample_interval_summary"), emit: bam_group_doc

    // Specifies script
    script:
    """
    # Launches GATK tool with 10 GB of memory
    gatk --java-options "-Xmx10G" \
        DepthOfCoverage \   # Use DepthOfCoverage tool
        -I ${list} \   # Specify input BAM files
        -L ${probes} \   # Specify the intervals to restrict the analysis to (i.e. target capture regions)
        -R ${ref} \   # Specify the reference genome
        --max-depth-per-sample 5000 \   # Sets maximum DOC per sample to 5000 reads
        --verbosity INFO \   # Sets verbosity level of logging to 'INFO', providing detailed info about execution
        --omit-depth-output-at-each-base true \   # Omits DOC output at each base, reducing output size
        --omit-locus-table true \   # Omits locus table, focusing only on DOC metrics
        --min-base-quality 0 \   # Specifies minimum base quality threshold for base calls
        --read-filter MappingQualityReadFilter \ # Applies a read filter based on mapping quality
        --minimum-mapping-quality 20 \   # Sets the minimum mapping quality threshold for reads to be included in analysis
        --start 1 --stop 5000 --nBins 200 \   # Defines the range and number of bins for DOC calculations
        --include-ref-n-sites true \   # Includes reference N sites in the coverage analysis
        --count-type COUNT_READS \   # Specifies the type of counting for depth calculations
        --output ${group}.DATA   # Specifies output file
        """
}

// PROCESS THAT COMBINES GATK DEPTH_OF_COVERAGE OUTPUTS FOR MULTIPLE SAMPLES (AT SAME LOCI):
process combineDOC {
    // Tags process as 'combine_doc'
    tag { 'combine_doc' }

    // Labels process as 'xhmm'
    label 'xhmm'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(list)

    // Specifies output
    output:
    path("DATA.RD.txt"), emit: combined_doc

    // Specifies script
    script:
    """
    # Replace commas with tabs in files
    for i in *.sample_interval_summary;
    do
        sed 's/,/   /g' \$i > \${i%.sample_interval_summary}.fixed_sample_interval_summary;
    done

    # List fixed files
    ls *.fixed_sample_interval_summary > the_list

    # Merge GATK depths with XHMM
    xhmm --mergeGATKdepths -o DATA.RD.txt --GATKdepthsList the_list
    """
}

// OPTIONALLY, RUN GATK TO CALCULATE THE PER_TARGET GC CONTENT AND CREATE A LIST OF THE TARGETS WITH EXTREME GC CONTENT:
process calcGC_XHMM {
    // Tags process as 'calc-dc'
    tag { 'calc_gc' }

    // Specifies the destination directory where output file are to be published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies output
    output:
    path("DATA.locus_GC.txt"), emit: gc_content
    path("extreme_gc_targets.txt"), emit: extreme_gc_targets

    // Specifies script
    script:
    """
    # Executes Java with 3072 MB of memory and a temporary directory, running the GATK JAR file
    java -Xmx3072m -Djava.io.tmpdir=TEMP -jar /home/phelelani/applications/gatk-2.1-9/GenomeAnalysisTK.jar \
        -T GCContentByInterval \   # Specifies the GATK tool to be 'GCContentByInterval', which caluclates GC content by genomic intervals
        -L ${probes} \   # Specifies the intervals (based on target regions) to compute GC content for
        -R ${ref} \   # Specifies the reference genome
        -o DATA.locus_GC.txt   # Specifies the output file

    # Filter based on GC content, with a GC content of under 0.1 or over 0.9 being considered extreme
    cat DATA.locus_GC.txt | awk '{ if (\$2 < 0.1 || \$2 > 0.9) print \$1 }' > extreme_gc_targets.txt
    """
}

// FILTERS SAMPLES AND TARGETS AND THEN MEAN_CENTERS THE TARGETS:
process filterSamples {
    // Tags process as 'filter_samples'
    tag { 'filter_samples' }

    // Labels process as 'xhmm'
    label 'xhmm'

    // Specifies the destination directory where output files are to be published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(combined_doc)
    path(extreme_gc_targets)

    // Specifies outputs
    output:
    path("DATA.filtered_centered.RD.txt"), emit: filtered_centered
    path("DATA.filtered_centered.RD.txt.filtered_targets.txt"), emit: excluded_filtered_targets
    path("DATA.filtered_centered.RD.txt.filtered_samples.txt"), emit: excluded_filtered_samples

    // Specifies script
    script:
    """
    xhmm --matrix -r DATA.RD.txt --centerData --centerType target \
        -o DATA.filtered_centered.RD.txt \
        --outputExcludedTargets DATA.filtered_centered.RD.txt.filtered_targets.txt
    """
}