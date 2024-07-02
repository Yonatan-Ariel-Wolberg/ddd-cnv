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
    # Calls XHMM tool
    # --matrix specifies operation will be on the matrix of read depths
    # -r specifies the input file
    # --centerData centers the data to adjust for biases
    # --centerType target specifies that centering should be done based on targets
    xhmm --matrix -r DATA.RD.txt --centerData --centerType target \
        -o DATA.filtered_centered.RD.txt \   # Specifies output file
        --outputExcludedTargets DATA.filtered_centered.RD.txt.filtered_targets.txt \   # Specifies output file for list of excluded targets
        --outputExcludedSamples DATA.filtered_centered.RD.txt.filtered_samples.txt \   # Specifies output file for list of excluded samples
        --excludeTargets extreme_gc_targets.txt \   # Excludeds targets listed in the given file from analysis
        --minTargetSize 10 --maxTargetSize 10000 \   # Specifies min and max target sizes
        --minMeanTargetRD 10 --maxMeanTargetRD 500 \   # Specifies min and max read depths for targets
        --minMeanSampleRD 25 --maxMeanSampleRD 200 \   # Specifies min and max read depths for samples
        --maxSdSampleRD 150   # Specifies the max standard deviation of read depth for samples
    """
}

// RUNS PCA ON MEAN_CENTERED DATA:
process runPCA {
    // Tags process as 'run_pca'
    tag { 'run_pca' }

    // Labels process as 'xhmm'
    label 'xhmm'

    // Specifies input
    input:
    path(filtered_centered)

    // Specifies output
    output:
    tuple val("pca_data"), path("DATA.RD_PCA*"), emit: pca_data

    // Specifies script
    script:
    """
    # xhmm --PCA runs PCA operation by XHMM tool
    # -r  specifies input data
    # --PCAfiles specifies the prefix for the output files related to PCA
    xhmm --PCA -r DATA.filtered_centered.RD.txt --PCAfiles DATA.RD_PCA
    """
}

// NORMALIZES MEAN-CENTERED DATA USING PCA INFORMATION:
process normalisePCA {
    // Tags process as 'norm_pca'
    tag { 'norm_pca' }

    // Label process as 'xhmm'
    label 'xhmm'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(filtered_centered)
    tuple val(pca), path(pca_data)

    // Specifies output
    output:
    path("DATA.PCA_normalized.txt"), emit: data_pca_norm

    // Specifies script
    script:
    """
    # xhmm --normalize specifies that the operation to be performed by xhmm is normalization
    # -r specifies the input data
    xhmm --normalize -r DATA.filtered_centered.RD.txt \
        --PCAfiles DATA.RD_PCA \   # Specifies the prefix for the PCA output files generated in the previous PCA step
        --normalizeOutput DATA.PCA_normalized.txt \   # Specifies the output file for the normalzied read depth data
        --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7 # Specifies normalization method to be used (Proportion of Variance Explained mean) and the factor to use for said method (0.7)
    """
}

// FILTERS AND Z_SCORE CENTERS (BY SAMPLE) THE PCA-NORMALIZED DATA:
process filterZScore {
    // Tags process as 'filter_zscore'
    tag { 'filter_zscore' }

    // Labels process as 'xhmm'
    label 'xhmm'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(data_pca_norm)

    // Specifies output
    output:
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt"), emit: pca_norm_zscore
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt"), emit: excluded_zscore_targets
    path("DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt"), emit: excluded_zscore_samples

    // Specifies script
    script:
    """
    xhmm --matrix -r DATA.PCA_normalized.txt   #  Calls xhmm tool and specifies that operation is matrix of read depths
        --centerData --centerType sample --zScoreData \   # Centers data based on samples  and Z-scores the data to standardize it
        -o DATA.PCA_normalized.filtered.sample_zscores.RD.txt \   # Specifies output file
        --outputExcludedTargets DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \   # Specifies output for list of excluded targets
        --outputExcludedSamples DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \   # Specifies output for list of excluded samples
        --maxSdTargetRD 30   # Excludes targets with a standard deviation greater than 30
    """
}

// FILTERS ORIGINAL READ-DEPTH DATA TO BE THE SAME AS FILTERED, NORMALIZED DATA:
process filterRD {
    // Tags process as 'filter_rd'
    tag { 'filter_rd' }

    // Labels process as 'xhmm'
    label 'xhmm'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(combined_doc)
    path(excluded_filtered_targets)
    path(excluded_zscore_targets)
    path(excluded_filtered_samples)
    path(excluded_zscore_samples)

    // Specifies output
    output:
    path("DATA.same_filtered.RD.txt"), emit: orig_filtered

    // Specifies script
    script:
    """
    xhmm --matrix -r DATA.RD.txt \   # Calls XHMM and specifies tool will be matrix of read depths
        --excludeTargets DATA.filtered_centered.RD.txt.filtered_targets.txt \   # Excludes targets listed in this file
        --excludeTargets DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \   # Excludes targets listed in this file
        --excludeSamples DATA.filtered_centered.RD.txt.filtered_samples.txt \   # Excludes samples listed in this file
        --excludeSamples DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \   # Excludes samples listed in this file
        -o DATA.same_filtered.RD.txt   # Specifies the output file
    """
}

// DISCOVERS CNVS IN NORMALIZED DATA:
process discoverCNVs {
    // Tags process as 'discover_cnvs'
    tag { 'discover_cnvs' }

    // Labels process as 'xhmm'
    label 'xhmm'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(orig_filtered)
    path(pca_norm_zscore)

    // Specifies output
    output:
    path("*"), emit: cnvs
    
    script:
    """
    xhmm --discover -p ${xhmm_conf} \   # Calls XHMM to discover CNVs with the required configuration file
        -r DATA.PCA_normalized.filtered.sample_zscores.RD.txt \   # Specifies input file containing the PCA-normalized and z-scored read depth data
        -R DATA.same_filtered.RD.txt \   # Specifies input file containing the same filtered read depth data
        -c DATA.xcnvv -a DATA.aux_xcnv -s DATA   # Specifies output files for the discovered CNVs, auxiliary CNV info and additional data files generated during the CNV discovery process with the prefix 'DATA'
    """
}

// GENOTYPES DISCOVERED CNVS IN ALL SAMPLES:
process genotypeCNVs {
    // Tags process as 'genotype_cnvs'
    tag { 'genotype_cnvs' }
    
    // Labels process as 'xhmm'
    label 'xhmm'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(orig_filtered)
    path(pca_norm_zscore)
    path(cnvs)

    // SPecifies output
    output:
    path("*"), emit: genotypes

    // Specifies output
    script:
    """
    xhmm --genotype -p ${xhmm_conf} \   # Calls XHMM to genotype CNVs with the required configuration file
        -r DATA.PCA_normalized.filtered.sample_zscores.RD.txt \   # Specifies input file containing the PCA_normalized and z_scored read depth data
        -R DATA.same_filtered.RD.txt \   # Specifies input file containing the same filtered read depth data
        -g DATA.xcnv -F ${ref} \   # Specifies input file containing the discovered CNVs to be genotyped and the reference genome
        -v DATA.vcf   # Specifies output file in VCF format containing the genotyped CNVs
    """
}

// FILTERS CNVS CALLED BY XHMM
process filterXHMMCNVs {
    // Tags process as 'filter_cnvs'
    tag { 'filter_cnvs' }

    // Labels process as 'xhmm'
    label 'xhmm'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_XHMM", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(cnvs)

    // Specifies output
    output:
    path("DATA_filtered.xcnv"), emit: filtered_cnvs

    // Specifies script
    script:
    """
    # Filters DATA.xcnv files that have Q_SOME = or > 60, Q_NON_DIPLOID = or > 60, and Q_START = or > 60
    awk '{ if( (\$9>=60) && (\$10>=60) && (\$11>=60) ) { print } }' DATA.xcnv > DATA_filtered.xcnv
    """
}