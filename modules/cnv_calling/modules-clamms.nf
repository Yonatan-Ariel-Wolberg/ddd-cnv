#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// INPUT PARAMETERS
// Reference Genome
ref         = file(params.ref, type: 'file')

// Probes
probes      = file(params.probes, type: 'file')

// Mappability
mappability = file(params.mappability, type: 'file')

// Special Regions
special_reg = file(params.special_reg, type: 'file')

// Output Directory
outdir      = file(params.outdir, type: 'dir')

// Sex Information
sexinfo     = file(params.sexinfo, type: 'file')


// PROCESS THAT GENERATES WINDOWS
process generateWindows {
    // Tags process as 'generate_windows'
    tag { 'generate_windows' }

    // Assigns 11 GB of memory to process
    memory '11 GB'

    // Loads 'bedtools/2.30.0' module
    module 'bedtools/2.30.0'

    // Specifies destination directory where output files are published
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    // Specifies output
    output:
    path("windows.bed"), emit: windows

    // Specifies script
    script:
    """
    # Sets environmental variable INSERT_SIZE to 200
    export INSERT_SIZE=200

    # Sets environmental variable CLAMMS_DIR to a specified directory path
    export CLAMMS_DIR=/home/phelelani/applications/clamms

    # Sorts 'probes' file by first column (chromosome) and then by second column (starting position) numerically
    sort -k1,1 -k2,2n ${probes} > targets_sorted.bed

    # Runs a script to annotate windows
    \$CLAMMS_DIR/annotate_windows.sh targets_sorted.bed ${ref} ${mappability} \
        \$INSERT_SIZE ${special_reg} > windows.bed
    """
}

// PROCESS THAT USES SAMTOOLS TO CALCULATE DOC
process samtoolsDOC {
    // Tags process with 'sample' variable (sample ID)
    tag { sample }

    // Sets maximum number of concurrent tasks to 20
    maxForks 20

    // Assigns 11 GB of memory to process
    memory '11 GB'

    // Loads 'samtools/1.15' module
    module 'samtools/1.15'

    // Specifies destination directory where output files are published
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    // Specifies input
    input:
    tuple val(sample), path(bam), path(bai)
    path(windows)

    // Specifies output
    tuple val(sample), path("${sample}.coverage.bed"), emit: coverage

    // Specifies script
    script:
    """
    # Calculates depth of coverage of regions specified in windows.bed file
    # -Q 30 sets minimum mapping quality threshold of 30
    # Processes output of samtools bedcov using awk
    # \$1: first field (chromosome)
    # \$2: second field (start position)
    # \$3: third field (end position)
    # \$NF: last filed (sum of coverage across the region)
    # \$NF/(\$3-\$2): calculates the average coverage across the region by dividing the sum of coverage by the length of the region
    
    samtools bedcov -Q 30 ${windows} ${bam} | awk '{ printf \"%s\\t%d\\t%d\\t%.6g\\n\", \$1, \$2, \$3, \$NF/(\$3-\$2); }' > ${sample}.coverage.bed
    """
}

// PROCESS THAT NORMALIZES DOC
process normalizeDOC {
    // Tags process with variable 'sample' (sample ID)
    tag { sample }

    // Assigns memory of 11 GB to process
    memory '11 GB'

    // Specifies destination directory where output files are published
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    // Specifies input
    input:
    tuple val(sample), path(coverage)
    path(windows)

    // Specifies output
    path("${sample}.norm.cov.bed"), emit: coverage_norm

    // Specifies script
    script:
    """
    # Sets environmental variable 'CLAMMS_DIR' to the specified directory path
    export CLAMMS_DIR=/home/phelelani/applications/clamms

    # Runs normalize_coverage script to normalize coverage data and uses sed to remove 'chr' prefix from chromosome names in the output
    $CLAMMS_DIR/normalize_coverage ${sample}.coverage.bed ${windows} | sed 's/^chr//g' > ${sample}.norm.cov.bed
    """
}

// PROCESS THAT TRAINS LATTICE_ALIGNED MIXTURE MODELS
process trainModel {
    // Tags process as 'train_models'
    tag { 'train_models' }

    // Assigns 11 GB of memory to the process
    memory '11 GB'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(coverage_norm)
    path(windows)

    // Specifies output
    output:
    path("models.bed"), emit: models

    // Specifies script
    script:
    """
    # Sets environmental variable CLAMMS_DIR to a specified directory path
    export CLAMMS_DIR=/home/phelelani/applications/clamms
    
    # Removes 'chr' prefix globally from chromosome names in windows.bed file
    sed 's/^chr//g' ${windows} > windows_new

    # Runs the fit_models script to train Lattice-Aligned Mixture Models
    # sed removes 'chr' prefixes from output
    $CLAMMS_DIR/fit_models ${sexinfo} windows_new | sed 's/^chr//g' > models.bed
    """
}

// CALL CNVS WITH CLAMMS
process callCNVs {
    // Tags process with variable 'sample' (sample ID)
    tag { sample }

    // Assigns 11 GB of memory to process
    memory '11 GB'

    // Instructs Nextflow to continue executing subsequent tasks even if errors occur with the current process
    errorStrategy 'ignore'

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    // Specifies input
    input:
    tuple val(sample), path(norm_cov)
    path(models)

    // Specifies output
    output:
    tuple val(sample), path("${sample}.cnv.bed"), emit: cnvs

    // Specifies script
    script:
    """
    # Sets the environmental variable CLAMMS_DIR to a specified directory path
    export CLAMMS_DIR=/home/phelelani/applications/clamms

    # Searches for sample ID in sex info file
    # cut extracts second field (sex info) from output of grep
    sex =`grep ${sample} ${sexinfo} | cut -f 2`

    # Executes call_cnv script to call CNVs 
    $CLAMMS_DIR/call_cnv ${norm_cov} ${model} --sex \$sex > ${sample}.cnv.bed
    """
}

// PROCESS THAT FILTERS CNV CALLS
process filterCLAMMSCNVs {
    // Tags process as 'filter_cnvs'
    tag { 'filter_cnvs' }

    // Specifies destination directory where output files are to be published
    publishDir "${outdir}/out_CLAMMS", mode: 'copy', overwrite: true

    // Specifies input
    input:
    path(cnvs)

    // Specifies output
    output:
    tuple path("samples.cnv.bed"), path("samples.cnv.filtered.bed"), emit: filtered_cnvs

    // Specifies script
    script:
    """
    # Writes out each .cnv.bed file and collates them in a single file
    cat *.cnv.bed > samples.cnv.bed

    # Prints lines where Q_SOME > 500 and Q_EXACT > 0
    awk '{ if( (\$9>=500) && (\$10>0) ) { print }}' samples.cnv.bed > samples.cnv.filtered.bed
    """
}