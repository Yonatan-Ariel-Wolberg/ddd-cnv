#!/usr/bin/env nextflow

// DOMAIN SPECIFIC LANGUAGE (DSL) VERSION 2
nextflow.enable.dsl=2

// INPUT FILES FOR CANOES & XHMM
// Creates a channel named bams from BAM files and index files (BAI)
Channel.fromFilePairs([params.bams + '/*{.bam,.bam.bai}'])
    .map { it -> [ it[0], it[1].find { it =~ '.bam$'}, it[1].find { it =~ '.bai$'}] } // Maps each pair to tuple with a BAM and BAI file
    .set { bams } // Set channel name to 'bams'

// DEFINE AND CREATE OUTPUT DIRECTORY
outdir = file(params.outdir, type: 'dir')
outdir.mkdir()

/* Include workflow modules for separate scripts for each tool:
* CANOES, XHMM, CLAMMS
*/

workflow = params.workflow

// INCLUDE MODULES FOR DIFFERENT FUNCTIONS OF EACH TOOL
// INCLUDE CANOES MODULE
include { genReadCounts; calcGC_CANOES; runCANOES; filterCANOESCNVs } from './modules/cnv_calling/modules-canoes.nf'

// INCLUDE XHMM MODULE
include { groupBAMs; gatkDOC; combineDOC; calcGC_XHMM; filterSamples; runPCA;
        normalisePCA; filterZScore; filterRD; discoverCNVs; genotypeCNVs; filterXHMMCNVs } from './modules/cnv_calling/modules-xhmm.nf'

// INCLUDE CLAMMS MODULE
include { generateWindows; samtoolsDOC; normalizeDOC; trainModels; callCNVs; filterCLAMMSCNVs } from './modules/cnv_calling/modules-clamms.nf'

// CANOES WORKFLOW
workflow RUN_CANOES {
    
    take:
    bams

    main:
    // Process BAM files to generate read counts
    genReadCounts(bams.collectFile () { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] })

    // Calculate GC content
    calcGC_CANOES()

    // Call CNVs with CANOES
    runCANOES(genReadCounts.out.canoes_reads, calcGC_CANOES.out.gc_content)

    // Filter resulting CNVs
    filterCANOESCNVs(runCANOES.out.cnvs)
}

// XHMM WORKFLOW
workflow RUN_XHMM {
    
    take:
    bams

    main:
    // Group BAM files
    groupBAMs(bams.collectFile () { item -> [ 'bam_list_unsorted.txt', "${item[1]}" + '\n' ] })

    // Generate DOC metrics
    gatkDOC(groupBAMs.out.bam_groups.flatMap().map { it -> [it.name[0..-6], it] })

    // Combine DOC metrics
    combineDOC(gatkDOC.out.bam_group_doc.collect { it -> it[1] })

    // Calculate GC content
    calcGC_XHMM()

    // Filter samples
    filterSamples(combineDOC.out.combined_doc, calcGC_XHMM.out.extreme_gc_targets)

    // Run PCA
    runPCA(filterSamples.out.filtered_centered)

    // Normalise PCA
    normalisePCA(filterSamples.out.filtered_centered, runPCA.out.pca_data)

    // Filter by Z Score
    filterZScore(normalisePCA.out.data_pca_norm)

    // Filter by Read Depth
    filterRD(combineDOC.out.combined_doc,
        filterSamples.out.excluded_filtered_targets, filterSamples.out.excluded_filtered_samples,
        filterZScore.out.excluded_zscore_targets, filterZScore.out.excluded_zscore_samples )

    // Call CNVs with XHMM
    discoverCNVs(filterRD.out.orig_filtered, filterZScore.out.pca_norm_zscore)

    // Genotype CNVs
    genotypeCNVs(filterRD.out.orig_filtered, filterZScore.out.pca_norm_zscore, discoverCNVs.out.cnvs)

    // Filter the final set of CNVs
    filterXHMMCNVs(discoverCNVs.out.cnvs.collect())
}

// CLAMMS WORKFLOW
workflow RUN_CLAMMS {
    take:
    bams

    main:
    // Process BAMs to generate windows
    generateWindows()

    // Use SAMtools to calculate DOC
    samtoolsDOC(bams, generateWindows.out.windows)

    // Normalises DOC
    normalizeDOC(samtoolsDOC.out.coverage, generateWindows.out.windows)

    // Training the statistical model (Lattice-Aligned Mixture Model)
    trainModels(normalizeDOC.out.coverage_norm.collect(), generateWindows.out.windows)

    // Call CNVs with CLAMMS
    callCNVs(normalizeDOC.out.coverage_norm_set, trainModels.out.models)

    // Filter the final set of CNVs
    filterCLAMMSCNVs(callCNVs.out.cnvs.map { it -> it[1] }.collect())
}

// Pick and Choose
workflow {
    // Workflow that can switch between using CANOES, XHMM and CLAMMS
    switch (workflow) {
        // Run CANOES
        case['canoes']:
        RUN_CANOES(bams)
        break
        // =====

        // Run XHMM
        case['xhmm']:
        RUN_XHMM(bams)
        break
        // =====

        // Run CLAMMS
        case['clamms']:
        RUN_CLAMMS(bams)
        break
        // =====

        // Run ALL
        case['all']:
        RUN_CANOES(bams)
        RUN_XHMM(bams)
        RUN_CLAMMS(bams)
        // =====

        default:
        exit 1, """
        OOPS! SEEMS LIKE WE HAVE A WORKFLOW ERROR!

        No workflow \'mode\' give! Please use one of the following options for workflows:
        --mode canoes // To run CANOES workflow
        --mode xhmm // To run the XHMM workflow
        --mode clamms // To run the CLAMMS workflow
        --mode all // To run all workflows
        """

        break
    }
}