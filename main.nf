#!/usr/bin/env nextflow
/*
===========================================================================================
Pipeline for pangenome intra-diversity analysis
===========================================================================================
Entrypoint 1. Raw fastq reads.
Optional stop A: After assembly and binning.
Entrypoint 2. Pre-existing bins.
Optional stop B: after constructing pangenomes.
Entrypoint 3. A reference genome in fasta format.
Completion: Reference/pangenome with SNVs in core and some population genomic calculations.
-------------------------------------------------------------------------------------------
*/

include { validateParameters ; paramsHelp ; paramsSummaryLog ; samplesheetToList } from 'plugin/nf-schema'

// import subworkflow modules
//include {  } from './modules/subworkflows/'
include { subsample_reads } from './modules/subworkflows/subsample_reads'
include { raw_to_bins } from './modules/subworkflows/raw_to_bins'
include { provided_bins } from './modules/subworkflows/provided_bins'
include { pangenome_assembly } from './modules/subworkflows/pangenome_assembly'
include { match_samps_to_pang } from './modules/subworkflows/match_samps_to_pang'
include { variant_calling } from './modules/subworkflows/variant_calling'

workflow {
    // Validate input parameters
    validateParameters(parameters_schema: 'nextflow_schema.json', monochrome_logs: true)

    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    //Check project parameter
    def badChars = ["^", "(", ")", "+", " ", "|"]
    //findAll goes through each character (pchar) of the project name
    //pchar is checked against each element (bchar) in the badChars list
    // any pchar that matches a bchar will be added to a list which is returned. No bad pchar = no list = false/no exception.
    if (params.project.findAll { pchar -> badChars.any { bchar -> pchar.contains(bchar) } }) {
        throw new Exception("Invalid project name. Special characters and whitespaces not allowed.")
    }

    //Gives a warning if project already exists
    if (!workflow.resume) {
        //Workflow was not resumed, checking project dir
        def Path projDir = new File(params.project).toPath()
        if (projDir.exists() == true) {
            throw new Exception("Project directory ${params.project} already exists, choose a new name or use the -resume flag. WARNING: Note that if you resume the wrong job, this might overwrite previous results.")
        }
    }

    //Check that subsampling conditions are correct
    println("The subsample parameter is set to ${params.subsample}")
    if (!params.subsample) {
        if (params.ref_genome == null) {
            if (params.bins != null) {
                if (params.readcount == null) {
                    throw new Exception("When skipping subsampling from bins entry it is assumed that you already have subsampled reads.\nTherefore a tab delimited readcount file with Sample, Nr_fastqs, Total_reads needs to be provided with --readcount <path/to/file>")
                } 
            }    
            else {
                throw new Exception("You cannot skip subsampling when running assembly and binning. For options run nextflow main.nf --help.")
            }
        }
    }

    if (params.bins != null && params.ref_genome != null) {
        throw new Exception("You can either provide pre-assembled bins or previously created reference genomes, but not both. Please use either the --bins flag or the --ref_genome flag. ")
    }

    if (params.only_bins && (params.bins != null || params.ref_genome != null)) {
        throw new Exception("Stopping after binning is only allowed for Entrypoint 1: raw reads.\nEither set --only_bins false or remove --bins/--ref_genome input.")
    }

    println("Starting. Your results will be published at ${params.project}.")

    //Subsampling
    sam_ch = channel.fromPath(params.samples, type: "file", checkIfExists: true)
    fastq_ch = channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
    subsample_reads(sam_ch, fastq_ch, params.project, params.subsample, params.nr_subsamp, params.readcount)
    readcounts = subsample_reads.out.readcounts 
    sub_reads = subsample_reads.out.sub_reads
    individual_samp_files = subsample_reads.out.individual_samp_files

    /*
    If the user provided a dir with reference genomes, the pipeline will only run 
    the map_and_detect_variants workflow. This is the third entrypoint.
    */
    if (params.ref_genome != null) {
        channel.fromPath("${params.ref_genome}", type: "file", checkIfExists: true).multiMap { ch -> core: NBPs: ch }.set { ref_gen }
        println("Entrypoint 3: Reference genome.")
        /*
        When using a reference genome we don't have core and consensus,
        therefore handling the reference as both.
        This means that the whole genome is used both for mapping a subset of the reads,
        and for the variance analysis. Contigs of interest can be provided in a separate file to avoid this.
        */
        core_ch = ref_gen.core
        NBPs_ch = ref_gen.NBPs

        if (params.contigs != null) {
            contigs_ch = channel.fromPath("${params.contigs}", type: "file", checkIfExists: true)
        }
        else {
            contigs_ch = null
        }
    } else {
        if (params.bins != null) {
            println("Entrypoint 2: Provided bins.")
            //If bins were provided we don't need to do assembly
            //Optional inputs can be null
            provided_bins(params.samples, params.bins, params.fastq, params.project, params.bintables)
            bins_ch = provided_bins.out.bins
            bintable_ch = provided_bins.out.bintable
        } else {
            println("Entrypoint 1: Raw reads to assemble and bin.")
            /*
    	    Runs assembly and binning.
    	    */
            raw_to_bins(params.fastq, individual_samp_files, params.project, params.binners, params.contig_len)
            bins_ch = raw_to_bins.out.bins
            bintable_ch = raw_to_bins.out.bintable
        } 
        /*
    	This workflow will cluster bins, create pangenomes, and send out core and NBPs fasta files for the pangenomes.
        Done for first and second entrypoints (raw and existing bins).
    	*/
        if (!params.only_bins) {
            pangenome_assembly(bins_ch, bintable_ch, params.project, params.MAGcomplete, params.MAGcontam, params.min_mOTU_MAGs,
                               params.min_contig_len, params.taxSort)

            core_ch = pangenome_assembly.out.core_fasta
            NBPs_ch = pangenome_assembly.out.NBPs_fasta
            contigs_ch = pangenome_assembly.out.contigs_tsv

            
        }
    }
    //If we have pangenomes or a reference genome and have subsampled reads we estimate coverages to determine which samples our genomes are present in.
    if (!params.only_bins) {
        if ( ( params.subsample == true ) || params.readcount != null ) {
            match_samps_to_pang(params.samples, core_ch, sub_reads, readcounts, params.project, params.min_cov, params.nr_samps_threshold, params.nr_subsamp)
            pang_samples = match_samps_to_pang.out.pang_samples    
        } else if (params.ref_genome != null) {
            pang_samples = channel.fromPath(params.samples, type: "file", checkIfExists: true)
        } else {
            throw new Exception("No subsampling, no readcount file AND no reference genome? This should never happen.")
        }

        if (params.run_VC) {
            fastq_dir = channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
            variant_calling(fastq_dir, params.subsample, NBPs_ch, pang_samples, contigs_ch, params.ref_genome,
             params.project, params.min_locus_cov, params.min_cov, params.min_breadth, params.min_contig_len, params.block_size )
        }
    }



    /* Currently there is a bug with the event handler:
    https://github.com/nextflow-io/nextflow/issues/5261
    */
    def proj_path = params.project
    workflow.onComplete {
        match_samps_to_pang.out.tot_nr_pangs.view{ nr -> "Out of $nr total pangenomes checked, " }
        match_samps_to_pang.out.passed.view{ nr -> "$nr passed the thresholds for variant calling."}
        println("If this is too few, consider lowering --min_cov and/or --nr_samps_threshold, increasing how many reads are subsampled or using more samples.")
        println("Your results can be found at ${proj_path}\nHave fun!")
    } 
}
