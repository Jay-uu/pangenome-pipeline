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

include { validateParameters ; paramsHelp ; paramsSummaryLog ; fromSamplesheet } from "plugin/nf-validation"


// import subworkflow modules
//include {  } from './modules/subworkflows/'
include { subsample_reads } from './modules/subworkflows/subsample_reads'
include { raw_to_bins } from './modules/subworkflows/raw_to_bins'
include { provided_bins } from './modules/subworkflows/provided_bins'
include { pangenome_assembly } from './modules/subworkflows/pangenome_assembly'
include { match_samps_to_pang } from './modules/subworkflows/match_samps_to_pang'
include { variant_calling } from './modules/subworkflows/variant_calling'

workflow {
    // Print help message, supply typical command line usage for the pipeline
    if (params.help) {
        log.info(paramsHelp("nextflow run main.nf --project <project_name> --samples <tsv.samples> --fastq <path/to/dir> --threads <nr> \n\n  To get more info about a specific parameter write nextflow run main.nf --help <parameter_name>"))
        log.info("You can supply a config file with the parameters using -params-file. For more info check the GitHub page ${workflow.manifest.homePage} ")
    }
    else {
        // Validate input parameters
        validateParameters()

        // Print summary of supplied parameters
        log.info(paramsSummaryLog(workflow))

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
        if (!params.subsample && params.ref_genome == null && params.readcount == null) {
            //If no reference genome: want to always subsample
            //But will allow skipping if readcount file available (entrypoint 1 and 2)
            throw new Exception("When skipping subsampling from raw reads entry or provided bins entry it is assumed that you already have subsampled reads.\nTherefore a tab delimited readcount file with Sample, Nr_fastqs, Total_reads needs to be provided with --readcount <path/to/file>")
        }

        if (params.bins != null && params.ref_genome != null) {
            throw new Exception("You can either provide pre-assembled bins or previously created reference genomes, but not both. Please use either the --bins flag or the --ref_genome flag. ")
        }

        println("Starting. Your results will be published at ${params.project}.")

        //Subsampling
        sam_ch = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
        fastq_ch = Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
        subsample_reads(sam_ch, fastq_ch, params.project, params.subsample, params.nr_subsamp, params.readcount)
        readcounts = subsample_reads.out.readcounts 
        sub_reads = subsample_reads.out.sub_reads
        individual_samp_files = subsample_reads.out.individual_samp_files

        /*
        If the user provided a dir with reference genomes, the pipeline will only run 
        the map_and_detect_variants workflow. This is the third entrypoint.
        */
        if (params.ref_genome != null) {
            Channel.fromPath("${params.ref_genome}", type: "file", checkIfExists: true).multiMap { ch -> core: NBPs: ch }.set { ref_gen }
            
            /*
            When using a reference genome we don't have core and consensus,
            therefore handling the reference as both.
            This means that the whole genome is used both for mapping a subset of the reads,
            and for the variance analysis. Contigs of interest can be provided in a separate file to avoid this.
            */
            core_ch = ref_gen.core
            NBPs_ch = ref_gen.NBPs

            if (params.contigs != null) {
                contigs_ch = Channel.fromPath("${params.contigs}", type: "file", checkIfExists: true)
            }
            else {
                contigs_ch = Channel.empty()
            }
        } else {
            if (params.bins != null) {
            //If bins were provided we don't need to do assembly
            //Optional inputs can be null
            provided_bins(params.samples, params.bins, params.fastq, params.project)
            bins_ch = provided_bins.out.bins
            bintable_ch = provided_bins.out.bintable
            } else {
            /*
    	    Runs assembly and binning.
    	    */
            raw_to_bins(params.fastq, individual_samp_files, params.project, params.binners)
            bins_ch = raw_to_bins.out.bins
            bintable_ch = raw_to_bins.out.bintable
            } 
            /*
    	    This workflow will cluster bins, create pangenomes, and send out core and NBPs fasta files for the pangenomes.
            Done for first and second entrypoints (raw and existing bins).
    	    */
            pangenome_assembly(bins_ch, bintable_ch, params.project, params.MAGcomplete, params.MAGcontam, params.min_mOTU_MAGs,
                               params.min_contig_len, params.taxSort)

            core_ch = pangenome_assembly.out.core_fasta
            NBPs_ch = pangenome_assembly.out.NBPs_fasta
            contigs_ch = pangenome_assembly.out.contigs_tsv
        }
        if (params.subsample == true || params.readcount != null ) {
            match_samps_to_pang(params.samples, core_ch, sub_reads, readcounts, params.project, params.min_cov, params.nr_samps_threshold,
                                params.nr_subsamp)
            pang_samples = match_samps_to_pang.out.pang_samples
        } else if (params.ref_genome != null) {
            pang_samples = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
        } else {
            throw new Exception("No subsampling, no readcount file AND no reference genome? This should never happen.")
        }

        if (params.run_VC == true) {
            variant_calling(params.fastq, params.subsample, NBPs_ch, pang_samples, contigs_ch, params.ref_genome, params.contigs,
                            params.project, params.min_locus_cov, params.min_cov, params.min_breadth, params.min_contig_len, params.block_size)
        }
    }

        /* Currently there is a bug with the event handler:
        https://github.com/nextflow-io/nextflow/issues/5261
        */
        def proj_path = params.project
        workflow.onComplete {
            println("Your results can be found at ${proj_path}\nHave fun!")
        }
        
    }
