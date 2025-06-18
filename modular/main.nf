#!/usr/bin/env nextflow
/*
========================================================================================
Pipeline for pangenome intra-diversity analysis
========================================================================================
Modular attempt
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from "plugin/nf-validation" 


// import modules
//maybe later I will move the workflows into separate files and only import the necessary modules for that workflow
include { format_samples } from './modules/format_samples'
include { fastq_to_bins } from './modules/fastq_to_bins'
include { subsample_fastqs } from './modules/subsample_fastqs'
include { parse_taxonomies } from './modules/parse_taxonomies'
include { bins_to_mOTUs } from './modules/bins_to_mOTUs'
include { create_mOTU_dirs } from './modules/create_mOTU_dirs'
include { mOTUs_to_pangenome } from './modules/mOTUs_to_pangenome'
//include { checkm_pangenomes } from './modules/checkm_pangenomes'
include { checkm2_pangenomes } from './modules/checkm2_pangenomes'
include { index_pangenomes } from './modules/index_pangenomes'
include { index_coreref } from './modules/index_coreref'
include { map_subset } from './modules/map_subset'
include { cov_to_pang_samples } from './modules/cov_to_pang_samples'
include { pang_to_bams } from './modules/pang_to_bams'
include { downsample_bams_merge } from './modules/downsample_bams_merge'
include { detect_variants } from './modules/detect_variants'
include { classify_bins } from './modules/classify_bins.nf'
include { calc_pang_div } from './modules/calc_pang_div.nf'
include { tuplify_samp_fastqs } from './modules/tuplify_samp_fastqs.nf'

def summarize_bintables(bintable_ch) {
	def bintable_header = Channel.value( "Bin ID	Completeness	Contamination	Tax GTDB-Tk" )
	def bintable_rows = bintable_ch.collectFile(keepHeader: true, skip: 2)
    			.splitCsv(header: true, skip: 1, sep: "\t")
    			.map { row -> "${row.'Bin ID'}	${row.Completeness}	${row.Contamination}	${row.'Tax GTDB-Tk'}" }
    	bintable_header.concat(bintable_rows).collectFile(name: "summarized_bintable.tsv", newLine: true, sort: false, storeDir: "${params.project}/bins")
}

def concat_readcounts(readcounts_ch) {
	readcounts_ch.collectFile(keepHeader: true, name: "original_readcounts.tsv", newLine: false, sort: false, storeDir: "${params.project}/subsamples")
}

def concat_subsamp_samples(sample_ch) {
	sample_ch.collectFile(name: "${params.project}.subsampled.samples", newLine: true, storeDir: "${params.project}/subsamples")
}

workflow subsample_reads {
    take:
    samples_files
    fastq_ch
    main:
    	//Subsample fastq files
    	subsample_fastqs(samples_files, fastq_ch.first())
        
        //Mapping channel to be able to concatenate the readcounts and publish them, as well as send to next WF.
        subsample_fastqs.out.readcount.multiMap { chan -> to_emit: to_concat: chan }.set { ch_readcounts }
        
        //Publish concatenated files
        concat_readcounts(ch_readcounts.to_concat)
    	concat_subsamp_samples(subsample_fastqs.out.sample_file)

    emit:
    	sub_reads = subsample_fastqs.out.sub_reads //channel: [val("ID"), path("sub_ID.fq.gz")]
    	readcounts = ch_readcounts.to_emit //channel: path(ID_readcount.txt)
}


workflow raw_to_bins {  
    main:
    	//The dir with all the fastqs
    	Channel.fromPath(params.fastq, type: "dir", checkIfExists: true).multiMap { it -> to_format: to_assembly: to_subsamp: it }.set { fastq_ch }
     
    	//File with which fastq files belong to which samples. Tab delimited with sample-name, fastq file name and pair.
    	sam_ch = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
    
    	//Runs the process that creates individual samples files
    	format_samples(sam_ch, fastq_ch.to_format)
    	format_samples.out.flatten().multiMap { it -> to_subsamp: to_assembly: it }.set { samples_files }
    
        //Binning
    	fastq_to_bins(samples_files.to_assembly, fastq_ch.to_assembly.first())
    	
    	//Summarizing bintables into one file and only printing certain columns
    	fastq_to_bins.out.bintable.multiMap { chan -> to_emit: to_summarize: chan }.set { ch_bintables }
    	summarize_bintables(ch_bintables.to_summarize)

    	//Concatenating fastqs and subsampling for later mapping for each singles sample
    	subsample_reads(samples_files.to_subsamp, fastq_ch.to_subsamp)


    emit:
    	bins = fastq_to_bins.out.bins //channel: path(ID/results/bins/*.fa)
    	bintable = ch_bintables.to_emit //channel: path(ID/results/18.ID.bintable)
    	sub_reads = subsample_reads.out.sub_reads //channel: [val("ID"), path("sub_ID.fq.gz")]
    	readcounts = subsample_reads.out.readcounts //channel: path(ID_readcount.txt)
    	
}

workflow provided_bins {
    main:
        sample_file = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
        //there should be a check that there's fastas in the dir too, maybe in the workflow or the process?
        bins_dir = Channel.fromPath(params.bins, type: "dir", checkIfExists: true)
        fastq_dir = Channel.fromPath(params.fastq, type: "dir", checkIfExists: true) //should be subsampled fastqs provided by user
	//if taxonomy and completeness already provided, don't need to run this. Might add something for it in the future.
        classify_bins(sample_file, bins_dir, fastq_dir.first())
        classify_bins.out.bintable.multiMap { chan -> to_emit: to_summarize: chan }.set { ch_bintables }
        summarize_bintables(ch_bintables.to_summarize)
        
        sam_ch = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
        fastq_ch = Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
        fastq_ch.multiMap { chan -> to_format: to_sub: chan }.set { fastq_ch }
        
        format_samples(sam_ch, fastq_ch.to_format)
        samples_files = format_samples.out.flatten()
        
        if ( params.subsample == true ) {
            subsample_reads(samples_files, fastq_ch.to_sub)           
    	    readcounts = subsample_reads.out.readcounts
    	    sub_reads = subsample_reads.out.sub_reads

        }
        else {
            //tuplify_samp_fastqs output has the same format as subsample_reads.out.sub_reads tuple(val("${sample.baseName}"), path("sub_*.fq.gz"), emit: sub_reads)
            tuplify_samp_fastqs(samples_files, fastq_ch.to_sub.first())
            sub_reads = tuplify_samp_fastqs.out.reads
            readcounts = Channel.fromPath( "${params.readcount}", type: "file", checkIfExists: true )
        
        }
    emit:
    	bins = classify_bins.out.bins //channel: path(ID/results/bins/*.fa)
    	bintable = ch_bintables.to_emit //channel: path(ID/results/18.ID.bintable)
    	sub_reads = sub_reads //channel: [val("ID"), path("sub_ID.fq.gz")]
    	readcounts = readcounts //channel: path(ID_readcount.txt)

}


workflow pangenome_assembly {
    take:
        bins		//channel: path(bin.fa)
    	bintable	//channel: path(sample.bintable)
    
    main:    
    	/*
    	Before running mOTUlizer, the checkM and GTDB-Tk outputs (bintables) need to be parsed.
    	All bintables and all bins from different samples need to be collected so the taxonomy_parser 
    	process can run once with all data.
    	*/

    	bintable.collect().set { all_bintables }
    	bins.collect().multiMap { it -> to_tax_parser: to_mOTU_dirs: it }.set { all_bins }
    
    	parse_taxonomies(all_bins.to_tax_parser, all_bintables) 
    
    	/*
    	Clustering of bins, if they've been presorted to lower taxonomic ranks this can spawn parallell processes
    	*/
    	bins_to_mOTUs(parse_taxonomies.out.tax_bin_dirs.flatten())

    	/*
    	Creating dirs for the mOTUs by sorting based on the mOTUlizer output,
    	so each mOTU directory has the correct bins.
    	*/
    	create_mOTU_dirs(bins_to_mOTUs.out.mOTUs_tuple, all_bins.to_mOTU_dirs)

    	/*
    	Running SuperPang, creating pangenomes. Transpose makes it so that each mOTU from the same grouping within
    	the taxonomy selection will be sent individually to the process together with the matching bintable.
    	*/
    	mOTUs_to_pangenome(create_mOTU_dirs.out.transpose())
    	
    	/*
    	The checkm2 process will need to be updated when checkm2 is in the SQM env.
    	*/
        checkm2_pangenomes(mOTUs_to_pangenome.out.pangenome_dir)
    
    emit:
    	core_fasta = mOTUs_to_pangenome.out.core_fasta //channel: path(pangenomes/${mOTU_dir}/*.core.fasta)
    	NBPs_fasta = mOTUs_to_pangenome.out.NBPs_fasta //channel: path(pangenomes/${mOTU_dir}/*.NBPs.fasta)
    	contigs_tsv = mOTUs_to_pangenome.out.contigs_tsv //channel: path(pangenomes/${mOTU_dir}_core_contigs.tsv)
    
}

workflow match_samps_to_pang {
    take:
    core_fasta		//channel: path(core.fasta)
    sub_reads		//channel: path()
    readcounts		//channel: path()
    main:    
    /*
    Index genomes for read mapping
    */
    index_coreref(core_fasta)
    /*
    map subset reads to pangenome and get coverage information
    */
    map_subset(index_coreref.out.fasta_index_id.combine(sub_reads))
    /*
    Using the coverage from the mapping, decides which reads "belong" to which pangenome and creates new .samples files
    */
    sample_file = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
    cov_to_pang_samples(map_subset.out.coverage.collect(), sample_file.first(), readcounts.collect())
    cov_to_pang_samples.out.not_passed_message.map { it.text.strip() }.view() //This file only gets created if not enough samples, meaning the text only gets printed if pipeline stops here.
    cov_to_pang_samples.out.pang_samples.flatten().map { [it.getSimpleName(), it] }.set { pang_samples } //if no pang_samples pipeline stops here, but results from cov_to_pang_samples still get published
    
    emit:
    pang_samples = pang_samples //channel: [val(ID), path(ID.samples)]

}

workflow variant_calling {
    take:
    NBPs_fasta		//channel: path(NBPs.fasta)
    pang_samples	//channel: [val(ID), path(ID.samples)] or if subsample == false path(project.samples)
    contigs_tsv	//channel: path(contigs.tsv)
    main:    
    //Going to mutliple processes
    fastq_dir = Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
    
    NBPs_fasta
		.map { [it.getSimpleName(), it] }
		.set { NBPs_fasta }
    
    /*
    Using the generated samples files for the pangenome, the raw reads and the pangenome assembly to map reads using SqueezeMeta.
    */
    if ( params.subsample == true ) {
    	pang_to_bams(NBPs_fasta.combine(pang_samples, by: 0), fastq_dir.first())
	}
    else {
    	//Here the full samples file will be used for each fasta, since no subsampling was done to test which have good coverage
    	pang_to_bams(NBPs_fasta.combine(pang_samples), fastq_dir.first())
    }
    /*
    Checking the breadth and the coverage of bams on the pangenome/ref-genome. Downsampling to even coverage and merging into one bam-file.
    */
    pang_to_bams.out.pang_sqm
		.map { [it.getSimpleName(), it] }
		.set { pang_sqm }

    //If using reference genome, either give empty file if no contigs were provided, or use the contigs_tsv
    if ( params.ref_genome != null ) {
        if ( params.contigs == null ) {
    	//Optional inputs for processes not yet available: https://github.com/nextflow-io/nextflow/issues/1694
        NO_FILE = file("no_file.text")
        contigs_tsv = Channel.fromPath(NO_FILE, type: "file")
        }
        downsample_bams_merge(pang_sqm.combine(contigs_tsv))
    }
    else {
	contigs_tsv
                .map { [it.getSimpleName(), it] }
                .set { contigs_to_downsample }
        downsample_bams_merge(pang_sqm.combine(contigs_to_downsample, by: 0))
    }

    downsample_bams_merge.out.not_passed_message.map { it.text.strip() }.view()
    
    /*
    Running freebayes on the merged bam to get a filtered vcf file.
    */
    detect_variants(downsample_bams_merge.out.ref_merged)
    
    //Need to match vcf with right gff
    vcf_gff_ch = detect_variants.out.filt_vcf.combine(pang_to_bams.out.id_gff_genome, by: 0)
    /*
    Run pogenom
    */
    calc_pang_div(vcf_gff_ch)
    calc_pang_div.out.success_message.map { it.text.strip() }.view()
    
}


workflow {
    // Print help message, supply typical command line usage for the pipeline
    if (params.help) {
        log.info paramsHelp("nextflow run main.nf --project <project_name> --samples <tsv.samples> --fastq <path/to/dir> --threads <nr> \n\n  To get more info about a specific parameter write nextflow run main.nf --help <parameter_name>")
        log.info """\
        You can supply a config file with the parameters using -c. For more info see the usr_template.config or
        ${workflow.manifest.homePage}
        """

        exit 0
    }

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

//Check project parameter
def badChars = ["^","(",")","+", " ", "|"]
if ( params.project.findAll { a -> badChars.any { a.contains(it) } } ) {
	throw new Exception("Invalid project name. Special characters and whitespaces not allowed.")
}

//Gives a warning if project already exists
if (workflow.resume == false) {
	//Workflow was not resumed, checking project dir
	Path projDir = new File(params.project).toPath()
	if (projDir.exists() == true) {
		throw new Exception("Project directory $params.project already exists, choose a new name or use the -resume flag. WARNING: Note that if you resume the wrong job, this might overwrite previous results.")
	}
}

println "The subsample parameter is set to ${params.subsample}"
if (params.subsample == false) {
	//Currently we want to allow skipping subsampling for any entry
	//if (params.bins == null) {
	//    throw new Exception("Skipping subsampling is only allowed for the bin entry. Please provide a directory with --bins <path/to/dir> or set --subsample <true>.")
	//}
	if (params.bins != null && params.readcount == null) {
	    throw new Exception("When skipping subsampling and using already constructed bins a tab delimited readcount file with Sample, Nr_fastqs, Total_reads needs to be provided with --readcount <path/to/file>")
	    //file format Sample  Nr_fastqs       Total_reads
	    //MAYBE UPDATE TO NOT REQUIRE NR_FASTQS BECAUSE IT'S OBSOLETE NOW
	}
    }
    
if (params.bins != null && params.ref_genome != null) {
	throw new Exception("You can either provide pre-assembled bins or previously created reference genomes, but not both. Please use either the --bins flag or the --ref_genome flag. ")
}

    println "Starting. Your results will be published at ${params.project}."
    
    /*
    If the user provided a dir with reference genomes, the pipeline will only run 
    the map_and_detect_variants workflow.
    */
    if ( params.ref_genome != null ) {
        Channel.fromPath( "${params.ref_genome}" , type: "file", checkIfExists: true ).multiMap { it -> core: NBPs: it }.set { ref_gen }
        /*
        When using a reference genome we don't have core and consensus,
        therefore handling the reference as both.
        This means that the whole genome is used both for mapping a subset of the reads,
        and for the variance analysis.
        */
        core_ch = ref_gen.core
        NBPs_ch = ref_gen.NBPs
        
        if ( params.contigs != null ) {
            contigs_ch = Channel.fromPath( "${params.contigs}" , type: "file", checkIfExists: true ) //IF THIS WASNT PROVIDED, USE WHOLE FASTA.
        }
        else {	
            contigs_ch = Channel.empty() //Channel.fromPath("no_file.txt", type: "file")
        }
        
        if ( params.subsample == true ) {
            sam_ch = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
            fastq_ch = Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
            
            format_samples(sam_ch, fastq_ch)
            samples_files = format_samples.out.flatten()
            subsample_reads(samples_files, fastq_ch)
                        
    	    readcounts = subsample_reads.out.readcounts
    	    sub_reads = subsample_reads.out.sub_reads
    	    
    	    match_samps_to_pang(core_ch, sub_reads, readcounts) 
    	}
    }
    /*
    If no reference genome directory was provided, pangenomes will be constructed.
    */
    else {
        //If bins were provided we don't need to do assembly
        if ( params.bins != null ) {
        provided_bins()
        bins_ch = provided_bins.out.bins
        bintable_ch = provided_bins.out.bintable
        sub_reads = provided_bins.out.sub_reads
        readcounts = provided_bins.out.readcounts
        }
        
        else {
    	    /*
    	    Runs assembly and binning.
    	    */
    	    raw_to_bins()
    	    bins_ch = raw_to_bins.out.bins
            bintable_ch = raw_to_bins.out.bintable
            sub_reads = raw_to_bins.out.sub_reads
            readcounts = raw_to_bins.out.readcounts
    	    }
    	/*
    	This workflow will cluster bins, create pangenomes, and send out core and NBPs fasta files for the pangenomes.
    	*/
    	pangenome_assembly(bins_ch, bintable_ch)
    	
    	core_ch = pangenome_assembly.out.core_fasta
    	NBPs_ch = pangenome_assembly.out.NBPs_fasta
    	contigs_ch = pangenome_assembly.out.contigs_tsv
	
	core_ch.multiMap { it -> to_map: to_variants: it }.set { core_ch }
    	match_samps_to_pang(core_ch.to_map, sub_reads, readcounts)
    }
    if ( params.run_VC == true ) {
    	if ( params.subsample == true ) {
    		//Use only the samples with estimated good coverage from the subsampled reads.
    		variant_calling(NBPs_ch, match_samps_to_pang.out.pang_samples, contigs_ch)
    		}
    	else if (params.force_variant_calling == true) {
    		//Using the full read files/all samples to map reads.
    		variant_calling(NBPs_ch, Channel.fromPath(params.samples, type: "file", checkIfExists: true), contigs_ch)
    		}
	else {
		throw new Exception("It seems you're trying to run the variant calling workflow without subsampling. If you're sure about doing this, make sure to set --force_variant_calling to true.")
		}
    	}

    //It should be possible to add a message for when the pipeline finishes.
    workflow.onComplete {
        println "Your results can be found at ${params.project}\nHave fun!"
    }
}


