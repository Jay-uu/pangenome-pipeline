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

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run main.nf --project <project_name> --samples <tsv.samples> --fastq <path/to/dir> --threads <nr> \n\n  To get more info about a specific parameter write nextflow run main.nf --help <parameter_name>")
   log.info """\
   You can supply a config file with the parameters using -c. For more info see the usr_template.config or
   ${params.manifest.homePage}
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

if (workflow.resume == false) {
	//Workflow was not resumed, checking project dir
	Path projDir = new File(params.project).toPath()
	if (projDir.exists() == true) {
		throw new Exception("Project directory $params.project already exists, choose a new name or use the -resume flag. WARNING: Note that if you resume the wrong job, this might overwrite previous results.")
	}
}

// import modules
//maybe later I will move the workflows into separate files and only import the necessary modules for that workflow
include { format_samples } from './modules/format_samples'
include { fastq_to_bins } from './modules/fastq_to_bins'
include { subsample_fastqs } from './modules/subsample_fastqs'
include { parse_taxonomies } from './modules/parse_taxonomies'
include { bins_to_mOTUs } from './modules/bins_to_mOTUs'
include { create_mOTU_dirs } from './modules/create_mOTU_dirs'
include { mOTUs_to_pangenome } from './modules/mOTUs_to_pangenome'
include { checkm_pangenomes } from './modules/checkm_pangenomes'
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
//include {  } from './modules/'


workflow raw_to_bins {
    take:
    	samples_files    
    main:
        fastq_dir = Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
    	fastq_to_bins(samples_files, fastq_dir.first())
    emit:
    	bins = fastq_to_bins.out.bins
    	bintable = fastq_to_bins.out.bintable
    	
}

workflow provided_bins {
    main:
        sample_file = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
        //there should be a check that there's fastas in the dir too, maybe in the workflow or the process?
        bins_dir = Channel.fromPath(params.bins, type: "dir", checkIfExists: true)
        fastq_dir = Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
        classify_bins(sample_file, bins_dir, fastq_dir.first())
    emit:
    	bins = classify_bins.out.bins
    	bintable = classify_bins.out.bintable

}


workflow pangenome_assembly {
    take:
        bins
    	bintable
    
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
    	FOR NOW, RUNNING CHECKM HERE BECAUSE THE INPUT NEEDS A DIR. WILL PROBABLY CHANGE THAT SO CHECKM CAN RUN WITH THE REF GENOME START TOO
    	BUT WE HAVEN'T YET DECIDED WHICH CHECKM VERSION TO USE, SO I WON'T EDIT YET CAUSE THIS IS EASIER.
    	*/
    	mOTUs_to_pangenome.out.pangenome_dir.multiMap { dir -> to_checkm: to_checkm2: dir }.set { pang_dirs }
    	checkm_pangenomes(pang_dirs.to_checkm)
        checkm2_pangenomes(pang_dirs.to_checkm2)
    
    emit:
    	core_fasta = mOTUs_to_pangenome.out.core_fasta
    	NBPs_fasta = mOTUs_to_pangenome.out.NBPs_fasta
    
}

workflow variant_calling {
    take:
    single_samples
    core_fasta
    NBPs_fasta
    
    main:
    //Going to mutliple processes
    Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
    		.multiMap { it -> to_subsamp: to_pang_to_bams: it }.set { fastq_dir }
    
    //Concatenating fastqs and subsampling for later mapping for each singles sample
    subsample_fastqs(single_samples, fastq_dir.to_subsamp.first())
    
    core_fasta.multiMap { it -> to_coreref: to_downsample: it }.set { core_fasta }
    /*
    Index genomes for read mapping
    */
    index_coreref(core_fasta.to_coreref)

    /*
    map subset reads to pangenome and get coverage information
    */
    map_subset(index_coreref.out.fasta_index_id.combine(subsample_fastqs.out.sub_reads))

    /*
    Using the coverage from the mapping, decides which reads "belong" to which pangenome and creates new .samples files
    */
    sample_file = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
    cov_to_pang_samples(map_subset.out.coverage.collect(),sample_file.first(), subsample_fastqs.out.readcount.collect())

    /*
    Create keys to match right samples file to right NBPs fasta (from same pangenome) as input to pang_to_bams.
    */
    cov_to_pang_samples.out.pang_samples
		.flatten()
		.map { [it.getSimpleName(), it] }
		.set { pang_samples }

    NBPs_fasta
		.map { [it.getSimpleName(), it] }
		.set { NBPs_fasta }

    /*
    Using the generated samples files for the pangenome, the raw reads and the pangenome assembly to map reads using SqueezeMeta.
    */
    pang_to_bams(pang_samples.combine(NBPs_fasta, by: 0), fastq_dir.to_pang_to_bams.first())

    /*
    Checking the breadth and the coverage of bams on the pangenome/ref-genome. Downsampling to even coverage and merging into one bam-file.
    */
    //ADD THE CORE FASTA HERE, AND MAYBE CHANGE OUTPUT HANDLING TO MATCH CORRECT PANGENOMES
    pang_to_bams.out.pang_sqm
		.map { [it.getSimpleName(), it] }
		.set { pang_sqm }

    core_fasta.to_downsample
		.map { [it.getSimpleName(), it] }
		.set { core_fasta.to_downsample }
    
    pang_sqm.combine(core_fasta.to_downsample, by: 0)).view()

    downsample_bams_merge(pang_sqm.combine(core_fasta.to_downsample, by: 0))

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
    
}


workflow {
    //Maybe I should have a check for incompatible parameters? For example if both ref_genomes and bins were provided?
    
    /*The fastq_dir is needed for:
	- Formating the individual sample files
	- Assembly (if starting from raw reads)
	- The variant calling workflow
    */
    fastq_chan = Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)
     
    //File with which fastq files belong to which samples. Tab delimited with sample-name, fastq file name and pair.
    sam_chan = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
    
    /*Runs the process that creates individual samples files and creates two output channels:
    Would like to find a way to not unneccesarily create the to_assembly channel, if assembly is skipped.
	- For assembly
	- for variant calling
    */
    format_samples(sam_chan, fastq_chan)
    format_samples.out.flatten().multiMap { it -> to_assembly: to_variants: it }.set { single_samps }
    
    /*
    If the user provided a dir with reference genomes, the pipeline will only run 
    the map_and_detect_variants workflow.
    */
    if ( params.ref_genomes != null ) {
        Channel.fromPath( "${params.ref_genomes}/*.{fasta, fa}" , type: "file", checkIfExists: true ).multiMap { it -> core: NBPs: it }.set { ref_gens }
        /*
        When using a reference genome we don't have core and consensus,
        therefore handling the reference as both.
        This means that the whole genome is used both for mapping a subset of the reads,
        and for the variance analysis.
        */
        core_ch = ref_gens.core
        NBPs_ch = ref_gens.NBPs   
    }
    /*
    If no reference genome directory was provided, pangenomes will be constructed.
    */
    else {
        //If bins were provided we don't need to do assembly, and only need the singles .samples files for the variant calling workflow.
        if ( params.bins != null ) {
        provided_bins() //need to test if the output looks as I want it.
        bins_ch = provided_bins.out.bins
        bintable_ch = provided_bins.out.bintable
        }
        
        else {
    	/*
    	Runs assembly and binning.
    	*/
    	raw_to_bins(single_samps.to_assembly)
    	bins_ch = raw_to_bins.out.bins
        bintable_ch = raw_to_bins.out.bintable
    	}
    	/*
    	This workflow will cluster bins, create pangenomes, and send out core and NBPs fasta files for the pangenomes.
    	*/
    	pangenome_assembly(bins_ch, bintable_ch)
    	
    	core_ch = pangenome_assembly.out.core_fasta
    	NBPs_ch = pangenome_assembly.out.NBPs_fasta

    }
    
    /*
    This workflow maps a subset of the reads to each pangenome to estimate which would pass the coverage checks.
    The samples that pass the initial coverage checks are then used as samples for the pangenome/ref genome for
    variance analysis, by mapping the reads, estimating coverage and breadth, downsampling etc.
    Creates VCF files. Will add so it also creates pogenom results.
    */
    variant_calling(single_samps.to_variants, core_ch, NBPs_ch)
    	
    //It should be possible to add a message for when the pipeline finishes.
    
}

