#!/usr/bin/env nextflow
/*
========================================================================================
Pipeline for pangenome intra-diversity analysis
========================================================================================
Modular attempt
----------------------------------------------------------------------------------------
*/

include { validateParameters ; paramsHelp ; paramsSummaryLog ; fromSamplesheet } from "plugin/nf-validation"


// import modules
//maybe later I will move the workflows into separate files and only import the necessary modules for that workflow
include { format_samples } from './modules/format_samples'
include { fastq_to_bins } from './modules/fastq_to_bins'
include { subsample_fastqs } from './modules/subsample_fastqs'
include { parse_taxonomies } from './modules/parse_taxonomies'
include { bins_to_mOTUs } from './modules/bins_to_mOTUs'
include { create_mOTU_dirs } from './modules/create_mOTU_dirs'
include { mOTUs_to_pangenome } from './modules/mOTUs_to_pangenome'
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

def summarize_bintables(bintable_ch, proj_name) {
    def bintable_header = Channel.value("Bin ID	Completeness	Contamination	Tax GTDB-Tk")
    def bintable_rows = bintable_ch
        .collectFile(keepHeader: true, skip: 2)
        .splitCsv(header: true, skip: 1, sep: "\t")
        .map { row -> "${row.'Bin ID'}	${row.Completeness}	${row.Contamination}	${row.'Tax GTDB-Tk'}" }
    bintable_header.concat(bintable_rows).collectFile(name: "summarized_bintable.tsv", newLine: true, sort: false, storeDir: "${proj_name}/bins")
}

def concat_readcounts(readcounts_ch, proj_name) {
    readcounts_ch.collectFile(keepHeader: true, name: "original_readcounts.tsv", newLine: false, sort: false, storeDir: "${proj_name}/subsamples")
}

def concat_subsamp_samples(sample_ch, proj_name) {
    sample_ch.collectFile(name: "${proj_name}.subsampled.samples", newLine: true, storeDir: "${proj_name}/subsamples")
}

workflow subsample_reads {
    take:
    samples_files //samples path
    fastq_ch // fastq dir path
    proj_name //project name string

    main:
    //Subsample fastq files
    subsample_fastqs(samples_files, fastq_ch.first())

    //Mapping channel to be able to concatenate the readcounts and publish them, as well as send to next WF.
    subsample_fastqs.out.readcount.multiMap { ch -> to_emit: to_concat: ch }.set { ch_readcounts }

    //Publish concatenated files
    concat_readcounts(ch_readcounts.to_concat, proj_name)
    concat_subsamp_samples(subsample_fastqs.out.sample_file, proj_name)

    emit:
    sub_reads = subsample_fastqs.out.sub_reads //channel: [val("ID"), path("sub_ID.fq.gz")]
    readcounts = ch_readcounts.to_emit //channel: path(ID_readcount.txt)
}


workflow raw_to_bins {
    take:
    fq_path //fastq param
    samp_file //samples file param
    proj_name //project name param
    main:
    //The dir with all the fastqs
    Channel.fromPath(fq_path, type: "dir", checkIfExists: true).multiMap { ch -> to_format: to_assembly: to_subsamp: ch }.set { fastq_ch }

    //File with which fastq files belong to which samples. Tab delimited with sample-name, fastq file name and pair.
    sam_ch = Channel.fromPath(samp_file, type: "file", checkIfExists: true)

    //Runs the process that creates individual samples files
    format_samples(sam_ch, fastq_ch.to_format)
    format_samples.out.flatten().multiMap { ch -> to_subsamp: to_assembly: ch }.set { samples_files }

    //Binning
    fastq_to_bins(samples_files.to_assembly, fastq_ch.to_assembly.first())

    //Summarizing bintables into one file and only printing certain columns
    fastq_to_bins.out.bintable.multiMap { ch -> to_emit: to_summarize: ch }.set { ch_bintables }
    summarize_bintables(ch_bintables.to_summarize, proj_name)

    //Concatenating fastqs and subsampling for later mapping for each singles sample
    subsample_reads(samples_files.to_subsamp, fastq_ch.to_subsamp, proj_name)

    emit:
    bins = fastq_to_bins.out.bins //channel: path(ID/results/bins/*.fa)
    bintable = ch_bintables.to_emit //channel: path(ID/results/18.ID.bintable)
    sub_reads = subsample_reads.out.sub_reads //channel: [val("ID"), path("sub_ID.fq.gz")]
    readcounts = subsample_reads.out.readcounts //channel: path(ID_readcount.txt)
}

workflow provided_bins {
    take:
    sample_path // samples file path
    bins_path // path to dir with bins in fasta format
    fastq_path // fastq dir path
    to_subsamp // boolean to subsample or not
    proj_name // if subsamp need proj_name
    readcount // if not subsamp, need readcount file path. Can be null if subsamp is true.
    main:
    sam_ch = Channel.fromPath(sample_path, type: "file", checkIfExists: true)
    sam_ch.multiMap { ch -> to_classify: to_format: ch }.set { sam_ch }

    //there should be a check that there's fastas in the dir too, maybe in the workflow or the process?
    bins_dir = Channel.fromPath(bins_path, type: "dir", checkIfExists: true)

    fastq_ch = Channel.fromPath(fastq_path, type: "dir", checkIfExists: true)
    fastq_ch.multiMap { ch -> to_classify: to_format: to_sub: ch }.set { fastq_ch }
    
    //if taxonomy and completeness already provided, don't need to run this. Might add something for it in the future.
    classify_bins(sam_ch.to_classify, bins_dir, fastq_ch.to_classify.first())
    classify_bins.out.bintable.multiMap { ch -> to_emit: to_summarize: ch }.set { ch_bintables }
    summarize_bintables(ch_bintables.to_summarize, proj_name)

    format_samples(sam_ch.to_format, fastq_ch.to_format)
    samples_files = format_samples.out.flatten()

    if (to_subsamp == true) {
        subsample_reads(samples_files, fastq_ch.to_sub, proj_name)
        readcounts = subsample_reads.out.readcounts
        sub_reads = subsample_reads.out.sub_reads
    }
    else {
        //If skipping subsampling for bin entry, it's assumed that the fastq files are already subsampled
        //that's why a file with original readcounts is needed.
        //tuplify_samp_fastqs output has the same format as subsample_reads.out.sub_reads tuple(val("${sample.baseName}"), path("sub_*.fq.gz"), emit: sub_reads)
        tuplify_samp_fastqs(samples_files, fastq_ch.to_sub.first())
        sub_reads = tuplify_samp_fastqs.out.reads
        readcounts = Channel.fromPath("${readcount}", type: "file", checkIfExists: true)
    }

    emit:
    bins = classify_bins.out.bins //channel: path(ID/results/bins/*.fa)
    bintable = ch_bintables.to_emit //channel: path(ID/results/18.ID.bintable)
    sub_reads = sub_reads //channel: [val("ID"), path("sub_ID.fq.gz")]
    readcounts = readcounts //channel: path(ID_readcount.txt)
}


workflow pangenome_assembly {
    take:
    bins //channel: path(bin.fa)
    bintable //channel: path(sample.bintable)
    proj_name //project name string
    MAGcomplete // Minimum MAG completeness threshold
    MAGcontam //Maximum MAG contamination threshold

    main:
    /*
    Before running mOTUlizer, the checkM and GTDB-Tk outputs (bintables) need to be parsed.
    All bintables and all bins from different samples need to be collected so the taxonomy_parser 
    process can run once with all data.
    */
    proj_ch = Channel.value(proj_name)
    MAGcomplete_ch = Channel.value(MAGcomplete)
    MAGcontam_ch = Channel.value(MAGcontam)

    bintable.collect().set { all_bintables }
    bins.collect().multiMap { ch -> to_tax_parser: to_mOTU_dirs: ch }.set { all_bins }

    parse_taxonomies(all_bins.to_tax_parser, all_bintables)

    /*
    	Clustering of bins, if they've been presorted to lower taxonomic ranks this can spawn parallell processes
    	*/
    bins_to_mOTUs(parse_taxonomies.out.tax_bin_dirs.flatten(), proj_ch, MAGcomplete_ch, MAGcontam_ch)

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
    checkm2_pangenomes(mOTUs_to_pangenome.out.pangenome_dir, proj_ch)

    emit:
    core_fasta = mOTUs_to_pangenome.out.core_fasta //channel: path(pangenomes/${mOTU_dir}/*.core.fasta)
    NBPs_fasta = mOTUs_to_pangenome.out.NBPs_fasta //channel: path(pangenomes/${mOTU_dir}/*.NBPs.fasta)
    contigs_tsv = mOTUs_to_pangenome.out.contigs_tsv //channel: path(pangenomes/${mOTU_dir}_core_contigs.tsv)
}

workflow match_samps_to_pang {
    take:
    sample_path // samples file path
    core_fasta //channel: path(core.fasta)
    sub_reads //channel: path()
    readcounts //channel: path()

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
    sample_file = Channel.fromPath(sample_path, type: "file", checkIfExists: true)
    cov_to_pang_samples(map_subset.out.coverage.collect(), sample_file.first(), readcounts.collect())
    cov_to_pang_samples.out.not_passed_message.map { msg -> msg.text.strip() }.view()
    //This file only gets created if not enough samples, meaning the text only gets printed if pipeline stops here.
    cov_to_pang_samples.out.pang_samples.flatten().map { psam -> [psam.getSimpleName(), psam] }.set { pang_samples }

    emit:
    pang_samples = pang_samples //channel: [val(ID), path(ID.samples)]
}

workflow variant_calling {
    take:
    fastq_path //fastq dir path
    subsample //boolean to subsample or not
    NBPs_fasta //channel: path(NBPs.fasta)
    pang_samples //channel: [val(ID), path(ID.samples)] or if subsample == false path(project.samples)
    contigs_tsv //channel: path(contigs.tsv)
    ref_genome
    contigs_path

    main:
    //Going to mutliple processes
    fastq_dir = Channel.fromPath(fastq_path, type: "dir", checkIfExists: true)

    NBPs_fasta
        .map { Nfasta -> [Nfasta.getSimpleName(), Nfasta] }
        .set { NBPs_fasta }

    /*
    Using the generated samples files for the pangenome, the raw reads and the pangenome assembly to map reads using SqueezeMeta.
    */
    if (subsample == true) {
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
        .map { pangsqm -> [pangsqm.getSimpleName(), pangsqm] }
        .set { pang_sqm }

    //If using reference genome, either give empty file if no contigs were provided, or use the contigs_tsv
    if (ref_genome != null) {
        if (contigs_path == null) {
            //Optional inputs for processes not yet available: https://github.com/nextflow-io/nextflow/issues/1694
            NO_FILE = file("no_file.text")
            contigs_tsv = Channel.fromPath(NO_FILE, type: "file")
        }
        downsample_bams_merge(pang_sqm.combine(contigs_tsv))
    }
    else {
        contigs_tsv
            .map { contsv -> [contsv.getSimpleName(), contsv] }
            .set { contigs_to_downsample }
        downsample_bams_merge(pang_sqm.combine(contigs_to_downsample, by: 0))
    }

    downsample_bams_merge.out.not_passed_message.map { msg -> msg.text.strip() }.view()

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
    calc_pang_div.out.success_message.map { msg -> msg.text.strip() }.view()
}


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
        if (workflow.resume == false) {
            //Workflow was not resumed, checking project dir
            def Path projDir = new File(params.project).toPath()
            if (projDir.exists() == true) {
                throw new Exception("Project directory ${params.project} already exists, choose a new name or use the -resume flag. WARNING: Note that if you resume the wrong job, this might overwrite previous results.")
            }
        }

        println("The subsample parameter is set to ${params.subsample}")
        if (params.subsample == false) {
            //Currently we want to allow skipping subsampling for any entry
            //if (params.bins == null) {
            //    throw new Exception("Skipping subsampling is only allowed for the bin entry. Please provide a directory with --bins <path/to/dir> or set --subsample <true>.")
            //}
            if (params.bins != null && params.readcount == null) {
                throw new Exception("When skipping subsampling and using already constructed bins a tab delimited readcount file with Sample, Nr_fastqs, Total_reads needs to be provided with --readcount <path/to/file>")
            }
        }

        if (params.bins != null && params.ref_genome != null) {
            throw new Exception("You can either provide pre-assembled bins or previously created reference genomes, but not both. Please use either the --bins flag or the --ref_genome flag. ")
        }

        println("Starting. Your results will be published at ${params.project}.")

        /*
        If the user provided a dir with reference genomes, the pipeline will only run 
        the map_and_detect_variants workflow.
        */
        if (params.ref_genome != null) {
            Channel.fromPath("${params.ref_genome}", type: "file", checkIfExists: true).multiMap { ch -> core: NBPs: ch }.set { ref_gen }
            
        /*
        When using a reference genome we don't have core and consensus,
        therefore handling the reference as both.
        This means that the whole genome is used both for mapping a subset of the reads,
        and for the variance analysis.
        */
            core_ch = ref_gen.core
            NBPs_ch = ref_gen.NBPs

            if (params.contigs != null) {
                contigs_ch = Channel.fromPath("${params.contigs}", type: "file", checkIfExists: true)
            }
            else {
                contigs_ch = Channel.empty()
            }

            if (params.subsample == true) {
                sam_ch = Channel.fromPath(params.samples, type: "file", checkIfExists: true)
                fastq_ch = Channel.fromPath(params.fastq, type: "dir", checkIfExists: true)

                format_samples(sam_ch, fastq_ch)
                samples_files = format_samples.out.flatten()
                subsample_reads(samples_files, fastq_ch, params.project)

                readcounts = subsample_reads.out.readcounts
                sub_reads = subsample_reads.out.sub_reads

                match_samps_to_pang(params.samples, core_ch, sub_reads, readcounts)
            }
        }
        else {
            //If bins were provided we don't need to do assembly
            if (params.bins != null) {
                //How to deal with readcount optional input?
                provided_bins(params.samples, params.bins, params.fastq, params.subsample, params.project, params.readcount)
                bins_ch = provided_bins.out.bins
                bintable_ch = provided_bins.out.bintable
                sub_reads = provided_bins.out.sub_reads
                readcounts = provided_bins.out.readcounts
            }
            else {
                /*
    	        Runs assembly and binning.
    	        */
                raw_to_bins(params.fastq, params.samples, params.project)
                bins_ch = raw_to_bins.out.bins
                bintable_ch = raw_to_bins.out.bintable
                sub_reads = raw_to_bins.out.sub_reads
                readcounts = raw_to_bins.out.readcounts
            }
            /*
    	This workflow will cluster bins, create pangenomes, and send out core and NBPs fasta files for the pangenomes.
    	*/
            pangenome_assembly(bins_ch, bintable_ch, params.project, params.MAGcomplete, params.MAGcontam)

            core_ch = pangenome_assembly.out.core_fasta
            NBPs_ch = pangenome_assembly.out.NBPs_fasta
            contigs_ch = pangenome_assembly.out.contigs_tsv

            core_ch.multiMap { ch -> to_map: to_variants: ch }.set { core_ch }
            match_samps_to_pang(params.samples, core_ch.to_map, sub_reads, readcounts)
        }
        if (params.run_VC == true) {
            if (params.subsample == true) {
                //Use only the samples with estimated good coverage from the subsampled reads.
                variant_calling(params.fastq, params.subsample, NBPs_ch, match_samps_to_pang.out.pang_samples, contigs_ch, params.ref_genome, params.contigs)
            }
            else if (params.force_variant_calling == true) {
                //Using the full read files/all samples to map reads.
                variant_calling(params.fastq, params.subsample, NBPs_ch, Channel.fromPath(params.samples, type: "file", checkIfExists: true), contigs_ch, params.ref_genome, params.contigs)
            }
            else {
                throw new Exception("It seems you're trying to run the variant calling workflow without subsampling. If you're sure about doing this, make sure to set --force_variant_calling to true.")
            }
        }

        //It should be possible to add a message for when the pipeline finishes.
        workflow.onComplete {
            println("Your results can be found at ${params.project}\nHave fun!")
        }
    }
}
