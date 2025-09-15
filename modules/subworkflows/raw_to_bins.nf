#!/usr/bin/env nextflow
include { fastq_to_bins } from './../processes/fastq_to_bins'
include { summarize_bintables } from './../helper_functions/summarize_bintables'

workflow raw_to_bins {
    take:
    fq_path //fastq param
    individual_samp_files //channel of individual samples files
    proj_name //project name param
    binners //params.binners
    contig_len //params.contig_len
    main:
    //The dir with all the fastqs
    Channel.fromPath(fq_path, type: "dir", checkIfExists: true).multiMap { ch -> to_format: to_assembly: to_subsamp: ch }.set { fastq_ch }
    //File with which fastq files belong to which samples. Tab delimited with sample-name, fastq file name and pair.
    //sam_ch = Channel.fromPath(samp_file, type: "file", checkIfExists: true)

    /*Runs the process that creates individual samples files
    format_samples(sam_ch, fastq_ch.to_format)
    format_samples.out.flatten().multiMap { ch -> to_subsamp: to_assembly: ch }.set { samples_files }
    */

    //Binning
    fastq_to_bins(individual_samp_files, fastq_ch.to_assembly.first(), proj_name, binners, contig_len)

    //Summarizing bintables into one file and only printing certain columns
    fastq_to_bins.out.bintable.multiMap { ch -> to_emit: to_summarize: ch }.set { ch_bintables }
    summarize_bintables(ch_bintables.to_summarize, proj_name)

    emit:
    bins = fastq_to_bins.out.bins //channel: path(ID/results/bins/*.fa)
    bintable = ch_bintables.to_emit //channel: path(ID/results/18.ID.bintable)
}
