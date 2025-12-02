#!/usr/bin/env nextflow
//include { fastq_to_bins } from './../processes/fastq_to_bins'
include { assemble_metagenome } from './../processes/assemble_metagenome'
include { map_samples } from './../processes/map_samples'
include { binning } from './../processes/binning'
include { dastool } from './../processes/dastool'
include { checkbins } from './../processes/checkbins'
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
    channel.fromPath(fq_path, type: "dir", checkIfExists: true).multiMap { ch -> to_format: to_assembly: to_subsamp: ch }.set { fastq_ch }

    //Assembly and Binning
    //fastq_to_bins(individual_samp_files, fastq_ch.to_assembly.first(), proj_name, binners, contig_len)
    assemble_metagenome(individual_samp_files, fastq_ch.to_assembly.first(), binners, contig_len)
    map_samples(assemble_metagenome.out.sqm_dir)
    binning(map_samples.out.sqm_dir)
    dastool(binning.out.sqm_dir, proj_name)
    
    checkbins(dastool.out.sqm_dir, proj_name)
    checkbins.out.bintable.multiMap { ch -> to_emit: to_summarize: ch }.set { ch_bintables }

    summarize_bintables(ch_bintables.to_summarize, proj_name)

    emit:
    //bins = fastq_to_bins.out.bins //channel: path(ID/results/bins/*.fa)
    bins = dastool.out.bins
    bintable = ch_bintables.to_emit //channel: path(ID/results/18.ID.bintable)
}
