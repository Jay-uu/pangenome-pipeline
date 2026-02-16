#!/usr/bin/env nextflow

//Import modules
include { concat_readcounts } from './../helper_functions/concat_readcounts'
include { concat_subsamp_samples } from './../helper_functions/concat_subsamp_samples'
include { subsample_fastqs } from './../processes/subsample_fastqs'
include { format_samples } from './../processes/format_samples'
include { tuplify_samp_fastqs } from './../processes/tuplify_samp_fastqs'

workflow subsample_reads {
    take:
    sample_file //samples path as channel
    fastq_dir // fastq dir path
    proj_name //project name string
    to_subsamp //boolean params.subsample
    nr_subsamp //params.nr_subsamp
    readcount_file //null or file with Sample, Nr_fastqs, Total_reads

    main:
    individual_samp_files = null //set default values for outputs
    sub_reads = null
    readcounts = null

    //fastqs are needed for formatting and subsampling/tuplifying
    fastq_dir.multiMap { ch -> to_sub_or_tup: to_format: ch }.set { fastq_ch }
    //Creating indivudal sample files from the one given
    format_samples(sample_file, fastq_ch.to_format)
    format_samples.out.flatten().multiMap { ch -> to_emit: to_sub_or_tup: ch }.set { samples_ch }
    individual_samp_files = samples_ch.to_emit

    if (to_subsamp) {
        //Subsample fastq files
        subsample_fastqs(samples_ch.to_sub_or_tup, fastq_ch.to_sub_or_tup.first(), proj_name, nr_subsamp)
        sub_reads = subsample_fastqs.out.sub_reads

        //Mapping channel to be able to concatenate the readcounts and publish them, as well as send to next WF.
        subsample_fastqs.out.readcount.multiMap { ch -> to_emit: to_concat: ch }.set { ch_readcounts }
        readcounts = ch_readcounts.to_emit

        //Publish concatenated files
        concat_readcounts(ch_readcounts.to_concat, proj_name)
        concat_subsamp_samples(subsample_fastqs.out.sample_file, proj_name)

    } else if ( readcount_file != null ) {
        //If skipping subsampling for bin entry, it's assumed that the fastq files are already subsampled
        //that's why a file with original readcounts is needed.
        //tuplify_samp_fastqs output has the same format as subsample_reads.out.sub_reads tuple(val("${sample.baseName}"), path("sub_*.fq.gz"), emit: sub_reads)
        tuplify_samp_fastqs(samples_ch.to_sub_or_tup, fastq_ch.to_sub_or_tup.first())
        sub_reads = tuplify_samp_fastqs.out.reads //channel: [val("ID"), path("sub_ID.fq.gz")]
        readcounts = Channel.fromPath(readcount_file, type: "file", checkIfExists: true)        
    }

    emit:
    sub_reads = sub_reads //channel: [val("ID"), path("sub_ID.fq.gz")]
    readcounts = readcounts //channel: path(ID_readcount.txt)
    individual_samp_files = individual_samp_files
}