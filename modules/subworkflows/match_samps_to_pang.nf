#!/usr/bin/env nextflow

//Import process modules
include { index_coreref } from './../processes/index_coreref'
include { map_subset } from './../processes/map_subset'
include { cov_to_pang_samples } from './../processes/cov_to_pang_samples'

workflow match_samps_to_pang {
    take:
    sample_path // samples file path param
    core_fasta //channel: path(core.fasta)
    sub_reads //channel: path()
    readcounts //channel: path()
    proj_name //project path param
    minimum_coverage //min_cov param
    nr_samps_threshold // nr_samps_threshold param
    nr_subsamp //nr_subsamp param

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
    sample_file = channel.fromPath(sample_path, type: "file", checkIfExists: true)

    cov_to_pang_samples(map_subset.out.coverage.collect(), map_subset.out.stats.collect(), sample_file.first(), readcounts.collect(), proj_name, minimum_coverage, nr_samps_threshold, nr_subsamp)
    cov_to_pang_samples.out.not_passed_message.map { msg -> msg.text.strip() }.view()
    //This file only gets created if not enough samples, meaning the text only gets printed if pipeline stops here.

    cov_to_pang_samples.out.pang_samples.flatten().map { psam -> [psam.getSimpleName(), psam] }.set { pang_samples }

    emit:
    pang_samples = pang_samples //channel: [val(ID), path(ID.samples)]
}
