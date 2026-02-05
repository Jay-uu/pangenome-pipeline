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
    pang_cov_stats = map_subset.out.pang_cov_stats.groupTuple(by: 0) //group by pang_ID, size would be nr of samples. Can I get that? It would be nr of readcounts files?
    /*
    Using the coverage from the mapping, decides which reads "belong" to which pangenome and creates new .samples files
    */
    sample_file = channel.fromPath(sample_path, type: "file", checkIfExists: true)
    //
    cov_to_pang_samples(pang_cov_stats, sample_file.first(), readcounts.collect(), proj_name, minimum_coverage, nr_samps_threshold, nr_subsamp)
    //cov_to_pang_samples.out.not_passed_message.map { msg -> msg.text.strip() }.view()
    // how to manage if null?
    //not_passed = cov_to_pang_samples.out.not_passed.count() //.collectFile(name: 'NOT_PASSED.txt', newLine: false)
    cov_to_pang_samples.out.passed.multiMap{ nr -> to_tot: to_emit: nr }.set{ passed } //.collectFile(name: 'PASSED.txt', newLine: false)
    tot_nr_pangs = passed.to_tot.mix(cov_to_pang_samples.out.not_passed).count()
    //flatten might not be neccessary anymore
    cov_to_pang_samples.out.pang_samples.flatten().map { psam -> [psam.getSimpleName(), psam] }.set { pang_samples }

    emit:
    pang_samples = pang_samples //channel: [val(ID), path(ID.samples)]
    passed = passed.to_emit.count()
    tot_nr_pangs = tot_nr_pangs

}
