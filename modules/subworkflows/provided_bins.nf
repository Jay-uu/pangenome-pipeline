#!/usr/bin/env nextflow
//Import process modules
include { classify_bins } from './../processes/classify_bins'
include { checkbins} from './../processes/checkbins'
include { format_samples } from './../processes/format_samples'
include { summarize_bintables } from './../helper_functions/summarize_bintables'

workflow provided_bins {
    take:
    sample_path // samples file path
    bins_path // path to dir with bins in fasta format
    fastq_path // fastq dir path
    proj_name // needed for publishdir
    main:
    sam_ch = channel.fromPath(sample_path, type: "file", checkIfExists: true)
    sam_ch.multiMap { ch -> to_classify: to_format: ch }.set { sam_ch }

    //there should be a check that there's fastas in the dir too, maybe in the workflow or the process?
    bins_dir = channel.fromPath(bins_path, type: "dir", checkIfExists: true)

    fastq_ch = channel.fromPath(fastq_path, type: "dir", checkIfExists: true)
    fastq_ch.multiMap { ch -> to_classify: to_format: to_sub: ch }.set { fastq_ch }
    
    //if taxonomy and completeness already provided, don't need to run this. Might add something for it in the future.
    classify_bins(sam_ch.to_classify, bins_dir, fastq_ch.to_classify.first(), proj_name)
    checkbins(classify_bins.out.sqm_dir, proj_name)
    checkbins.out.bintable.multiMap { ch -> to_emit: to_summarize: ch }.set { ch_bintables }

    summarize_bintables(ch_bintables.to_summarize, proj_name)

    emit:
    bins = classify_bins.out.bins //channel: path(ID/results/bins/*.fa)
    bintable = ch_bintables.to_emit //channel: path(ID/results/18.ID.bintable)
}
