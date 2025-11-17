#!/usr/bin/env nextflow

//Import modules
include { pang_to_bams } from './../processes/pang_to_bams'
include { downsample_bams_merge } from './../processes/downsample_bams_merge'
include { detect_variants } from './../processes/detect_variants'
include { calc_pang_div } from './../processes/calc_pang_div.nf'

workflow variant_calling {
    take:
    fastq_dir //fastq dir as channel
    subsample //boolean to subsample or not
    NBPs_fasta //channel: path(NBPs.fasta)
    pang_samples //channel: [val(ID), path(ID.samples)] or if subsample == false path(project.samples)
    contigs_tsv //channel: path(contigs.tsv) or null? Does this work?
    ref_genome // null or string
    proj_name //project name string
    min_locus_cov   //minimum locus coverage parameter for pogenom
    min_cov //params.min_cov
    min_breadth //params.min_breadth
    min_contig_len //params.min_contig_len
    block_size //params.block_size

    main:
    NBPs_fasta
        .map { Nfasta -> [Nfasta.getSimpleName(), Nfasta] }
        .set { NBPs_fasta }

    /*
    Using the generated samples files for the pangenome, the raw reads and the pangenome assembly to map reads using SqueezeMeta.
    */
    if (subsample == true) {
        pang_to_bams(NBPs_fasta.combine(pang_samples, by: 0), fastq_dir.first(), proj_name, block_size)
    }
    else {
        //Here the full samples file will be used for each fasta, since no subsampling was done to test which have good coverage
        pang_to_bams(NBPs_fasta.combine(pang_samples), fastq_dir.first(), proj_name, block_size)
    }
    /*
    Checking the breadth and the coverage of bams on the pangenome/ref-genome. Downsampling to even coverage and merging into one bam-file.
    */
    pang_to_bams.out.pang_sqm
        .map { pangsqm -> [pangsqm.getSimpleName(), pangsqm] }
        .set { pang_sqm }

    //If using reference genome, either give empty file if no contigs were provided, or use the contigs_tsv
    if (ref_genome != null) {
        if (contigs_tsv == null) { //Can a channel be compared to null?
            //Optional inputs for processes not yet available: https://github.com/nextflow-io/nextflow/issues/1694
            NO_FILE = file("no_file.text")
            contigs_tsv = Channel.fromPath(NO_FILE, type: "file")
        }
        downsample_bams_merge(pang_sqm.combine(contigs_tsv), proj_name, min_cov, min_breadth, min_contig_len)
    }
    else {
        //If no reference genome, the contigs_tsvs get matched to the right pangenomes based on filenames.
        contigs_tsv
            .map { contsv -> [contsv.getSimpleName(), contsv] }
            .set { contigs_to_downsample }
        downsample_bams_merge(pang_sqm.combine(contigs_to_downsample, by: 0), proj_name, min_cov, min_breadth, min_contig_len)
    }

    downsample_bams_merge.out.not_passed_message.map { msg -> msg.text.strip() }.view()

    /*
    Running freebayes on the merged bam to get a filtered vcf file.
    */
    detect_variants(downsample_bams_merge.out.ref_merged, proj_name)

    //Need to match vcf with right gff
    vcf_gff_ch = detect_variants.out.filt_vcf.combine(pang_to_bams.out.id_gff_genome, by: 0)
    /*
    Run pogenom
    */
    calc_pang_div(vcf_gff_ch, proj_name, min_locus_cov)
    calc_pang_div.out.success_message.map { msg -> msg.text.strip() }.view()
    calc_pang_div.out.pog_path.collectFile(name: "completed_pogs.txt", newLine: true, sort: false, storeDir: proj_name)
}
