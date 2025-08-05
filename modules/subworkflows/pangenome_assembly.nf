#!/usr/bin/env nextflow
//Import process modules
include { parse_taxonomies } from './../processes/parse_taxonomies'
include { bins_to_mOTUs } from './../processes/bins_to_mOTUs'
include { create_mOTU_dirs } from './../processes/create_mOTU_dirs'
include { mOTUs_to_pangenome } from './../processes/mOTUs_to_pangenome'
include { checkm2_pangenomes } from './../processes/checkm2_pangenomes'

workflow pangenome_assembly {
    take:
    bins //channel: path(bin.fa)
    bintable //channel: path(sample.bintable)
    proj_name //project name string
    MAGcomplete // Minimum MAG completeness threshold
    MAGcontam //Maximum MAG contamination threshold
    min_mOTU // min_mOTU_MAGs parameter
    min_contig_len //min_contig_len parameter
    tax_sort //params.taxSort
    main:
    /*
    Before running mOTUlizer, the checkM and GTDB-Tk outputs (bintables) need to be parsed.
    All bintables and all bins from different samples need to be collected so the taxonomy_parser 
    process can run once with all data.
    */
    bintable.collect().set { all_bintables }
    bins.collect().multiMap { ch -> to_tax_parser: to_mOTU_dirs: ch }.set { all_bins }

    parse_taxonomies(all_bins.to_tax_parser, all_bintables, tax_sort, MAGcomplete, MAGcontam)

    /*
    	Clustering of bins, if they've been presorted to lower taxonomic ranks this can spawn parallell processes
    	*/
    bins_to_mOTUs(parse_taxonomies.out.tax_bin_dirs.flatten(), proj_name, MAGcomplete, MAGcontam)

    /*
    	Creating dirs for the mOTUs by sorting based on the mOTUlizer output,
    	so each mOTU directory has the correct bins.
    	*/
    create_mOTU_dirs(bins_to_mOTUs.out.mOTUs_tuple, all_bins.to_mOTU_dirs, MAGcontam, min_mOTU)

    /*
    	Running SuperPang, creating pangenomes. Transpose makes it so that each mOTU from the same grouping within
    	the taxonomy selection will be sent individually to the process together with the matching bintable.
    	*/
    mOTUs_to_pangenome(create_mOTU_dirs.out.transpose(), proj_name, min_contig_len)

    /*
    	The checkm2 process will need to be updated when checkm2 is in the SQM env.
    	*/
    checkm2_pangenomes(mOTUs_to_pangenome.out.pangenome_dir, proj_name)

    emit:
    core_fasta = mOTUs_to_pangenome.out.core_fasta //channel: path(pangenomes/${mOTU_dir}/*.core.fasta)
    NBPs_fasta = mOTUs_to_pangenome.out.NBPs_fasta //channel: path(pangenomes/${mOTU_dir}/*.NBPs.fasta)
    contigs_tsv = mOTUs_to_pangenome.out.contigs_tsv //channel: path(pangenomes/${mOTU_dir}_core_contigs.tsv)
}
