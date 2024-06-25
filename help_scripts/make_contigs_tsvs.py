#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Provided a dir with fasta files, will create tsv for each with matching names
#that have the contig names and lengths.
from glob import glob
from Bio import SeqIO

genomes_path = "input" #path to dir with genomes you want to analyse
outdir_path = "output" #where you want the resulting .contigs files

min_contig_len = 1000
genomes = glob(genomes_path+"/*.fa") + glob(genomes_path+"/*.fasta")

for fasta in genomes:
    fname = fasta.split("/")[-1].split(".")[0]
    with open(fasta, "rt") as pg:
        all_contigs = list(SeqIO.parse(pg, "fasta"))
        with open(f"{outdir_path}/{fname}.contigs.tsv", "w") as outfile:
            for i,s in enumerate(all_contigs):
                c_len = len(s.seq)
                if c_len >= min_contig_len:
                    outfile.write(f"{s.id}\t{c_len}\n")
