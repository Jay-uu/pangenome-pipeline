#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 09:22:37 2024

@author: jay
"""

#Filter the header of a bam file to only include contigs that are
#also in the contigs.tsv
import os

workdir = "/home/jay/squeezemeta/pangenome-pipeline/results/test_bam_filtering"
os.chdir(f'{workdir}')

sam = "only_core.sam"
out = "filtered_header.sam"
contigs = "contigs.tsv"

contig_list = []
with open(contigs) as file:
    for line in file:
        contig_list.append(line.split("\t")[0])

with open(out, "w") as outfile:    
    with open(sam) as file:
        for line in file:
            if line.startswith("@SQ"):
                test = line
                contig = line.split("\t")[1].split(":")[1]
                if contig not in contig_list:
                    continue
            outfile.write(line)
        