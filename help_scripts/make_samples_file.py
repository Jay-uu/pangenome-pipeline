#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: jay
"""

#Create a samples file with format sample_name, readfile, pair1 or pair2
#This script only works on files with format similar to SampleLake_fwd.fastq.gz
import os
from glob import glob

OUTPUT_DIR = "<path/to/place>"
SAMPLES_NAME = "<sample_name>.samples"
FASTQ_DIR = "<path/to/read/files>"
REV = "rev" #how you mark in the name of the read file which direction of reads it is
FWD = "fwd"
DIR_DIVIDER = "_" #direction divider. What's between sample name and direction in name of file
FQ_FORMAT = ".fastq.gz"

read_files = glob(FASTQ_DIR + "/*.fastq.gz") + glob(FASTQ_DIR + "/*.fq.gz")
read_files.sort()

with open(OUTPUT_DIR + "/" + SAMPLES_NAME, "w") as outfile:
    for rf in read_files:
        file_name = rf.split("/")[-1]
        name,pair = file_name.split(FQ_FORMAT)[0].split(DIR_DIVIDER)
        if pair == REV:
            pair = "pair2"
        elif pair == FWD: 
            pair = "pair1"
        else:
            raise ValueError("Could not identify direction of file {file_name}.")
        outfile.write(f"{name}\t{file_name}\t{pair}\n")
        

