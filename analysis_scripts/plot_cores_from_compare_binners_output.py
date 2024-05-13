#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 14:20:13 2024

@author: jay
"""

#using the output from fer_compare_binners.py


from plotnine import (ggplot, aes, geom_boxplot, ggtitle, labs, geom_point)
import pandas as pd
from os import listdir

data_file = "/domus/h1/jay/squeezemeta/pangenome-pipeline/results/binners_eval/compare_pangenomes_output_90ANI.tsv"
#sample_pangs_dir = "/domus/h1/jay/squeezemeta/pangenome-pipeline/results/loclat_202400426_100K_all_binners/pangenomes/samples" 
sample_pangs_dir = "/domus/h1/jay/squeezemeta/pangenome-pipeline/results/loclat_20240503_100K_metabat/pangenomes/samples"
comp_thrsh = 80
#Goals for script
#make boxplots of nr_bins per core mOTU to completeness.

df = pd.read_csv(data_file, sep="\t")


#

#%% make it easier to change method
method = "metabat2"

to_plot = df[[f"mOTU_{method}", f"core_completeness_{method}",f"nr_bins_{method}", \
                  f"NBPs_completeness_{method}"]].copy()
#remove na columns
to_plot = to_plot.dropna()
to_plot["core_comp_group"] = "na"

#sort into two completeness groups
for genome_type in ["core", "NBPs"]:
   to_plot.loc[to_plot[f"{genome_type}_completeness_{method}"] >= comp_thrsh, f"{genome_type}_comp_group"] = "high"
   to_plot.loc[to_plot[f"{genome_type}_completeness_{method}"] < comp_thrsh, f"{genome_type}_comp_group"] = "low"

ggplot(to_plot)+ geom_boxplot(aes(x="core_comp_group", y=f"nr_bins_{method}")) + \
ggtitle(f"Core completeness to nr of bins, all mOTUs\nMethod: {method}") + labs(y="Nr bins", x=f"Completeness, threshold: {comp_thrsh}")

to_plot_no_singlemOTUs = to_plot[to_plot[f"nr_bins_{method}"] > 1].copy()

ggplot(to_plot_no_singlemOTUs)+ geom_boxplot(aes(x='core_comp_group', y=f"nr_bins_{method}")) + \
ggtitle(f"Core completeness to nr of bins\n(excluding mOTUs with only one bin)\nMethod: {method}") + labs(y="Nr bins", x=f"Completeness, threshold: {comp_thrsh}")

#%% I want to make a scatter plot instead. I feel like that will be clearer
ggplot(to_plot_no_singlemOTUs, aes(x=f"core_completeness_{method}", y=f"nr_bins_{method}"))+ \
    geom_point() + ggtitle("Core completeness to nr of bins\n(excluding mOTUs with only one bin)\nMethod: {method}") + labs(y="Nr bins", x="Completeness")
    
 
ggplot(to_plot_no_singlemOTUs)+ geom_boxplot(aes(x='NBPs_comp_group', y=f"nr_bins_{method}")) + \
ggtitle(f"NBPs completeness to nr of bins\n(excluding mOTUs with only one bin)\nMethod: {method}") + labs(y="Nr bins", x=f"Completeness, threshold: {comp_thrsh}")

ggplot(to_plot_no_singlemOTUs, aes(x=f"core_completeness_{method}", y=f"nr_bins_{method}"))+ \
    geom_point() + ggtitle(f"NBPs completeness to nr of bins\n(excluding mOTUs with only one bin)\nMethod: {method}") + labs(y="Nr bins", x="Completeness")
    
ggplot(to_plot) + geom_boxplot(aes(x="NBPs_comp_group", y=f"nr_bins_{method}")) + \
ggtitle(f"NBPs completeness to nr of bins\n(including mOTUs with only one bin)\nMethod: {method}") + labs(y="Nr bins", x=f"Completeness, threshold: {comp_thrsh}")

#%%
samps = [file.split(".")[0] for file in listdir(sample_pangs_dir) if file.endswith(".samples")]

sample_pangs = to_plot[to_plot[f"mOTU_{method}"].isin(samps)].copy()