#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 13:48:28 2024

@author: jay
"""
"""
This script is for investigating the core qualities.
"""
import os
import pandas as pd
import numpy as np
import glob
#After giving Superpang checkm res:
#checkm_sum = "/crex/proj/fume/nobackup/private/jay/test_pipeline/loclat_20240110/checkm_pangenomes/combined_cm1_cm2_summary.tsv"
#Without checkm res:
#checkm_sum = "/crex/proj/fume/nobackup/private/jay/test_pipeline/loclat_20240110/checkm_pangenomes/combined_cm1_cm2_summary_20240110.tsv"
#workdir = "/crex/proj/fume/nobackup/private/jay/test_pipeline/loclat_20240110"
#motu_dir = "/crex/proj/fume/nobackup/private/jay/test_pipeline/loclat_20240110/mOTUs"

#projname = "loclat_20240202_all_binners_w_superpangcheckm"
#projdir = "/crex/proj/fume/nobackup/private/jay/test_pipeline/" + projname

projname = "loclat_202400426_100K_all_binners"
projdir = "/domus/h1/jay/squeezemeta/pangenome-pipeline/results/" + projname

checkm_sum = projdir + "/checkm_pangenomes/" + projname + "_combined_cm1_cm2_summary.tsv"
workdir = projdir
motu_dir = projdir + "/mOTUs"
os.chdir(f'{workdir}')

#%%
for line in open(checkm_sum):
    if line.startswith('Bin Id'):
        continue
    if 'Actinomycetia' in line:
        print(line)
    Bin_Id, cM1_Completeness, cM1_Contamination, cM2_Completeness, cM2_Contamination = line.strip().split("\t")

assert False

#%%
checkm = pd.read_csv(checkm_sum, sep="\t")
#filter res to only keep mOTUs of interest
#dont want the ones with only a single genome (singlemOTU in name)
#also need to remove the assembly I think?
#index_names = checkm[ (checkm['cM2_Completeness'] <= 80) | (checkm["Bin Id"].str.contains("singlemOTU") == True)
#                     | (checkm["Bin Id"].str.contains(".assembly") == True) ].index 
index_names = checkm[ (checkm["Bin Id"].str.contains("singlemOTU") == True) | (checkm["Bin Id"].str.contains(".assembly") == True) ].index 
checkm.drop(index_names, inplace = True)

#keep only pangenomes where at least either the core or .NBPs has more than 80% completeness
#get all "bin Ids" then remove endings, then go through each and if no row with that has good completeness, add to a list both NBPs and .core
nbps = checkm[ (checkm["Bin Id"].str.endswith("NBPs") == True) ]
pangnames = nbps["Bin Id"].to_list()
passed = []
failed = []
for pang in pangnames:
    print(f"Checking {pang}")
    if ((len(checkm[checkm["Bin Id"] == pang])) + 
         (len(checkm[checkm["Bin Id"] == pang + ".core"])) == 2):
        print (f"{pang} has both NBPs and core")
        if (checkm[checkm["Bin Id"] == pang]["cM2_Completeness"].item() >= 80 
            | (checkm[checkm["Bin Id"] == pang + ".core"]["cM2_Completeness"].item() > 80)):
            print(f"{pang} NBPs or core has good completeness. Adding to list.")
            passed.append(pang)
            passed.append(pang+".core")
        else:
            print(f"{pang} does not pass completeness criteria.")
            failed.append(pang)
            failed.append(pang+".core")
    else:
        print(f"{pang} is missing a core")
        failed.append(pang)
        failed.append(pang+".core")

#Actually I should do this on uppmax
#The pangenomes we want to investigate:
print("Copying pangenomes to new df")
df = checkm[checkm["Bin Id"].isin(passed)].copy()
#need to add number of bins per pangenome also
print("Adding mOTU names")    
df["mOTU_name"] = df['Bin Id'].str.split('.', n=1, expand=True)[0]
uniq_motu_names = df["mOTU_name"].unique()
print("Adding nr bins")
df["Nr_bins"] = ""
for m in uniq_motu_names:
    bins = len(os.listdir(motu_dir+"/"+m))
    df.loc[df['mOTU_name'] == m, 'Nr_bins'] = bins
print("Adding the new core column")    
#might need to add a core column?
df["Type"] = df['Bin Id'].str.split('.', expand=True)[2]
df["Type"] = df["Type"].fillna("NBPs")
#and a qual_group column
df['Completeness_group'] = pd.cut(df['cM2_Completeness'], bins=range(0, 121, 60))



#%% Need to get all bins completeness, and then add those for a total mOTU bins completeness
#this code doesn't seem to work anymore because sqm no longer gives the bintable
sqm_dirs = projdir + "/pangenomes/sqm/"

bintable_list = []
for file in glob.glob(sqm_dirs + "*/results/18.*bintable"):
    print(file)
    tmp = pd.read_csv(file, sep="\t", skiprows=1)
    bintable_list.append(tmp[["Bin ID", "Completeness"]].copy())
    
bintables = pd.concat(bintable_list, axis=0, ignore_index=True)
bintables["Bin ID"] = bintables["Bin ID"].astype(str) + ".fa"

#I guess now I check which bins are in which mOTU/pangenome and add the total completeness?
#loop-di-loop
df["Tot_Comp"] = 0
for m in uniq_motu_names:
    bins = os.listdir(motu_dir+"/"+m)
    tot_comp = 0
    for b in bins:
        tot_comp = tot_comp + bintables[bintables["Bin ID"]==b]["Completeness"].item() #CONTINUE HERE
    df.loc[df['mOTU_name'] == m, 'Tot_Comp'] = tot_comp   
  
#%%
print("Time to plot")
bp = df.boxplot("Nr_bins", by=["Type", "Completeness_group"], ylabel = "Nr bins")

only_core = df[df["Type"]=="core"].copy()
only_core.boxplot("Nr_bins", by="Completeness_group", ylabel = "Nr bins", xlabel="Core Completeness")
only_nbps = df[df["Type"]=="NBPs"].copy()
only_nbps.boxplot("Nr_bins", by=["Type", "Completeness_group"], ylabel = "Nr bins NBPs" )

only_core.groupby("Completeness_group").count()

only_core.boxplot("Tot_Comp", by="Completeness_group", ylabel = "Total bin completeness", xlabel="Core Completeness")
only_nbps = df[df["Type"]=="NBPs"].copy()

#%%
#columns motuname, nbps completeness, core completeness, , nr bins
