#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 14:55:29 2024

@author: jay
This script is for evaluating the difference in results between using one
binner (metabat2) or 3 binners when running SqueezeMeta
for the pangenome pipeline
python v.3.11.5
fastANI v.1.34
"""
#These variables need to be changed for your purpose:
metabat_res = "/home/jay/data/mock_20240409_metabat"
all_res = "/home/jay/data/mock_20240409_all_binners"
workdir = "/home/jay/data/compare_binners" #output ends up here
ani_threshold = 95

#%%
import glob
import os
import subprocess
import shutil
from Bio import SeqIO
import pandas as pd

os.chdir(workdir)

to_comp = [metabat_res, all_res]
for run in to_comp:
    print(f"Basic stats for {run}:")
    nr_bins = len(glob.glob(f"{run}/bins/*/*.fa"))
    print(f"Total bins in all samples: {nr_bins}")
    nr_motus = len(glob.glob(f"{run}/mOTUs/*_mOTU_*"))
    print(f"Total mOTUs calculated: {nr_motus}")
    pass_cov = len(glob.glob(f"{run}/pangenomes/samples/*.samples"))
    print(f"Pangenomes that have enough samples that passed the cov check: {pass_cov}\n")
    

#Compare both the mOTUs groups, how many overlap etc, and the pangenome ones.
def get_motugroups(res_path):
    groups = []
    motus = glob.glob(f"{res_path}/mOTUs/*_mOTU_*")
    for m in motus:
        group = m.split("c__", 1)[1].split("_mOTU", 1)[0]
        groups.append(group)
    return groups

metabat_groups = get_motugroups(metabat_res)
all_groups = get_motugroups(all_res)

smetabat = set(metabat_groups)
sall = set(all_groups)
motus_union = sall | smetabat
motus_intersect = list(sall & smetabat)

print("All_binners has", len(sall), "taxonomic groups represented among mOTUs")
print("Metabat has", len(smetabat), "taxonomic groups represented among mOTUs")

if motus_intersect != motus_union:
    print("Taxonomic groups present in all binners but not with metabat:", 
          list(sall - smetabat))
    print("Taxonomic groups present using metabat, but not with all binners:", 
          list(smetabat - sall))

#%%
def get_pangsamples(res_path):
    groups = []
    samples = glob.glob(f"{res_path}/pangenomes/samples/*.samples")
    for s in samples:
        group = s.split("c__", 1)[1].split("_mOTU", 1)[0]
        groups.append(group)
    return groups

metabat_pang_groups = get_pangsamples(metabat_res)
all_pang_groups = get_pangsamples(all_res)

smetpang = set(metabat_pang_groups)
sallpang = set(all_pang_groups)
pang_union = list(sallpang | smetpang)
pang_intersect = list(sallpang & smetpang)

#pangenome stats = motus that passed coverage and nr samples checks

if pang_intersect != pang_union:
    print("All_binners has", len(sallpang), "taxonomic groups represented among pangenomes")
    print("Metabat has", len(smetpang), "taxonomic groups represented among pangenomes")
    print("Taxonomic groups present in all binners but not with metabat:", 
          list(sallpang - smetpang))
    print("Taxonomic groups present using metabat, but not with all binners:", 
          list(smetpang - sallpang))
    print("Only the groups present using both methods will be used for ANI computation.")
else:
    print("The same taxonomic groups are represented among pangenomes using both methods.")

#%%run fastANI between pangenomes with same tax-group?
#create files for each tax-group
#I might run into issues if there's a group present in one method and not the other
#so using intersect
print("Creating reference and query files for fastANI.")
for group in pang_intersect:
    for run in to_comp:
        r = run.split("/")[-1]
        samps = glob.glob(f"{run}/mOTUs/pangenomes/*{group}*/*.core.fasta")
        with open(f'{group}_{r}.txt', 'w') as fp:
            fp.write('\n'.join(samps))
    #fastANI --ql [QUERY_LIST] --rl [REFERENCE_LIST] -o [OUTPUT_FILE]
    query = to_comp[0].split("/")[-1]
    ref = to_comp[1].split("/")[-1]
    subprocess.run(["fastANI", "--ql", f"{group}_{query}.txt", "--rl",
                    f"{group}_{ref}.txt", "-o", f"{group}_ANI.out"])
    
#create one tsv of the results
print("Combining fastANI results into one file.")
ani_res = glob.glob("*_ANI.out")
with open('mOTU_identities.txt','wb') as wfd:
    for f in ani_res:
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd)

#%%
#compare sizes of pangenomes
#maybe count number of contigs and the lengths too
with open("pangenome_sizes.tsv", "w") as outfile:
    outfile.write("Pangenome\tLength\tContigs\n")
    for run in to_comp:
        method = run.split("/")[-1]
        nbps = (glob.glob(f"{run}/mOTUs/pangenomes/*_mOTU_*/*.NBPs.fasta"))
        for bin in nbps:
            #pang_name = bin.split("/")[-1]
            pang_name = bin
            all_contigs = list(SeqIO.parse(bin, "fasta"))
            nr_contigs = len(all_contigs)
            pang_length = 0
            for contig in all_contigs:
                pang_length = pang_length + len(all_contigs[0].seq)
            outfile.write(pang_name + "\t" + str(pang_length) + "\t" + str(nr_contigs) + "\n")
            
    
#%%
print("Reading ANI results and sizes into dataframes")
ani_df = pd.read_csv('mOTU_identities.txt', sep = '\t',
                     names=["Metabat", "All_binners", "ANI", "ortho_matches", "Total_fragments"])
ani_df = ani_df[ani_df["ANI"] > ani_threshold ].copy()
size_df = pd.read_csv("pangenome_sizes.tsv", sep="\t")

metabat = metabat_res.split("/")[-1]
binners = all_res.split("/")[-1]
size_df["Method"] = ["Metabat" if metabat in x else "All_binners" for x in size_df["Pangenome"]]
size_df["Pang_name"] = [x.split("/")[-1].split(".")[0] for x in 
                        size_df[size_df["Method"]=="Metabat"]["Pangenome"]] + [
                            x.split("/")[-1].split(".")[0] for x in 
                            size_df[size_df["Method"]=="All_binners"]["Pangenome"]]
                            
metabat_sizes = size_df[size_df["Method"]=="Metabat"]
binners_sizes = size_df[size_df["Method"]=="All_binners"]

#okay so I need to match the right NBPs to the right core. Maybe would have been easier to do earlier?
#for the ANI df, add method column and pang-name column. Then use that to add the checkm and pang size columns.
print("Adding pang_names column per method")
ani_df["Metabat_pang_name"] = [x[-1].split(".")[0] for x in ani_df["Metabat"].str.split("/")]
ani_df["All_binners_pang_name"] = [x[-1].split(".")[0] for x in ani_df["All_binners"].str.split("/")]

def set_sizes(row):
    msize = metabat_sizes[metabat_sizes["Pang_name"]==row["Metabat_pang_name"]]["Length"].item()
    bsize = binners_sizes[binners_sizes["Pang_name"]==row["All_binners_pang_name"]]["Length"].item()
    return msize, bsize

def read_cm1_core(path):
    header = ["Bin Id", "Marker", "Lineage", "# genomes", "# markers",
              "# marker sets", "0", "1", "2", "3", "4", "5+",
              "cM1_Completeness", "cM1_Contamination", "Strain heterogeneity"]
    cm1 = pd.DataFrame(columns=header)
    for f in glob.glob(path+"/checkm_pangenomes/*_cM1_summary.txt"):
        #print(f)
        with open(f) as file:
            for line in file:
                if "c__" in line:
                    #print(line)
                    fields = line.split()
                    cm1.loc[len(cm1)] = fields

    cm1 = cm1[["Bin Id", "cM1_Completeness", "cM1_Contamination"]]
    cm_ret = cm1[cm1["Bin Id"].str.contains("core")].copy()
    cm_ret["Pang_name"] = [x.split(".")[0] for x in cm_ret["Bin Id"]]
    return cm_ret

print("Getting checkm results into dfs")
metabat_core_cm = read_cm1_core(metabat_res)
binners_core_cm = read_cm1_core(all_res)

print("Adding pangenome sizes and checkm results to fastANI results.")
for index, row in ani_df.iterrows():
    msize, bsize = set_sizes(row)
    ani_df.loc[index, "Metabat_bases"] = msize
    ani_df.loc[index, "All_binners_bases"] = bsize
    ani_df.loc[index, ["Metabat_Completeness", "Metabat_Contamination"]] = \
        [float(x) for x in metabat_core_cm[metabat_core_cm["Pang_name"] == \
        row["Metabat_pang_name"]][["cM1_Completeness",
                                   "cM1_Contamination"]].values.tolist()[0]]
    ani_df.loc[index,["All_binners_Completeness", 
                      "All_binners_Contamination"]] = \
        [float(x) for x in binners_core_cm[binners_core_cm["Pang_name"] == \
        row["All_binners_pang_name"]][["cM1_Completeness",
                                      "cM1_Contamination"]].values.tolist()[0]]

#to_add list of pangneome names
#method string of the method used
def create_row_by_method(to_add=[], method="Metabat"):
    ret_df = pd.DataFrame()
    for motu in to_add:
        size = size_df[(size_df["Method"] == method) & (size_df["Pang_name"] == motu)]["Length"].item()
        if method == "Metabat":
            completeness = metabat_core_cm[metabat_core_cm["Pang_name"]==motu]["cM1_Completeness"].item()
            contam = metabat_core_cm[metabat_core_cm["Pang_name"]==motu]["cM1_Contamination"].item()
        elif method == "All_binners":
            completeness = binners_core_cm[binners_core_cm["Pang_name"]==motu]["cM1_Completeness"].item()
            contam = binners_core_cm[binners_core_cm["Pang_name"]==motu]["cM1_Contamination"].item()
        ret_df = ret_df._append({f"{method}_pang_name": motu, f"{method}_bases": size,
                          f"{method}_Completeness": completeness,
                          f"{method}_Contamination": contam}, ignore_index=True)
    return ret_df
    
print("Creating output df")
#make a new df with only necessary columns
output_df = ani_df[["Metabat_pang_name","Metabat_bases","Metabat_Completeness",
                   "Metabat_Contamination","All_binners_pang_name",
                   "All_binners_bases", "All_binners_Completeness",
                   "All_binners_Contamination", "ANI"]].copy()

metabat_motus = [x.split("/")[-1] for x in glob.glob(f"{metabat_res}/mOTUs/*_mOTU_*")]
all_binners_motus = [x.split("/")[-1] for x in glob.glob(f"{all_res}/mOTUs/*_mOTU_*")]
#if a motu is present in metabat_motus but not ani_df, add it to ani_df.

print("Filling output df with the pangenomes without fastANI results")
new_metabat = list(set(metabat_motus) - set(ani_df["Metabat_pang_name"].to_list()))
output_df = output_df._append(create_row_by_method(new_metabat, "Metabat"), ignore_index=True)
new_binners = list(set(all_binners_motus) - set(ani_df["All_binners_pang_name"].to_list()))
output_df = output_df._append(create_row_by_method(new_binners, "All_binners"), ignore_index=True)
output_df = output_df.fillna("NA")

output_df.to_csv("binners_comparison.tsv", sep = '\t', index=False)
                                

