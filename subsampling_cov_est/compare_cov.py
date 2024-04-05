#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 14:08:55 2024

@author: jay
Excuse the very messy script... it's an ongoing process that will likely not
be re-adjusted to be pretty!
"""
#to compare coverages, but I also need to compare the mOTUs from the different runs.

import os
import pandas as pd
import glob

import numpy as np
import matplotlib.pyplot as plt
os.chdir('/home/jay/compare_cov')
file1M = '/home/jay/compare_cov/1000000_cov.tsv'
file100K = '/home/jay/compare_cov/100000_cov.tsv'
file10K = '/home/jay/compare_cov/10000_cov.tsv'

cov = pd.DataFrame()
for file in glob.glob("*_cov.tsv"):
    df = pd.read_csv(file, sep="\t")
    df["nr_subsampling"] = file.split("_", 1)[0]
    cov = pd.concat([cov, df], ignore_index=True)
    
#if this is right none should be NA
cov.isnull().any().any()
    

df_100K = pd.read_csv(file100K, sep="\t")
df_1M = pd.read_csv(file1M, sep="\t")
df_10K = pd.read_csv(file10K, sep="\t")


#%%
df_1M["Pangenome"].equals(df_100K["Pangenome"]) #true
df_1M.columns.equals(df_100K.columns) #true
df_1M.columns.equals(df_10K.columns)
df_1M["Pangenome"].equals(df_10K["Pangenome"])

oneM = df_1M.drop("Pangenome", axis = 1).to_numpy().flatten()
hundredK = df_100K.drop("Pangenome", axis = 1).to_numpy().flatten()
tenK = df_10K.drop("Pangenome", axis = 1).to_numpy().flatten()

plt.figure(dpi=1200)
plt.scatter(oneM, hundredK, s = 1, linewidths=0)
plt.title("Coverage") 
plt.xlabel("1M sub-reads") 
plt.ylabel("100K sub-reads")
plt.axvline(x=20, color='r', linestyle='-')
plt.axhline(y=20, color='r', linestyle='-')
plt.show()

#%%
plt.figure(dpi=1200)
plt.scatter(oneM, tenK, s = 1, linewidths=0)
plt.title("Coverage") 
plt.xlabel("1M sub-reads") 
plt.ylabel("10K sub-reads")
plt.axvline(x=20, color='r', linestyle='-')
plt.axhline(y=20, color='r', linestyle='-')
plt.show()

#%% log transform
a = oneM
b = hundredK
c = tenK
minval = np.min(a[np.nonzero(a)])
min_add = 0.001
maxval = np.max(a[np.nonzero(a)])

a = np.log(a + min_add)
b = np.log(b + min_add)
c = np.log(c + min_add)

#%%
plt.figure(dpi=1200)
plt.scatter(a, b, s = 1, linewidths=0)
plt.title("Coverage, log transformed") 
plt.xlabel("1M sub-reads") 
plt.ylabel("100K sub-reads")
plt.axvline(x=np.log(20), color='r', linestyle='-')
plt.axhline(y=np.log(20), color='r', linestyle='-')
plt.show()

#%%
plt.figure(dpi=1200)
plt.scatter(a, c, s = 1, linewidths=0)
plt.title("Coverage, natural log transformed") 
plt.xlabel("1M sub-reads") 
plt.ylabel("10K sub-reads")
plt.axvline(x=np.log(20), color='r', linestyle='-')
plt.axhline(y=np.log(20), color='r', linestyle='-')
plt.show()

#%% check how many samples are above 20 cov for each mOTU, and if it is the same for both methods.
#maybe I can just make a new column with how many of the others passed the cov criteria?
cols = df_1M.columns[1:13]
df_1M["passed"] = 0
df_100K["passed"] = 0
df_10K["passed"] = 0
cov_threshold = 10

#do it like this for now and see if I can use list comprehension and/or apply later
for motu in df_1M["Pangenome"].unique():
    passed1M = 0
    passed100K = 0
    passed10K = 0
    for col in cols:
        if df_1M[df_1M["Pangenome"]==motu][col].item() >= cov_threshold:
            passed1M = passed1M + 1
        if df_100K[df_100K["Pangenome"]==motu][col].item() >= cov_threshold:
            passed100K = passed100K + 1
        if df_10K[df_10K["Pangenome"]==motu][col].item() >= cov_threshold:
            passed10K = passed10K +1
    df_1M.loc[df_1M[df_1M["Pangenome"]==motu].index, "passed"] = passed1M
    df_100K.loc[df_100K[df_100K["Pangenome"]==motu].index, "passed"] = passed100K
    df_10K.loc[df_10K[df_10K["Pangenome"]==motu].index, "passed"] = passed10K
    
df_1M["Pangenome"].equals(df_100K["Pangenome"]) #true
df_1M.columns.equals(df_100K.columns) #true
df_1M["passed"].equals(df_100K["passed"])
(df_1M["passed"] == df_100K["passed"]).all()
#some of them differ, find out how many, which ones etc.
print("The treshold is: ", cov_threshold)
print("Total passed 1M: ", df_1M["passed"].sum())
print("Total passed 100K: ", df_100K["passed"].sum())
print("Total passed 10K: ", df_10K["passed"].sum())


#to 1M df add columns for mean cov and std

#binarize coverage matrices, 1 if cov >= 20 otherwise 0. Then calculate jaccard dissimilarity.
#Probably use apply right?
cols = df_1M.columns[1:13]
df = df_1M
df["mean_cov"] = df.iloc[:, 1:13].mean(axis=1)
df["std"] = df.iloc[:, 1:13].std(axis=1)

bdf_1M = df.copy()
bdf_1M[cols] = df[cols].gt(cov_threshold).astype(int)
bdf_100K = df_100K.copy()
bdf_100K[cols] = df_100K[cols].gt(cov_threshold).astype(int)
bdf_10K = df_10K.copy()
bdf_10K[cols] = df_10K[cols].gt(cov_threshold).astype(int)

#add dissimilarity scores to bdf_1M
bdf_1M["100K_similarity"] = 0.0
bdf_1M["10K_similarity"] = 0.0
for index, row in bdf_1M.iterrows():
    s1M = row[cols].to_numpy().flatten()
    s100K = bdf_100K.iloc[index,:][cols].to_numpy().flatten()
    s10K = bdf_10K.iloc[index,:][cols].to_numpy().flatten()
    if (sum(s1M==1)+sum(s100K==1)==0):
        bdf_1M.loc[index,"100K_similarity"] = 1
    else:
        insect = sum((s1M==s100K)[s1M==1])
        union = len([1 for i in range(0, len(s1M)) if s1M[i]==1 or s100K[i]==1])
        sim = (float(insect)/float(union))
        bdf_1M.loc[index,"100K_similarity"] = sim
    if (sum(s1M==1)+sum(s10K==1)==0):
        bdf_1M.loc[index,"10K_similarity"] = 1
    else:
        insect = sum((s1M==s10K)[s1M==1])
        union = len([1 for i in range(0, len(s1M)) if s1M[i]==1 or s10K[i]==1])
        sim = (float(insect)/float(union))
        bdf_1M.loc[index,"10K_similarity"] = sim


jacc100 = bdf_1M["100K_similarity"].to_numpy().flatten()
mean_cov = bdf_1M["mean_cov"].to_numpy().flatten()
jacc10 = bdf_1M["10K_similarity"].to_numpy().flatten()

len100 = len(bdf_1M[bdf_1M["100K_similarity"] < 1]["Pangenome"].to_list())
len10 = len(bdf_1M[bdf_1M["10K_similarity"] < 1]["Pangenome"].to_list())
print("Cov threshold is:", cov_threshold)
print("The number of motus that differ from the 1M coverage estimation: ")
print("100K to 1M: ", len100)
print("10K to 1M: ", len10)

plt.figure(dpi=1200)
plt.scatter(jacc100, mean_cov, s = 2, linewidths=0)
plt.title("1M to 100K jaccard similarity") 
plt.xlabel("Jaccard similarity") 
plt.ylabel("Mean coverage")
plt.show()

plt.figure(dpi=1200)
plt.scatter(jacc10, mean_cov, s = 2, linewidths=0)
plt.title("1M to 10K jaccard similarity") 
plt.xlabel("Jaccard similarity") 
plt.ylabel("Mean coverage")
plt.show()

#%%
bdf_1M.to_csv("cov20_comparison1M.csv", sep="\t")

#%%manually checking some
genomes = bdf_1M[bdf_1M["10K_similarity"] < 1]["Pangenome"].to_list()
len(bdf_1M[bdf_1M["100K_similarity"] < 1]["Pangenome"].to_list())

sub_1M = bdf_1M[bdf_1M["Pangenome"].isin(genomes)]
sub_100K = bdf_100K[bdf_100K["Pangenome"].isin(genomes)]
sub_10K = bdf_10K[bdf_10K["Pangenome"].isin(genomes)]

#from these I should be able to compare the passed columns to see which ones have more/less?
#if it's even relevant.
#%% this is probably the wrong way since he wants it for each mOTU, not overall. But I could get an overall too right?
def binarize20(num):
    if num >= 20:
        return 1
    else:
        return 0
    
binarize = np.vectorize(binarize20)
boneM = binarize(oneM)
bhundredK = binarize(hundredK)
btenK = binarize(tenK)


#%%
A = {1,1,0,1}
B = {1,1,1,1,1}
# Intersection and Union of two sets can also be done using & and | operators.
C = A.intersection(B)
D = A.union(B)
print('AnB = ', C)
print('AUB = ', D)
print('J(A,B) = ', float(len(C))/float(len(D)))

#%%
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)

ax1.plot(df_1M['Loc090402-8m'])
ax1.plot(df_100K['Loc090402-8m'])

df_1M.plot(subplots = True)

group = cov.groupby(df_1M["Pangenome"])

pd.plotting.scatter_matrix(df_1M)

cov.set_index("Pangenome").plot()
cov.set_index("nr_subsampling").plot()
#%% Fernandos suggestion:
#present1M = {sample for i,sample in enumerate(bdf_1M.columns) if bdf_1M[pan,i] if "Loc" in sample}

