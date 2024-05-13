#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:34:25 2023

@author: jay
"""

#combine checkm results
import glob
import pandas as pd
#path = "/crex/proj/fume/nobackup/private/jay/test_pipeline/loclat_20240110/checkm_pangenomes"
#path = "/crex/proj/fume/nobackup/private/jay/test_pipeline/loclat_20240110/checkm_pangenomes"
#path = "/crex/proj/fume/nobackup/private/jay/test_pipeline/loclat_20240202_all_binners_w_superpangcheckm/checkm_pangenomes"

path = "/crex/proj/fume/nobackup/private/jay/test_pipeline/loclat_202400426_100K_all_binners/checkm_pangenomes"
projname = "loclat_202400426_100K_all_binners"

#path = "/domus/h1/jay/squeezemeta/pangenome-pipeline/results/loclat_202400426_100K_metabat/checkm_pangenomes"
#projname = "loclat_202400426_100K_metabat"

#checkm 1 files *_cM1_summary.txt
#checkm2 *cM2/quality_report.tsv

df = pd.DataFrame()
for cm2 in glob.glob(path+"/*/quality_report.tsv"):
    print(cm2)
    name = cm2.split("/")[-2]
    df = pd.concat([pd.read_csv(cm2, sep='\t'),df], ignore_index=True)

cm2 = df[["Name","Completeness","Contamination"]]
cm2 = cm2.rename(columns={"Name": "Bin Id", "Completeness": "cM2_Completeness", "Contamination": "cM2_Contamination"})


header = ["Bin Id", "Marker", "Lineage", "# genomes", "# markers", "# marker sets",
          "0", "1", "2", "3", "4", "5+", "cM1_Completeness", "cM1_Contamination", "Strain heterogeneity"]
cm1 = pd.DataFrame(columns=header)
for f in glob.glob(path+"/*_cM1_summary.txt"):
    print(f)
    with open(f) as file:
        for line in file:
            if "c__" in line:
                print(line)
                fields = line.split()
                cm1.loc[len(cm1)] = fields

cm1 = cm1[["Bin Id", "cM1_Completeness", "cM1_Contamination"]]

out = pd.merge(cm1, cm2, on='Bin Id')

out.to_csv(path+ "/" + projname +"_combined_cm1_cm2_summary.tsv", sep = '\t', index=False)
