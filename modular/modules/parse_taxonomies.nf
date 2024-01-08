/*
Reads the GTDB-Tk results from the bintable files then sorts bins based on the group selected with params.taxSort.
This can for example lead to species-directories, with all bins that were identified a specific species within the same directory.
Input is all bins, and all bintables.
Output is the new directories which contains bins and a tsv with completeness and contamination.
*/
process parse_taxonomies {
    label "low_cpu"
    input:
    path(bins)
    path(all_bintables)
    output:
    path("*_bins", emit: tax_bin_dirs)
    shell:
    $/
    #!/usr/bin/env python
    import os
    import glob
    import pandas as pd
    import shutil
    ranks = ["root","auto","domain", "phylum","class","order", "family","genus","species"]
    rank_param = "!{params.taxSort}"
    rank_param = rank_param.lower() #allow the user to spell with big letters
    if rank_param not in ranks:
        raise Exception(f"The provided sorting rank is not valid. Please use a taxonomic rank from the following list {ranks}")
    #then read files to df
    print("Reading bintables to dataframe")
    all_bintable_files = glob.glob("*.bintable")
    all_bintables = pd.DataFrame()
    for file in all_bintable_files:
        #read next file to be added, remove unneeded columns to save memory
        next_table = pd.read_csv(file, sep='\t', header=1)
        c_list = next_table.columns.to_list()
        c_list = [e for e in c_list if e not in 
                  ("Bin ID", "Tax GTDB-Tk", "Completeness", "Contamination")]
        next_table = next_table.drop(labels=c_list, axis=1)
        all_bintables = pd.concat((next_table, all_bintables), ignore_index=True)
    
    #rename bin id an taxonomy cols
    print("Renaming columns")
    all_bintables.rename(columns={"Bin ID": "Bin Id"}, inplace=True)
    all_bintables.rename(columns={"Tax GTDB-Tk": "root"}, inplace=True)
    #Convert to categorical data to save memory
    print("Convert taxonomy column to categorical")
    all_bintables["root"] = all_bintables["root"].astype("category")
    #too small genomes won't give any completeness and contamination data, which is necessary for mOTUlizer
    
    all_bintables.fillna(value={"Completeness" : 1, "Contamination" : 1}, inplace = True)
    #split Tax GTDB-Tk/root column so I can easily search ranks.
    print("Splitting taxonomy, and fixing Unclassified data.")
    all_bintables[ranks[2:]] = all_bintables["root"].str.split(";", expand=True)
    #making sure each rank column has either a classification or is unclassified
    rank_cols = (dict([[e,"Unclassified"] for e in ranks[2:]]))
    all_bintables.fillna(value=rank_cols, inplace=True)
    #it's also possible for bins that were only classified at lower levels
    #to have example s__ instead of a species. This should be unclassified instead.
    tax_short = ["d__", "p__","c__","o__","f__","g__","s__"]
    all_bintables.replace(to_replace=tax_short, value="Unclassified", inplace=True)
    if rank_param == "auto":
        #find lowest rank with more than 90% classified bins and change rank_param to that
        print("Selected rank: auto. Finding lowest well classified rank.")
        rank_not_selected = True
        i = 8
        tot_bins = len(all_bintables)
        while rank_not_selected:
            print(f"Checking unclassified proportion in {ranks[i]}.")
            if i == 2:
                rank_param = "root" #no lower classification was good enough quality
                rank_not_selected = False
            elif all_bintables[ranks[i]].value_counts()["Unclassified"]/tot_bins < 0.1:
                rank_param = ranks[i]
                rank_not_selected = False
            else:
                i = i-1 
    #start sorting based on rank_param
    #create dirs and taxonomic-sorted bintables 
    print(f"Sorting based on {rank_param}.")
    if rank_param == "root":
        all_bintables["root"] = "root"
    for group in all_bintables[rank_param].unique():
        print(f"Writing {group} to dirs")
        g = all_bintables[all_bintables[rank_param] == group]
        #check that the group has at least one bin of good enough quality for mOTUlizer. Otherwise skip.
        if len(g.loc[(g["Completeness"] > !{params.MAGcomplete}) & (g["Contamination"] < !{params.MAGcontam})]) == 0:
            print(f"No {group} bin is good enough quality for clustering. Exlcuding from further analysis.")
        else:
            os.makedirs(group+"_bins")
            g.to_csv(f"{group}_bins/{group}.bintable", index=None,
                     sep='\t', columns = ["Bin Id","Completeness","Contamination"])
            for bin_id in g['Bin Id']:
                shutil.copy2(bin_id+".fa", group+"_bins")
    /$
    
}
