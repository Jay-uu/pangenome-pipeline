/*
Takes the output from mOTUlizer and the bins and sorts them into new directories based on which mOTU the belong to.
Input: Tuple with: group, name of the taxonomic classification/prefix for filenames. motus_file, has information about which bins belong to which
mOTU, and the bintable with bins quality data. Also takes the bins as input.
Output is a tuple with the new mOTU directory and the bintable.
*/
process create_mOTU_dirs {
    label "short_time"
    publishDir "${params.project}/mOTUs", mode: "copy", pattern: "${group}_mOTU_*"
    input:
    tuple(val(group), path(motus_file), path(bintable))
    path(bins)
    output:
    tuple(path("${group}_mOTU_*", type: "dir"), path("${bintable}"), optional: true)
    shell:
    '''
    #!/usr/bin/env python
    import os
    import shutil
    import pandas as pd
    import glob
    
    min_genomes = !{params.min_mOTU_MAGs} #nextflow param
    max_contam = 0.5 #change to nf param
    with open("!{motus_file}") as infile:
        for line in open("!{motus_file}"):
            if line.startswith("mOTU_"):
                fields = line.strip().split("\t")
                mOTU, rep, mean_ANI, min_ANI, missing_edges, MAGs, *SUBs = fields
                MAGs = MAGs.split(";")
                if len(MAGs) >= min_genomes:
                    bintdf = pd.read_csv("!{bintable}", sep = '\t')
                    if len(SUBs) > 0:
                        SUBs = SUBs[0].split(";")
                        MAGs.extend(SUBs)
                    #check bintdf if any MAG/SUB passes the contam filter, otherwise this mOTU will not be used.
                    if (bintdf["Contamination"] < max_contam).any():
                        os.mkdir("!{group}_" + mOTU)
                        for genome in MAGs:
                            #move file to mOTU directory
                            #Only bins with up to the maximum contamination threshold should be used. Any completeness is fine.
                            contam = bintdf.loc[bintdf["Bin Id"]==genome]["Contamination"].item()
                            if contam <= max_contam:
                                shutil.copy2(genome + ".fa", "!{group}_" + mOTU + "/")
                                print(f"{genome} being copied into !{group}_{mOTU} directory")
    '''
}
