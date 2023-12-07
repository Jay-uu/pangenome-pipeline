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
    tuple(path("${group}_mOTU_*", type: "dir"), path("${bintable}"))
    shell:
    '''
    #!/usr/bin/env python
    import os
    import shutil
    min_genomes = !{params.min_mOTU_MAGs} #nextflow param
    with open("!{motus_file}") as infile:
        for line in open("!{motus_file}"):
            if line.startswith("mOTU_"):
                fields = line.strip().split("\t")
                mOTU, rep, mean_ANI, min_ANI, missing_edges, nb_MAGs, nb_SUBS, MAGs, *SUBs = fields
                MAGs = MAGs.split(";")
                if int(nb_MAGs) >= min_genomes:
                    os.mkdir("!{group}_" + mOTU)
                    #will add here to make subs optional
                    if int(nb_SUBS) > 0:
                        SUBs = SUBs[0].split(";")
                        MAGs.extend(SUBs)
                    for genome in MAGs:
                        #move file to mOTU directory
                        shutil.copy2(genome + ".fa", "!{group}_" + mOTU + "/")
                        print(f"{genome} being copied into !{group}_{mOTU} directory")
    '''
}
