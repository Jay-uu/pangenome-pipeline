/*
Clusters bins based on similarity.
Input is a directory containing bins and a file with bin-name, completeness and contamination.
Output is the name of the taxonomic classification of the bins (unless taxSort = root),
a tsv with which bins belong to which mOTU, and the bintable file with quality data for the bins.
*/
process bins_to_mOTUs {
    //This process sometimes fails on UPPMAX due to lack of memory. Might implement something like this: https://www.nextflow.io/docs/latest/process.html#dynamic-computing-resources
    publishDir "${params.project}/mOTUs/simi_files", mode: "copy", pattern: "*_similarities.txt"
    label "low_cpu"
    tag "low_cpu"
    input:
    path(tax_dir)
    output:
    tuple(env(group), path("*_mOTUs.tsv"), path("${tax_dir}/*.bintable"), emit: mOTUs_tuple)
    path("*_similarities.txt", emit: simi_file)
    shell:
    '''
    #!/bin/bash
    group="!{tax_dir}"
    group=${group%"_bins"}
    echo $group
    mOTUlize.py --fnas !{tax_dir}/*.fa --checkm !{tax_dir}/*.bintable --MAG-completeness !{params.MAGcomplete} --MAG-contamination !{params.MAGcontam} --threads !{task.cpus} --keep-simi-file ${group}_similarities.txt -o ${group}_mOTUs.tsv
    '''
}
