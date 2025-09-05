/*
Clusters bins based on similarity using https://github.com/moritzbuck/mOTUlizer
Input is a directory containing bins and a file with bin-name, completeness and contamination.
Output is the name of the taxonomic classification of the bins (unless taxSort = root),
a tsv with which bins belong to which mOTU, and the bintable file with quality data for the bins.
*/
process bins_to_mOTUs {
    publishDir "${project_path}/mOTUs/mOTUlizer", mode: "copy", pattern: "*_similarities.txt", saveAs: { filename -> "${tax_dir}" - "_bins" + "/" + filename }
    publishDir "${project_path}/mOTUs/mOTUlizer", mode: "copy", pattern: "*_mOTUs.tsv", saveAs: { filename -> "${tax_dir}" - "_bins" + "/" + filename }
    label "low_cpu"
    label "bins_to_mOTUs"
    tag "${tax_dir.simpleName}"
    input:
    path(tax_dir)
    val(project_path)
    val(MAGcomplete)
    val(MAGcontam)
    output:
    tuple(env('group'), path("*_mOTUs.tsv"), path("${tax_dir}/*.bintable"), emit: mOTUs_tuple)
    path("*_similarities.txt", emit: simi_file)
    script:
    """
    #!/usr/bin/env bash
    group="${tax_dir}"
    group=\${group%"_bins"}
    echo \$group
    mOTUlize.py --fnas ${tax_dir}/*.fa --checkm ${tax_dir}/*.bintable --MAG-completeness ${MAGcomplete} --MAG-contamination ${MAGcontam} --threads ${task.cpus} --keep-simi-file \${group}_similarities.txt -o \${group}_mOTUs.tsv
    """
}
