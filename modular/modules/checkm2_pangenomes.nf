/*
This process might be removed if we decide to use checkM instead
*/
process checkm2_pangenomes {
    publishDir "${params.project}/mOTUs/results/${pangenome_dir.baseName}/pangenome", mode: "copy", pattern: "${pangenome_dir.baseName}_cM2/quality_report.tsv", saveAs: { "${pangenome_dir.baseName}_cM2_summary.txt" }
    label "checkm2_pangenomes"
    label "high_mem"
    tag "${pangenome_dir.baseName}"
    input:
    path(pangenome_dir)
    output:
    path("${pangenome_dir.baseName}_cM2/quality_report.tsv", type: "path")
    script:
    """
    #!/usr/bin/env bash
    pang_id="${pangenome_dir.baseName}"
    echo "Running checkM2 on \${pang_id} fastas."
    checkm2 predict --threads ${task.cpus} --input ${pangenome_dir} -x fasta --output-directory \${pang_id}_cM2
    """
}
