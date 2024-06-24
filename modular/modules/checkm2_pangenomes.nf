/*
This process might be removed if we decide to use checkM instead
*/
process checkm2_pangenomes {
    publishDir "${params.project}/mOTUs/results/${pangenome_dir.baseName}/pangenome", mode: "copy", pattern: "${pangenome_dir.baseName}_cM2/quality_report.tsv", saveAs: { filename -> "${pangenome_dir.baseName}_cM2_summary.txt" }
    publishDir "${params.project}/checkm_pangenomes", mode: "copy" //remove when finalizing pipeline
    //temporary until checkm2 in SQM env
    conda '/home/jay/mambaforge/envs/checkm2'
    //conda '/crex/proj/fume/nobackup/private/jay/mamba_envs/checkm2'
    label "no_label"
    tag "${pangenome_dir.baseName}"
    input:
    path(pangenome_dir)
    output:
    path("${pangenome_dir.baseName}_cM2/quality_report.tsv", type: "path")
    shell:
    '''
    pang_id=!{pangenome_dir.baseName}
    echo "Running checkM2 on !{pangenome_dir.baseName} fastas."
    checkm2 predict --threads !{task.cpus} --input !{pangenome_dir} -x fasta --output-directory ${pang_id}_cM2
    '''
}
