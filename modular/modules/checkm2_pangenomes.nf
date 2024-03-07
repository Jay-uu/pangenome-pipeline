/*
This process might be removed if we decide to use checkM instead
*/
process checkm2_pangenomes {
    publishDir "${params.project}/mOTUs/${pangenome_dir.baseName}/checkm", mode: "copy", pattern: "*_cM2/quality_report.tsv", saveAs: { filename -> "${pangenome_dir.baseName}_cM2_summary.txt" }

    publishDir "${params.project}/checkm_pangenomes", mode: "copy" //remove when finalizing pipeline
    //temporary until we decide which checkM version to use
    //conda '/home/jay/mambaforge/envs/checkm2'
    conda '/crex/proj/fume/nobackup/private/jay/mamba_envs/checkm2'
    tag "no_label"
    input:
    path(pangenome_dir)
    output:
    path("*_cM2", type: "dir")
    shell:
    '''
    pang_id=!{pangenome_dir.baseName}
    
    echo "Running checkM2 on !{pangenome_dir.baseName} fastas."
    checkm2 predict --threads !{params.threads} --input !{pangenome_dir} -x fasta --output-directory ${pang_id}_cM2

    '''
}
