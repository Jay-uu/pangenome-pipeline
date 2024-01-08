/*
This process might be removed if we decide to use checkM instead
*/
process checkm2_pangenomes {
    publishDir "${params.project}/checkm_pangenomes", mode: "copy"
    //temporary until we decide which checkM version to use
    //conda '/home/jay/mambaforge/envs/checkm2'
    conda '/crex/proj/fume/nobackup/private/jay/mamba_envs/checkm2'
    input:
    path(pangenome_dir)
    output:
    path("*_cM2", type: "dir")
    shell:
    '''
    pang_id=!{pangenome_dir.baseName}
    
    echo "Running checkM2 on core genome."
    checkm2 predict --threads !{params.threads} --input !{pangenome_dir} -x fasta --output-directory ${pang_id}_core_cM2 --output-directory ${pang_id}_cM2

    '''
}
